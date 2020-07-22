#include <stdio.h>
#include <math.h>

//#define J0_SMF(omega) (((1 - S2) * tau) / (1 + pow(omega * tau, 2)))

//#define J0_SMF_num(omega) ((1 - S2) * tau)
//#define J0_SMF_den(omega) (1 + pow(omega * tau, 2))

long double J0_SMF(double omega, long double * tau, long double * S2) {
	return (((1 - (*S2)) * (*tau)) / (1 + powl((omega) * (*tau), 2)));
}

long double J0_EMF(double omega, long double * taus, long double * S2s, long double * tauf, long double * S2f) {
	long double fast, slow;
	fast = (((1 - (*S2f)) * (*tauf)) / (1 + powl((omega) * (*tauf), 2)));
	slow = (((*S2f) * (1 - (*S2s)) * (*taus)) / (1 + powl((omega) * (*taus), 2)));
	return fast + slow;
}

double SMF_R1(struct Residue *res, struct Relaxation* relax, long double tau, long double S2, int mode) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */

	long double field = relax->field * 1000000; // conversion to Hz
	
	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a float128, and
	 * software division) I've predefined it. Constants are in datatypes.c */

	long double omega_1H = 2 * M_PI * field;
	long double omega_L;
	long double d, d2x, d2y, d2xy, d2tot, d2;
	float *csa;

	if (mode == MODE_15N) {
		csa = res->csaN;
		omega_L = 2 * M_PI * field / 9.869683408806043;
		d = -D_NH;
	} else if (mode == MODE_13C) {
		csa = res->csaC;
		omega_L = 2 * M_PI * field / 3.976489314034722;
		d = -D_CH;
	}

	d2x = powl(((csa[2] - csa[0]) * powl(10.,-6.)) * omega_L, 2.);
	d2y = powl(((csa[1] - csa[0]) * powl(10.,-6.)) * omega_L, 2.);
	d2xy= powl(powl(10.,-6.) * omega_L, 2.) * (csa[2] - csa[0]) * (csa[1] - csa[0]);
	d2tot= (powl(csa[2], 2) + powl(csa[1], 2) + powl(csa[0], 2));
	d2tot+=(-(csa[2] * csa[1] + csa[1] * csa[0] + csa[2] * csa[0]));
	d2tot*=powl(powl(10., -6.) * omega_L, 2);
	
	double R1D = 0, R1CSA = 0;
	long double Jcomp = 0;
	Jcomp += J0_SMF(omega_1H - omega_L, &tau, &S2);
	Jcomp += 3 * J0_SMF(omega_L, &tau, &S2);
	Jcomp += 6 * J0_SMF(omega_1H + omega_L, &tau, &S2);
	R1D = 0.1 * d * d * Jcomp;
	R1CSA = (2/15.) * d2tot * J0_SMF(omega_L, &tau, &S2);

	return R1D + R1CSA;
}

double SMF_R2(struct Residue *res, struct Relaxation* relax, long double tau, long double S2, int mode) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */

	long double field = relax->field * 1000000; // conversion to Hz

	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a float128, and
	 * software division) I've predefined it. Bond length taken as 1.02 A */

	//long double d = -72084.44597;
	long double omega_1H = 2 * M_PI * field;
	//long double omega_15N = 2 * M_PI * field / 9.869683408806043;
	//long double d2=(-170)*(pow(10, -6))*(omega_15N);

	long double omega_L;
	long double d, d2x, d2y, d2xy, d2tot, d2;
	float *csa;

	if (mode == MODE_15N) {
		csa = res->csaN;
		omega_L = 2 * M_PI * field / 9.869683408806043;
		d = -D_NH;
	} else if (mode == MODE_13C) {
		csa = res->csaC;
		omega_L = 2 * M_PI * field / 3.976489314034722;
		d = -D_CH;
	}

	d2x = powl(((csa[2] - csa[0]) * powl(10.,-6.)) * omega_L, 2.);
	d2y = powl(((csa[1] - csa[0]) * powl(10.,-6.)) * omega_L, 2.);
	d2xy= powl(powl(10.,-6.) * omega_L, 2.) * (csa[2] - csa[0]) * (csa[1] - csa[0]);
	d2tot= (powl(csa[2], 2) + powl(csa[1], 2) + powl(csa[0], 2));
	d2tot+=(-(csa[2] * csa[1] + csa[1] * csa[0] + csa[2] * csa[0]));
	d2tot*=powl(powl(10., -6.) * omega_L, 2);
	
	long double w1 = relax->w1;
	long double wr = relax->wr;

	long double R2NH = 0, R2NCSA = 0;
	long double J0sum = 0;
	J0sum += (2/3.) * J0_SMF(2 * M_PI * (w1 - 2 * wr), &tau, &S2);
	J0sum += (2/3.) * J0_SMF(2 * M_PI * (w1 + 2 * wr), &tau, &S2);
	J0sum += (4/3.) * J0_SMF(2 * M_PI * (w1 - wr), &tau, &S2);
	J0sum += (4/3.) * J0_SMF(2 * M_PI * (w1 + wr), &tau, &S2);

	long double JNH = J0sum + 3 * J0_SMF(omega_L, &tau, &S2);
	JNH += J0_SMF(omega_1H - omega_L, &tau, &S2);
	JNH += 6 * J0_SMF(omega_1H, &tau, &S2);
	JNH += 6 * J0_SMF(omega_1H + omega_L, &tau, &S2);


	R2NH = (1/20.) * d * d * JNH;
	R2NCSA = (1/45.) * d2tot * (J0sum + 3 * J0_SMF(omega_L, &tau, &S2));
	return R2NH + R2NCSA;
}

double EMF_15NR1(struct Residue *res, struct Relaxation* relax, long double taus, long double S2s, long double tauf) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	long double S2f = res->S2_dipolar / S2s;
	long double field = relax->field * 1000000; // conversion to Hz

	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a float128, and
	 * software division) I've predefined it. Bond length taken as 1.02 A */

	long double d = -72084.44597;
	long double omega_1H = 2 * M_PI * field;
	long double omega_15N = 2 * M_PI * field / 9.869683408806043;
	long double d2=(-170)*(pow(10, -6))*(omega_15N);

	double R1NH = 0, R1NCSA = 0;
	long double Jcomp = 0;

	// long double * taus, long double * S2s, long double * tauf, long double * S2f)
	Jcomp += J0_EMF(omega_1H - omega_15N, &taus, &S2s, &tauf, &S2f);
	Jcomp += 3 * J0_EMF(omega_15N, &taus, &S2s, &tauf, &S2f);
	Jcomp += 6 * J0_EMF(omega_1H + omega_15N, &taus, &S2s, &tauf, &S2f);
	R1NH = 0.1 * d * d * Jcomp;
	R1NCSA = (2/15.) * d2 * d2 * J0_EMF(omega_15N, &taus, &S2s, &tauf, &S2f);


	return R1NH + R1NCSA;
}

double EMF_15NR2(struct Residue *res, struct Relaxation* relax, long double taus, long double S2s, long double tauf) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	long double S2f = res->S2_dipolar / S2s;
	long double field = relax->field * 1000000; // conversion to Hz
	//printf("Field: %0.1f MHz\nS2f: %Le\nS2s: %Le\nts: %Le\ntf: %Le\nS2d: %f\ncsiso: %f\n", relax->field, S2f, S2s, taus, tauf, res->S2_dipolar, res->csisoN);
	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a float128, and
	 * software division) I've predefined it. Bond length taken as 1.02 A */

	long double d = -72084.44597;
	long double omega_1H = 2 * M_PI * field;
	long double omega_15N = 2 * M_PI * field / 9.869683408806043;
	long double d2=(-170)*(pow(10, -6))*(omega_15N);

	long double w1 = relax->w1;
	long double wr = relax->wr;

	long double R2NH = 0, R2NCSA = 0;
	long double J0sum = 0;
	J0sum += (2/3.) * J0_EMF(2 * M_PI * (w1 - 2 * wr), &taus, &S2s, &tauf, &S2f);
	J0sum += (2/3.) * J0_EMF(2 * M_PI * (w1 + 2 * wr), &taus, &S2s, &tauf, &S2f);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 - wr), &taus, &S2s, &tauf, &S2f);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 + wr), &taus, &S2s, &tauf, &S2f);

	long double JNH = J0sum + 3 * J0_EMF(omega_15N, &taus, &S2s, &tauf, &S2f);
	JNH += J0_EMF(omega_1H - omega_15N, &taus, &S2s, &tauf, &S2f);
	JNH += 6 * J0_EMF(omega_1H, &taus, &S2s, &tauf, &S2f);
	JNH += 6 * J0_EMF(omega_1H + omega_15N, &taus, &S2s, &tauf, &S2f);


	R2NH = (1/20.) * d * d * JNH;
	R2NCSA = (1/45.) * d2 * d2 * (J0sum + 3 * J0_EMF(omega_15N, &taus, &S2s, &tauf, &S2f));
	return R2NH + R2NCSA;
}
