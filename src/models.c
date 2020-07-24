#include <stdio.h>
#include <math.h>
#include <complex.h>

#define J0_SMF(omega, tau, S2) ((((1 - (double) S2) * (long double) tau)) / (1 + ((double) omega * (long double) tau * (double) omega * (long double) tau)))
#define J0_EMF(omega, taus, S2s, tauf, S2f) (((((1 - (double) S2f)) * (long double) tauf)) / (1 + ((double) omega * (long double) tauf * (double) omega * (long double) tauf))) + ((((double) S2f) * (1 - (double) S2s) * (long double) taus) / (1 + ((double) omega * (long double) taus * (double) omega * (long double) taus)))

#define sq_i(x) (int) x * (int) x
#define sq(x) (double) x * (double) x
/*long double J0_SMF(double omega, long double * tau, long double * S2) {
	return (((1 - (*S2)) * (*tau)) / (1 + powl((omega) * (*tau), 2)));
}

long double J0_EMF(double omega, long double * taus, long double * S2s, long double * tauf, long double * S2f) {
	//long double fast, slow;
	//fast = (((1 - (*S2f)) * (*tauf)) / (1 + powl((omega) * (*tauf), 2)));
	//slow = (((*S2f) * (1 - (*S2s)) * (*taus)) / (1 + powl((omega) * (*taus), 2)));
	return (((1 - (*S2f)) * (*tauf)) / (1 + powl((omega) * (*tauf), 2))) + (((*S2f) * (1 - (*S2s)) * (*taus)) / (1 + powl((omega) * (*taus), 2)));
}*/

double SMF_R1(struct Residue *res, struct Relaxation* relax, long double tau, long double S2, int mode) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */

	long double field = relax->field * 1000000; // conversion to Hz

	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a float128, and
	 * software division) I've predefined it. Constants are in datatypes.c */

	long double omega_1H = 2 * M_PI * field;
	long double omega_L;
	long double d, d2tot;
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

	//d2x = powl(((csa[2] - csa[0]) * powl(10.,-6.)) * omega_L, 2.);
	//d2y = powl(((csa[1] - csa[0]) * powl(10.,-6.)) * omega_L, 2.);
	//d2xy= powl(powl(10.,-6.) * omega_L, 2.) * (csa[2] - csa[0]) * (csa[1] - csa[0]);
	d2tot= (sq(csa[2]) + sq(csa[1]) + sq(csa[0]));
	d2tot+=(-(csa[2] * csa[1] + csa[1] * csa[0] + csa[2] * csa[0]));
	d2tot*=sq(0.000001 * omega_L);

	double R1D = 0, R1CSA = 0;
	long double Jcomp = 0;
	Jcomp += J0_SMF(omega_1H - omega_L, tau, S2);
	Jcomp += 3 * J0_SMF(omega_L, tau, S2);
	Jcomp += 6 * J0_SMF(omega_1H + omega_L, tau, S2);
	R1D = 0.1 * d * d * Jcomp;
	R1CSA = (2/15.) * d2tot * J0_SMF(omega_L, tau, S2);

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
	long double d, d2tot;
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

	//d2x = powl(((csa[2] - csa[0]) * powl(10.,-6.)) * omega_L, 2.);
	//d2y = powl(((csa[1] - csa[0]) * powl(10.,-6.)) * omega_L, 2.);
	//d2xy= powl(powl(10.,-6.) * omega_L, 2.) * (csa[2] - csa[0]) * (csa[1] - csa[0]);
	d2tot= (sq(csa[2]) + sq(csa[1]) + sq(csa[0]));
	d2tot+=(-(csa[2] * csa[1] + csa[1] * csa[0] + csa[2] * csa[0]));
	d2tot*=sq(0.000001 * omega_L);

	long double w1 = relax->w1;
	long double wr = relax->wr;

	long double R2NH = 0, R2NCSA = 0;
	long double J0sum = 0;
	J0sum += (2/3.) * J0_SMF(2 * M_PI * (w1 - 2 * wr), tau, S2);
	J0sum += (2/3.) * J0_SMF(2 * M_PI * (w1 + 2 * wr), tau, S2);
	J0sum += (4/3.) * J0_SMF(2 * M_PI * (w1 - wr), tau, S2);
	J0sum += (4/3.) * J0_SMF(2 * M_PI * (w1 + wr), tau, S2);

	long double JNH = J0sum + 3 * J0_SMF(omega_L, tau, S2);
	JNH += J0_SMF(omega_1H - omega_L, tau, S2);
	JNH += 6 * J0_SMF(omega_1H, tau, S2);
	JNH += 6 * J0_SMF(omega_1H + omega_L, tau, S2);


	R2NH = (1/20.) * d * d * JNH;
	R2NCSA = (1/45.) * d2tot * (J0sum + 3 * J0_SMF(omega_L, tau, S2));
	return R2NH + R2NCSA;
}

double EMF_R1(struct Residue *res, struct Relaxation* relax, long double taus, long double S2s, long double tauf, int mode) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	long double S2f = res->S2_dipolar / S2s;
	long double field = relax->field * 1000000; // conversion to Hz

	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a float128, and
	 * software division) I've predefined it. Bond length taken as 1.02 A */

	long double omega_1H = 2 * M_PI * field;
	long double omega_L;
	long double d, d2tot;
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

	//d2x = powl(((csa[2] - csa[0]) * powl(10.,-6.)) * omega_L, 2.);
	//d2y = powl(((csa[1] - csa[0]) * powl(10.,-6.)) * omega_L, 2.);
	//d2xy= powl(powl(10.,-6.) * omega_L, 2.) * (csa[2] - csa[0]) * (csa[1] - csa[0]);
	d2tot= (sq(csa[2]) + sq(csa[1]) + sq(csa[0]));
	d2tot+=(-(csa[2] * csa[1] + csa[1] * csa[0] + csa[2] * csa[0]));
	d2tot*=sq(0.000001 * omega_L);

	double R1D = 0, R1CSA = 0;
	long double Jcomp = 0;

	// long double * taus, long double * S2s, long double * tauf, long double * S2f)
	Jcomp += J0_EMF(omega_1H - omega_L, taus, S2s, tauf, S2f);
	Jcomp += 3 * J0_EMF(omega_L, taus, S2s, tauf, S2f);
	Jcomp += 6 * J0_EMF(omega_1H + omega_L, taus, S2s, tauf, S2f);
	R1D = 0.1 * d * d * Jcomp;
	R1CSA = (2/15.) * d2tot * J0_EMF(omega_L, taus, S2s, tauf, S2f);


	return R1D + R1CSA;
}

double EMF_R2(struct Residue *res, struct Relaxation* relax, long double taus, long double S2s, long double tauf, int mode) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	long double S2f = res->S2_dipolar / S2s;
	long double field = relax->field * 1000000; // conversion to Hz
	//printf("Field: %0.1f MHz\nS2f: %Le\nS2s: %Le\nts: %Le\ntf: %Le\nS2d: %f\ncsiso: %f\n", relax->field, S2f, S2s, taus, tauf, res->S2_dipolar, res->csisoN);
	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a float128, and
	 * software division) I've predefined it. Bond length taken as 1.02 A */

	long double omega_1H = 2 * M_PI * field;
	long double omega_L;
	long double d, d2tot;
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

	//d2x = powl(((csa[2] - csa[0]) * powl(10.,-6.)) * omega_L, 2.);
	//d2y = powl(((csa[1] - csa[0]) * powl(10.,-6.)) * omega_L, 2.);
	//d2xy= powl(powl(10.,-6.) * omega_L, 2.) * (csa[2] - csa[0]) * (csa[1] - csa[0]);
	d2tot= (sq(csa[2]) + sq(csa[1]) + sq(csa[0]));
	d2tot+=(-(csa[2] * csa[1] + csa[1] * csa[0] + csa[2] * csa[0]));
	d2tot*=sq(0.000001 * omega_L);

	long double w1 = relax->w1;
	long double wr = relax->wr;

	long double R2D = 0, R2CSA = 0;
	long double J0sum = 0;
	J0sum += (2/3.) * J0_EMF(2 * M_PI * (w1 - 2 * wr), taus, S2s, tauf, S2f);
	J0sum += (2/3.) * J0_EMF(2 * M_PI * (w1 + 2 * wr), taus, S2s, tauf, S2f);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 - wr), taus, S2s, tauf, S2f);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 + wr), taus, S2s, tauf, S2f);

	long double JNH = J0sum + 3 * J0_EMF(omega_L, taus, S2s, tauf, S2f);
	JNH += J0_EMF(omega_1H - omega_L, taus, S2s, tauf, S2f);
	JNH += 6 * J0_EMF(omega_1H, taus, S2s, tauf, S2f);
	JNH += 6 * J0_EMF(omega_1H + omega_L, taus, S2s, tauf, S2f);


	R2D = (1/20.) * d * d * JNH;
	R2CSA = (1/45.) * d2tot * (J0sum + 3 * J0_EMF(omega_L, taus, S2s, tauf, S2f));
	return R2D + R2CSA;
}


long double GAF_S2(long double sig[3], struct Orient * A, struct Orient * B, int mode) {
	/* #define MODE_REAL		0
#define MODE_IMAG		1
#define MODE_COMP		2*/

	int l, m, mp, k, kp;
	double complex Amp;
	double complex temp;
	double expo = 0;

	int ksqsum;
	int lsqsum;
	int msqsum;
	//l = 2;k=0;kp=2; m=0;mp=0;
	for (l = -2; l <= 2; l++) {
		lsqsum = sq_i(l);
		for (m = -2; m <= 2; m++) {
			for (mp = -2; mp <= 2; mp++) {
				msqsum = sq_i(m) + sq_i(mp);
				for (k = -2; k <= 2; k++) {
					for (kp = -2; kp <= 2; kp++) {
						ksqsum = sq_i(k) + sq_i(kp);

						if (A->Y2[m+2] == A->Y2c[m+2] && B->Y2[mp+2] == B->Y2c[mp+2]) {
							/* in this case then Y2 and Y2c are real numbers and therefore the last step of this is real.
							 * If mode is MODE_REAL, then, we may safely ignore any case in which
							 *  cpowl(-1*I, k-kp)
							 * is imaginary (as this will give only an imaginary component to the Amp).
							 * This is the case for any (k-kp)%2 != 0.
							 * In this case, if (k-kp)%4 == 0 the cpowl function gives 1, else -1.
							 */
							if ((k-kp)%2 != 0)
								continue;
							if ((k - kp) % 4 == 0)
								temp = 1;
							else
								temp = -1;
						} else {
							/* Other wise we have to go the long route... */
							temp = cpowl(-1 * I, k - kp);
						}
						//temp = 0;
						//temp = cpowl(-1 * I, k - kp); // come back to this
						expo = 0;
						expo = -(sq(sig[0]) * ksqsum / 2.);
						expo -= (sq(sig[1]) * lsqsum);
						expo -= (sq(sig[2]) * msqsum / 2.);
						temp *= expl(expo);
						temp *= Dwig[k+2][l+2] * Dwig[kp+2][l+2] * Dwig[m+2][k+2] * Dwig[mp+2][kp+2];
						temp *= A->Y2[m+2] * B->Y2c[mp+2];
						Amp += temp;
					}
				}
			}
		}
	}
	Amp *= (4 * M_PI / 5.);
	switch (mode) {
		case MODE_REAL: return creall(Amp); break;
		case MODE_IMAG: return cimagl(Amp); break;
		default: break;
	}
	return -1;
}

/* This has not been directly tested against the MATLAB model as I don't really have an equivalent to it set up (NHCO_3GAFSquaredEa is likely the closest but I don't have the setup around it for the Legendre polynomials) */
double GAF_15NR1(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double * sigf) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	long double field = relax->field * 1000000; // conversion to Hz

	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a float128, and
	 * software division) I've predefined it. Bond length taken as 1.02 A */

	long double omega_1H = 2 * M_PI * field;
	long double omega_13C, omega_15N, wCOCa;
	long double d2x, d2y, d2xy;


	omega_13C = 2 * M_PI * field / 3.976489314034722;
	omega_15N = 2 * M_PI * field / 9.869683408806043;
	wCOCa = 120 * omega_13C * 0.000001;

	float *csa;

	csa = res->csaN;



	d2x = sq(((csa[2] - csa[0]) * 0.000001) * omega_15N);
	d2y = sq(((csa[1] - csa[0]) * 0.000001) * omega_15N);
	d2xy= sq(0.000001 * omega_15N) * (csa[2] - csa[0]) * (csa[1] - csa[0]);


	/* Calculate all order parameters */
	long double S2NHs, S2NCSAxs, S2NCSAys, S2NCSAxys, S2CNs, S2CaNs;
	long double S2NHf, S2NCSAxf, S2NCSAyf, S2NCSAxyf, S2CNf, S2CaNf;
	S2NHs    = GAF_S2(sigs, &(res->orients[OR_NH]), &(res->orients[OR_NH]), MODE_REAL);
	S2NCSAxs = GAF_S2(sigs, &(res->orients[OR_NCSAxx]), &(res->orients[OR_NCSAxx]), MODE_REAL);
	S2NCSAys = GAF_S2(sigs, &(res->orients[OR_NCSAyy]), &(res->orients[OR_NCSAyy]), MODE_REAL);
	S2NCSAxys= GAF_S2(sigs, &(res->orients[OR_NCSAxx]), &(res->orients[OR_NCSAyy]), MODE_REAL);
	S2CNs    = GAF_S2(sigs, &(res->orients[OR_CN]), &(res->orients[OR_CN]), MODE_REAL);
	S2CaNs   = GAF_S2(sigs, &(res->orients[OR_NCA]), &(res->orients[OR_NCA]), MODE_REAL);

	S2NHf    = GAF_S2(sigf, &(res->orients[OR_NH]), &(res->orients[OR_NH]), MODE_REAL);
	S2NCSAxf = GAF_S2(sigf, &(res->orients[OR_NCSAxx]), &(res->orients[OR_NCSAxx]), MODE_REAL);
	S2NCSAyf = GAF_S2(sigf, &(res->orients[OR_NCSAyy]), &(res->orients[OR_NCSAyy]), MODE_REAL);
	S2NCSAxyf= GAF_S2(sigf, &(res->orients[OR_NCSAxx]), &(res->orients[OR_NCSAyy]), MODE_REAL);
	S2CNf    = GAF_S2(sigf, &(res->orients[OR_CN]), &(res->orients[OR_CN]), MODE_REAL);
	S2CaNf   = GAF_S2(sigf, &(res->orients[OR_NCA]), &(res->orients[OR_NCA]), MODE_REAL);

	/* N CSA relaxation contribution */
	//double R1D = 0, R1CSA = 0;

	double R1CSAx, R1CSAy, R1CSAxy, R1CSA, R1NH, R1NHr, R1CN, R1CaN;
	long double J1 = 0;

	// X
	J1 = J0_EMF(omega_15N, taus, S2NCSAxs, tauf, S2NCSAxf);
	R1CSAx = (2/15.) * d2x * J1;

	// Y
	J1 = J0_EMF(omega_15N, taus, S2NCSAys, tauf, S2NCSAyf);
	R1CSAy = (2/15.) * d2y * J1;

	// XY
	J1 = J0_EMF(omega_15N, taus, S2NCSAxys, tauf, S2NCSAxyf);
	R1CSAxy = (2/15.) * d2xy * J1;

	R1CSA = R1CSAx + R1CSAy - R1CSAxy;

	/* N Dipolar Interactions Contributions */

	long double J0HC, J2HC;

	// NH

	J0HC = J0_EMF(omega_1H - omega_15N, taus, S2NHs, tauf, S2NHf);
	J1   = J0_EMF(omega_15N, taus, S2NHs, tauf, S2NHf);
	J2HC = J0_EMF(omega_1H + omega_15N, taus, S2NHs, tauf, S2NHf);
	R1NH = (0.1) * sq(D_NH) * (J0HC + 3 * J1 + 6 * J2HC);

	// NHr - the original MATLAB uses the S2CN values

	J0HC = J0_EMF(omega_1H - omega_15N, taus, S2CNs, tauf, S2CNf);
	J1   = J0_EMF(omega_15N, taus, S2CNs, tauf, S2CNf);
	J2HC = J0_EMF(omega_1H + omega_15N, taus, S2CNs, tauf, S2CNf);
	R1NHr = (0.1) * sq(D_HNr) * (J0HC + 3 * J1 + 6 * J2HC);

	// CN

	J0HC = J0_EMF(omega_13C - omega_15N, taus, S2CNs, tauf, S2CNf);
	J1   = J0_EMF(omega_15N, taus, S2CNs, tauf, S2CNf);
	J2HC = J0_EMF(omega_13C + omega_15N, taus, S2CNs, tauf, S2CNf);
	R1CN = (0.1) * sq(D_CN) * (J0HC + 3 * J1 + 6 * J2HC);

	// CaN

	J0HC = J0_EMF(omega_13C - omega_15N, taus, S2CaNs, tauf, S2CaNf);
	J1   = J0_EMF(omega_15N, taus, S2CaNs, tauf, S2CaNf);
	J2HC = J0_EMF(omega_13C + omega_15N, taus, S2CaNs, tauf, S2CaNf);
	R1CaN = (0.1) * sq(D_CaN) * (J0HC + 3 * J1 + 6 * J2HC);

	return R1CSA + R1NH + R1NHr + R1CN + R1CaN;
}


/* ditto as GAF_15NR1 */
double GAF_15NR2(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double * sigf) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	double field = relax->field * 1000000; // conversion to Hz
	float w1 = relax->w1, wr = relax->wr;
	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a float128, and
	 * software division) I've predefined it. Bond length taken as 1.02 A */

	double omega_1H = 2 * M_PI * field;
	double omega_13C, omega_15N, wCOCa;
	double d2x, d2y, d2xy;

	omega_13C = 2 * M_PI * field / 3.976489314034722;
	omega_15N = 2 * M_PI * field / 9.869683408806043;
	wCOCa = 120 * omega_13C * 0.000001;

	float *csa;

	csa = res->csaN;



	d2x = sq(((csa[2] - csa[0]) * 0.000001) * omega_15N);
	d2y = sq(((csa[1] - csa[0]) * 0.000001) * omega_15N);
	d2xy= sq(0.000001 * omega_15N) * (csa[2] - csa[0]) * (csa[1] - csa[0]);


	/* Calculate all order parameters */
	double S2NHs, S2NCSAxs, S2NCSAys, S2NCSAxys, S2CNs, S2CaNs;
	double S2NHf, S2NCSAxf, S2NCSAyf, S2NCSAxyf, S2CNf, S2CaNf;
	S2NHs    = GAF_S2(sigs, &(res->orients[OR_NH]), &(res->orients[OR_NH]), MODE_REAL);
	S2NCSAxs = GAF_S2(sigs, &(res->orients[OR_NCSAxx]), &(res->orients[OR_NCSAxx]), MODE_REAL);
	S2NCSAys = GAF_S2(sigs, &(res->orients[OR_NCSAyy]), &(res->orients[OR_NCSAyy]), MODE_REAL);
	S2NCSAxys= GAF_S2(sigs, &(res->orients[OR_NCSAxx]), &(res->orients[OR_NCSAyy]), MODE_REAL);
	S2CNs    = GAF_S2(sigs, &(res->orients[OR_CN]), &(res->orients[OR_CN]), MODE_REAL);
	S2CaNs   = GAF_S2(sigs, &(res->orients[OR_NCA]), &(res->orients[OR_NCA]), MODE_REAL);

	S2NHf    = GAF_S2(sigf, &(res->orients[OR_NH]), &(res->orients[OR_NH]), MODE_REAL);
	S2NCSAxf = GAF_S2(sigf, &(res->orients[OR_NCSAxx]), &(res->orients[OR_NCSAxx]), MODE_REAL);
	S2NCSAyf = GAF_S2(sigf, &(res->orients[OR_NCSAyy]), &(res->orients[OR_NCSAyy]), MODE_REAL);
	S2NCSAxyf= GAF_S2(sigf, &(res->orients[OR_NCSAxx]), &(res->orients[OR_NCSAyy]), MODE_REAL);
	S2CNf    = GAF_S2(sigf, &(res->orients[OR_CN]), &(res->orients[OR_CN]), MODE_REAL);
	S2CaNf   = GAF_S2(sigf, &(res->orients[OR_NCA]), &(res->orients[OR_NCA]), MODE_REAL);

	/* N CSA relaxation contribution */
	//double R1D = 0, R1CSA = 0;

	double R2CSAx, R2CSAy, R2CSAxy, R2CSA, R2NH, R2NHr, R2CN, R2CaN;
	long double J1 = 0, J0sum;

	// X
	J0sum = 0; J1 = 0;
	J0sum += (2/3.) * J0_EMF(2 * M_PI * (w1 - 2 * wr), taus, S2NCSAxs, tauf, S2NCSAxf);
	J0sum += (2/3.) * J0_EMF(2 * M_PI * (w1 + 2 * wr), taus, S2NCSAxs, tauf, S2NCSAxf);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 - wr), taus, S2NCSAxs, tauf, S2NCSAxf);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 + wr), taus, S2NCSAxs, tauf, S2NCSAxf);
	J1 = J0_EMF(omega_15N, taus, S2NCSAxs, tauf, S2NCSAxf);
	R2CSAx = (1/45.) * d2x * (J0sum + 3 * J1);

	// Y
	J0sum = 0; J1 = 0;
	J0sum += (2/3.) * J0_EMF(2 * M_PI * (w1 - 2 * wr), taus, S2NCSAys, tauf, S2NCSAyf);
	J0sum += (2/3.) * J0_EMF(2 * M_PI * (w1 + 2 * wr), taus, S2NCSAys, tauf, S2NCSAyf);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 - wr), taus, S2NCSAys, tauf, S2NCSAyf);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 + wr), taus, S2NCSAys, tauf, S2NCSAyf);
	J1 = J0_EMF(omega_15N, taus, S2NCSAys, tauf, S2NCSAyf);
	R2CSAy = (1/45.) * d2y * (J0sum + 3 * J1);

	// XY
	J0sum = 0; J1 = 0;
	J0sum += (2/3.) * J0_EMF(2 * M_PI * (w1 - 2 * wr), taus, S2NCSAxys, tauf, S2NCSAxyf);
	J0sum += (2/3.) * J0_EMF(2 * M_PI * (w1 + 2 * wr), taus, S2NCSAxys, tauf, S2NCSAxyf);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 - wr), taus, S2NCSAxys, tauf, S2NCSAxyf);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 + wr), taus, S2NCSAxys, tauf, S2NCSAxyf);
	J1 = J0_EMF(omega_15N, taus, S2NCSAxys, tauf, S2NCSAxyf);
	R2CSAxy = (1/45.) * d2xy * (J0sum + 3 * J1);

	R2CSA = R2CSAx + R2CSAy - R2CSAxy;

	/* N Dipolar Interactions Contributions */

	long double J0HC, J2HC, J1H;

	// NH
	J0HC = J0_EMF(omega_1H - omega_15N, taus, S2NHs, tauf, S2NHf);
	J1   = J0_EMF(omega_15N, taus, S2NHs, tauf, S2NHf);
	J2HC = J0_EMF(omega_1H + omega_15N, taus, S2NHs, tauf, S2NHf);
	J1H  = J0_EMF(omega_1H, taus, S2NHs, tauf, S2NHf);
	J0sum  = (2/3.) * J0_EMF(2 * M_PI * (w1 + 2 * wr), taus, S2NHs, tauf, S2NHf);
	J0sum += (2/3.) * J0_EMF(2 * M_PI * (w1 - 2 * wr), taus, S2NHs, tauf, S2NHf);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 + wr), taus, S2NHs, tauf, S2NHf);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 - wr), taus, S2NHs, tauf, S2NHf);
	R2NH   = (1/20.) * sq(D_NH) * (J0sum + 3 * J1 + J0HC + 6 * J1H + 6 * J2HC);

	// NHr - the MATLAB appears to use S2CN
	J0HC = J0_EMF(omega_1H - omega_15N, taus, S2CNs, tauf, S2CNf);
	J1   = J0_EMF(omega_15N, taus, S2CNs, tauf, S2CNf);
	J2HC = J0_EMF(omega_1H + omega_15N, taus, S2CNs, tauf, S2CNf);
	J1H  = J0_EMF(omega_1H, taus, S2CNs, tauf, S2CNf);
	J0sum  = (2/3.) * J0_EMF(2 * M_PI * (w1 + 2 * wr), taus, S2CNs, tauf, S2CNf);
	J0sum += (2/3.) * J0_EMF(2 * M_PI * (w1 - 2 * wr), taus, S2CNs, tauf, S2CNf);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 + wr), taus, S2CNs, tauf, S2CNf);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 - wr), taus, S2CNs, tauf, S2CNf);
	R2NHr   = (1/20.) * sq(D_HNr) * (J0sum + 3 * J1 + J0HC + 6 * J1H + 6 * J2HC);

	// CN
	J0HC = J0_EMF(omega_13C - omega_15N, taus, S2CNs, tauf, S2CNf);
	J1   = J0_EMF(omega_15N, taus, S2CNs, tauf, S2CNf);
	J2HC = J0_EMF(omega_13C + omega_15N, taus, S2CNs, tauf, S2CNf);
	J1H  = J0_EMF(omega_13C, taus, S2CNs, tauf, S2CNf);
	J0sum  = (2/3.) * J0_EMF(2 * M_PI * (w1 + 2 * wr), taus, S2CNs, tauf, S2CNf);
	J0sum += (2/3.) * J0_EMF(2 * M_PI * (w1 - 2 * wr), taus, S2CNs, tauf, S2CNf);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 + wr), taus, S2CNs, tauf, S2CNf);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 - wr), taus, S2CNs, tauf, S2CNf);
	R2CN   = (1/20.) * sq(D_CN) * (J0sum + 3 * J1 + J0HC + 6 * J1H + 6 * J2HC);

	// CaN
	J0HC = J0_EMF(omega_13C - omega_15N, taus, S2CaNs, tauf, S2CaNf);
	J1   = J0_EMF(omega_15N, taus, S2CaNs, tauf, S2CaNf);
	J2HC = J0_EMF(omega_13C + omega_15N, taus, S2CaNs, tauf, S2CaNf);
	J1H  = J0_EMF(omega_13C, taus, S2CaNs, tauf, S2CaNf);
	J0sum  = (2/3.) * J0_EMF(2 * M_PI * (w1 + 2 * wr), taus, S2CaNs, tauf, S2CaNf);
	J0sum += (2/3.) * J0_EMF(2 * M_PI * (w1 - 2 * wr), taus, S2CaNs, tauf, S2CaNf);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 + wr), taus, S2CaNs, tauf, S2CaNf);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 - wr), taus, S2CaNs, tauf, S2CaNf);
	R2CaN   = (1/20.) * sq(D_CaN) * (J0sum + 3 * J1 + J0HC + 6 * J1H + 6 * J2HC);

	return R2CSA + R2NH + R2NHr + R2CN + R2CaN;
}

/* 	double R1CSAx, R1CSAy, R1CSAxy, R1CSA;
	long double J0sum = 0, J1 = 0;

	// X
	J0sum = 0; J1 = 0;
	J0sum += (2/3.) * J0_EMF(2 * M_PI * (w1 - 2 * wr), taus, S2NCSAxs, tauf, S2NCSAxf);
	J0sum += (2/3.) * J0_EMF(2 * M_PI * (w1 + 2 * wr), taus, S2NCSAxs, tauf, S2NCSAxf);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 - wr), taus, S2NCSAxs, tauf, S2NCSAxf);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 + wr), taus, S2NCSAxs, tauf, S2NCSAxf);
	J1 = J0_EMF(omega_15N, taus, S2NCSAxs, tauf, S2NCSAxf);
	R1CSAx = (2/15.) * d2x * J1;
	R2CSAx = (1/45.) * d2x * (J0sum + 3 * J1);

	// Y
	J0sum = 0; J1 = 0;
	J0sum += (2/3.) * J0_EMF(2 * M_PI * (w1 - 2 * wr), taus, S2NCSAys, tauf, S2NCSAyf);
	J0sum += (2/3.) * J0_EMF(2 * M_PI * (w1 + 2 * wr), taus, S2NCSAys, tauf, S2NCSAyf);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 - wr), taus, S2NCSAys, tauf, S2NCSAyf);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 + wr), taus, S2NCSAys, tauf, S2NCSAyf);
	J1 = J0_EMF(omega_15N, taus, S2NCSAys, tauf, S2NCSAyf);
	R1CSAx = (2/15.) * d2y * J1;
	R2CSAx = (1/45.) * d2y * (J0sum + 3 * J1);

	// XY
	J0sum = 0; J1 = 0;
	J0sum += (2/3.) * J0_EMF(2 * M_PI * (w1 - 2 * wr), taus, S2NCSAxys, tauf, S2NCSAxyf);
	J0sum += (2/3.) * J0_EMF(2 * M_PI * (w1 + 2 * wr), taus, S2NCSAxys, tauf, S2NCSAxyf);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 - wr), taus, S2NCSAxys, tauf, S2NCSAxyf);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 + wr), taus, S2NCSAxys, tauf, S2NCSAxyf);
	J1 = J0_EMF(omega_15N, taus, S2NCSAxys, tauf, S2NCSAxyf);
	R1CSAx = (2/15.) * d2xy * J1;
	R2CSAx = (1/45.) * d2xy * (J0sum + 3 * J1);

	long double J0sum = 0;
	J0sum = 0;
	J0sum += (2/3.) * J0_EMF(2 * M_PI * (w1 + 2 * wr), taus, S2NHs, tauf, S2NHf);
	J0sum += (2/3.) * J0_EMF(2 * M_PI * (w1 - 2 * wr), taus, S2NHs, tauf, S2NHf);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 + wr), taus, S2NHs, tauf, S2NHf);
	J0sum += (4/3.) * J0_EMF(2 * M_PI * (w1 - wr), taus, S2NHs, tauf, S2NHf);


	*/
