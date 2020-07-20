#include <stdio.h>
#include <math.h>

//#define J0_SMF(omega) (((1 - S2) * tau) / (1 + pow(omega * tau, 2)))

//#define J0_SMF_num(omega) ((1 - S2) * tau)
//#define J0_SMF_den(omega) (1 + pow(omega * tau, 2))

long double J0_SMF(double omega, long double * tau, long double * S2) {
	return (((1 - (*S2)) * (*tau)) / (1 + pow((omega) * (*tau), 2)));
}

double SMF_15NR1(struct Residue *res, struct Relaxation* relax, long double tau, long double S2) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	
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
	Jcomp += J0_SMF(omega_1H - omega_15N, &tau, &S2);
	Jcomp += 3 * J0_SMF(omega_15N, &tau, &S2);
	Jcomp += 6 * J0_SMF(omega_1H + omega_15N, &tau, &S2);
	R1NH = 0.1 * d * d * Jcomp;
	R1NCSA = (2/15.) * d2 * d2 * J0_SMF(omega_15N, &tau, &S2);

	return R1NH + R1NCSA;
}

double SMF_15NR2(struct Residue *res, struct Relaxation* relax, long double tau, long double S2) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	
	long double field = relax->field * 1000000; // conversion to Hz
	
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
	J0sum += (2/3.) * J0_SMF(2 * M_PI * (w1 - 2 * wr), &tau, &S2);
	J0sum += (2/3.) * J0_SMF(2 * M_PI * (w1 + 2 * wr), &tau, &S2);
	J0sum += (4/3.) * J0_SMF(2 * M_PI * (w1 - wr), &tau, &S2);
	J0sum += (4/3.) * J0_SMF(2 * M_PI * (w1 + wr), &tau, &S2);
	
	long double JNH = J0sum + 3 * J0_SMF(omega_15N, &tau, &S2);
	JNH += J0_SMF(omega_1H - omega_15N, &tau, &S2);
	JNH += 6 * J0_SMF(omega_1H, &tau, &S2);
	JNH += 6 * J0_SMF(omega_1H + omega_15N, &tau, &S2);
	
	
	R2NH = (1/20.) * d * d * JNH;
	R2NCSA = (1/45.) * d2 * d2 * (J0sum + 3 * J0_SMF(omega_15N, &tau, &S2));	
	return R2NH + R2NCSA;
}
