#include <stdio.h>
#include "datatypes.c"
#include "read_data.c"
#include "models.c"
#include "chisq.c"
#include "crosen.c" // implementation of Nelder-Mead simplex algorithm
#include "errors.c"
#include <time.h>

double EMF_R1_f(struct Residue *res, struct Relaxation* relax, long double taus, long double S2s, long double tauf, long double S2f, int mode) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	//long double S2f = res->S2_dipolar / S2s;
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

    long double J1 = 0;
    J1 = J0_EMF(omega_L, taus, S2s, tauf, S2f);

    printf("\t\tEMF\td2tot: \t%Lf\n\t\t\tJ1: \t%0.36Lf\n", d2tot, J1);
	R1CSA = (2/15.) * d2tot * J1;


	return R1CSA + R1CSA;
}

double EMF_R2_f(struct Residue *res, struct Relaxation* relax, long double taus, long double S2s, long double tauf, long double S2f, int mode) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	//long double S2f = res->S2_dipolar / S2s;
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
	return R2CSA + R2D;
}

int main(int argc, char * argv[]) {

	/* Initialisation */
	srand(time(NULL));
	initialise_dwig();

	char system_file[255] = "system/model.dx";
	int i;
	int err_mod = 0;

	struct Model m;
	int ret = read_system_file(system_file, &m);
	m.error_mode = err_mod;

	if (ret == -1) {
		printf("Error found, crashing peacefully...\n");
		exit(-1);
	}


    /* Throw random numbers into R1_EMF, R2_EMF and compare with R1_GAF, R2_GAF */

    /*double GAF_15NR1(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double * sigf) {
    double GAF_15NR2(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double * sigf) {
    double GAF_S2(long double sig[3], struct Orient * A, struct Orient * B, int mode) {
    double EMF_R2(struct Residue *res, struct Relaxation* relax, long double taus, long double S2s, long double tauf, long double S2f, int mode) {
    double EMF_R1(struct Residue *res, struct Relaxation* relax, long double taus, long double S2s, long double tauf, long double S2f, int mode) {
    MODE_15N
    MODE_REAL
    OR_NH*/

    int residue = 40;
    int rel = 10;

    long double sigs[3], sigf[3];
    long double taus, tauf;
    double S2s_eff, S2f_eff;
    double R1EMF, R2EMF, R1GAF, R2GAF;
    float s, f;
    for (s = 0.1; s <= 1; s += 0.1) {
        for (f = 0.1; s <= 1; s += 0.1) {
            taus = s * powl(10, -8);
            tauf = f * powl(10, -11);
            for (sigs[0] = 0; sigs[0] <= 0.1; sigs[0] += 0.05) {
                for (sigs[1] = 0; sigs[1] <= 0.1; sigs[1] += 0.05) {
                    for (sigs[2] = 0; sigs[2] <= 0.1; sigs[2] += 0.05) {
                        for (sigf[0] = 0; sigf[0] <= 0.1; sigf[0] += 0.05) {
                            for (sigf[1] = 0; sigf[1] <= 0.1; sigf[1] += 0.05) {
                                for (sigf[2] = 0; sigf[2] <= 0.1; sigf[2] += 0.05) {
                                    S2s_eff = GAF_S2(sigs, &(m.residues[residue].orients[OR_NH]), &(m.residues[residue].orients[OR_NH]), MODE_REAL);
                                    S2f_eff = GAF_S2(sigf, &(m.residues[residue].orients[OR_NH]), &(m.residues[residue].orients[OR_NH]), MODE_REAL);
                                    printf("Tauf: %Le; Taus: %Le\n", tauf, taus);
                                    printf("Slow [%Lf, %Lf, %Lf] -> %f\n", sigs[0], sigs[1], sigs[2], S2s_eff);
                                    printf("Fast [%Lf, %Lf, %Lf] -> %f\n", sigf[0], sigf[1], sigf[2], S2f_eff);
                                    R1EMF = EMF_R1_f(&(m.residues[residue]), &(m.residues[residue].relaxation[rel]), taus, S2s_eff, tauf, S2f_eff, MODE_15N);
                                    R2EMF = EMF_R2_f(&(m.residues[residue]), &(m.residues[residue].relaxation[rel]), taus, S2s_eff, tauf, S2f_eff, MODE_15N);

                                    R1GAF = GAF_15NR1(&(m.residues[residue]), &(m.residues[residue].relaxation[rel]), taus, tauf, sigs, sigf);
                                    R2GAF = GAF_15NR2(&(m.residues[residue]), &(m.residues[residue].relaxation[rel]), taus, tauf, sigs, sigf);


                                    printf("R1\t\tEMF: %f\n\t\tGAF: %f\n", R1EMF, R1GAF);
                                    printf("R2\t\tEMF: %f\n\t\tGAF: %f\n", R2EMF, R2GAF);
                                }
                            }
                        }
                    }
                }
            }
        }
    }





    return 1;
}
