#include <stdio.h>
#include "datatypes.c"
#include "read_data.c"
#include "models.c"
#include "models_old.c"
//#include "chisq.c"
//#include "crosen.c" // implementation of Nelder-Mead simplex algorithm
//#include "errors.c"
#include <time.h>

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

    int residue = 20;
    int rel = 13;

    long double sigs[3], sigf[3];
    long double taus, tauf;
    double S2s_eff, S2f_eff;
    double R1GAF, R2GAF, R1GAF_old, R2GAF_old;
    float s, f;
    for (s = 0.1; s <= 0.3; s += 0.1) {
        for (f = 0.1; f <= 0.3; f += 0.1) {
            taus = s * powl(10, -8);
            tauf = f * powl(10, -11);
            for (sigs[0] = 0.05; sigs[0] <= 0.1; sigs[0] += 0.05) {
                for (sigs[1] = 0.05; sigs[1] <= 0.1; sigs[1] += 0.05) {
                    for (sigs[2] = 0.05; sigs[2] <= 0.1; sigs[2] += 0.05) {
                        for (sigf[0] = 0.05; sigf[0] <= 0.1; sigf[0] += 0.05) {
                            for (sigf[1] = 0.05; sigf[1] <= 0.1; sigf[1] += 0.05) {
                                for (sigf[2] = 0.05; sigf[2] <= 0.1; sigf[2] += 0.05) {
                                    /*S2s_eff = GAF_S2(sigs, &(m.residues[residue].orients[OR_NH]), &(m.residues[residue].orients[OR_NH]), MODE_REAL);
                                    S2f_eff = GAF_S2(sigf, &(m.residues[residue].orients[OR_NH]), &(m.residues[residue].orients[OR_NH]), MODE_REAL);
                                    printf("Tauf: %Le; Taus: %Le\n", tauf, taus);
                                    printf("Slow [%Lf, %Lf, %Lf] -> %f\n", sigs[0], sigs[1], sigs[2], S2s_eff);
                                    printf("Fast [%Lf, %Lf, %Lf] -> %f\n", sigf[0], sigf[1], sigf[2], S2f_eff);
                                    R1EMF = EMF_R1_f(&(m.residues[residue]), &(m.residues[residue].relaxation[rel]), taus, S2s_eff, tauf, S2f_eff, MODE_15N);
                                    R2EMF = EMF_R2_f(&(m.residues[residue]), &(m.residues[residue].relaxation[rel]), taus, S2s_eff, tauf, S2f_eff, MODE_15N);
                                    */
                                    R2GAF = GAF_15NR2(&(m.residues[residue]), &(m.residues[residue].relaxation[rel]), taus, tauf, sigs, sigf);
                                    //R2GAF_old = GAF_15NR2_old(&(m.residues[residue]), &(m.residues[residue].relaxation[rel]), taus, tauf, sigs, sigf);

                                    //if (R2GAF != R2GAF_old)
                                    //    printf("Oops\n");
                                    //R2GAF = GAF_15NR2(&(m.residues[residue]), &(m.residues[residue].relaxation[rel]), taus, tauf, sigs, sigf);



                                    //printf("R1\t\tnew: %f\n\t\told: %f\n", R1GAF, R1GAF_old);
                                //    printf("R2\t\tnew: %f\n\t\told: %f\n", R2GAF, R2GAF_old);


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
