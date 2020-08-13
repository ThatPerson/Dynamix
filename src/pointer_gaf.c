#include <stdio.h>
#include "datatypes.c"
#include "read_data.c"
#include "models.c"
//#include "chisq.c"
//#include "crosen.c" // implementation of Nelder-Mead simplex algorithm
//#include "errors.c"
#include <time.h>

double GAF_S2_current(long double sig[3], struct Orient * A, struct Orient * B, int mode) {
	int l, m, mp, k, kp;
	double complex Amp = 0;
	double complex temp;

	int ksqsum;
	int lsqsum;
	int msqsum;
	double lexp, kexp, mexp;
	for (l = -2; l <= 2; l++) {
		lsqsum = sq_i(l);
		lexp = -(sq(sig[1]) * lsqsum);
		for (m = -2; m <= 2; m++) {
			for (mp = -2; mp <= 2; mp++) {
				msqsum = sq_i(m) + sq_i(mp);
				mexp = -(sq(sig[2]) * msqsum/2.);
				for (k = -2; k <= 2; k++) {
					for (kp = -2; kp <= 2; kp++) {
						ksqsum = sq_i(k) + sq_i(kp);
						kexp = -(sq(sig[0]) * ksqsum / 2.);
						//if (A->Y2[m+2] == A->Y2c[m+2] && B->Y2[mp+2] == B->Y2c[mp+2]) {
						if (cimag(A->Y2[m+2]) < 0.00001 && cimag(B->Y2[mp+2]) < 0.00001 && mode==MODE_REAL) {
							/* In this case then Y2 and Y2c are real numbers and therefore the last step of this is real.
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
							temp = cpow(-1 * I, k - kp);
						}
						temp *= exp(lexp + kexp + mexp);
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
		case MODE_REAL: return creal(Amp); break;
		case MODE_IMAG: return cimag(Amp); break;
		default: break;
	}
	return -1;
}


int GAF_S2_new(long double sig[3], struct Orient ** A, struct Orient ** B, double * S2[], int length,  int mode) {
	int l, m, mp, k, kp, i;
	double complex * Amp = (double complex *) malloc(sizeof(double complex) * length);
	double complex temp, ttemp;
	for (i = 0; i < length; i++) {
		Amp[i] = 0;
	}
	
	int ksqsum;
	int lsqsum;
	int msqsum;
	double lexp, kexp, mexp;
	for (l = -2; l <= 2; l++) {
		lsqsum = sq_i(l);
		lexp = -(sq(sig[1]) * lsqsum);
		for (m = -2; m <= 2; m++) {
			for (mp = -2; mp <= 2; mp++) {
				msqsum = sq_i(m) + sq_i(mp);
				mexp = -(sq(sig[2]) * msqsum/2.);
				for (k = -2; k <= 2; k++) {
					for (kp = -2; kp <= 2; kp++) {
						ksqsum = sq_i(k) + sq_i(kp);
						kexp = -(sq(sig[0]) * ksqsum / 2.);
						//if (A->Y2[m+2] == A->Y2c[m+2] && B->Y2[mp+2] == B->Y2c[mp+2]) {
						temp = 1;
						
						temp *= exp(lexp + kexp + mexp);
						temp *= Dwig[k+2][l+2] * Dwig[kp+2][l+2] * Dwig[m+2][k+2] * Dwig[mp+2][kp+2];
						for (i = 0; i < length; i++) {
							ttemp = temp;
							if (cimag(A[i]->Y2[m+2]) < 0.00001 && cimag(B[i]->Y2[mp+2]) < 0.00001 && mode==MODE_REAL) {
								/* In this case then Y2 and Y2c are real numbers and therefore the last step of this is real.
								* If mode is MODE_REAL, then, we may safely ignore any case in which
								*  cpowl(-1*I, k-kp)
								* is imaginary (as this will give only an imaginary component to the Amp).
								* This is the case for any (k-kp)%2 != 0.
								* In this case, if (k-kp)%4 == 0 the cpowl function gives 1, else -1.
								*/
								if ((k-kp)%2 != 0)
									continue;
								if ((k - kp) % 4 == 0)
									ttemp *= 1;
								else
									ttemp *= -1.;
							} else {
								/* Other wise we have to go the long route... */
								ttemp *= cpow(-1 * I, k - kp);
							}
							
							
							
							ttemp *= A[i]->Y2[m+2] * B[i]->Y2c[mp+2];
							Amp[i] += ttemp;
						}
						//temp *= A->Y2[m+2] * B->Y2c[mp+2];
						//Amp += temp;
					}
				}
			}
		}
	}
	//Amp *= (4 * M_PI / 5.);
	for (i = 0;i < length; i++) {
		Amp[i] = Amp[i] * (4 * M_PI / 5.);
		switch (mode) {
			case MODE_REAL: *S2[i] = creal(Amp[i]); break;
			case MODE_IMAG: *S2[i] = creal(Amp[i]); break;
			default: break;
		}
	}
	free(Amp);
	return 1;
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

    int residue = 20;
    int rel = 13;
	struct Residue *res = &(m.residues[residue]);
    long double sigs[3], sigf[3];
    long double taus, tauf;
    double S2s_eff, S2f_eff;
    double R1GAF, R2GAF, R1GAF_old, R2GAF_old;
    float s, f;
	double S2NHs, S2NCSAxs, S2NCSAys, S2NCSAxys, S2CNs, S2CaNs;
	double S2NHf, S2NCSAxf, S2NCSAyf, S2NCSAxyf, S2CNf, S2CaNf;
	double *S2[] = {&S2NHf, &S2NCSAxf, &S2NCSAyf, &S2NCSAxyf, &S2CNf, &S2CaNf};
	struct Orient * As[] = {&(res->orients[OR_NH]), &(res->orients[OR_NCSAxx]), &(res->orients[OR_NCSAyy]), &(res->orients[OR_NCSAxx]), &(res->orients[OR_CN]), &(res->orients[OR_NCA])};
	struct Orient * Bs[] = {&(res->orients[OR_NH]), &(res->orients[OR_NCSAxx]), &(res->orients[OR_NCSAyy]), &(res->orients[OR_NCSAyy]), &(res->orients[OR_CN]), &(res->orients[OR_NCA])};
	double Sums = 0, Sumf = 0;
	int count = 0;
	for (sigs[0] = 0.00; sigs[0] <= 0.1; sigs[0] += 0.01) {
		for (sigs[1] = 0.00; sigs[1] <= 0.1; sigs[1] += 0.01) {
			for (sigs[2] = 0.00; sigs[2] <= 0.2; sigs[2] += 0.01) {
				//printf("(%Lf, %Lf, %Lf)\n", sigs[0], sigs[1], sigs[2]);
				/*S2NHs    = GAF_S2_current(sigs, &(res->orients[OR_NH]), &(res->orients[OR_NH]), MODE_REAL);
				S2NCSAxs = GAF_S2_current(sigs, &(res->orients[OR_NCSAxx]), &(res->orients[OR_NCSAxx]), MODE_REAL);
				S2NCSAys = GAF_S2_current(sigs, &(res->orients[OR_NCSAyy]), &(res->orients[OR_NCSAyy]), MODE_REAL);
				S2NCSAxys= GAF_S2_current(sigs, &(res->orients[OR_NCSAxx]), &(res->orients[OR_NCSAyy]), MODE_REAL);
				S2CNs    = GAF_S2_current(sigs, &(res->orients[OR_CN]), &(res->orients[OR_CN]), MODE_REAL);
				S2CaNs   = GAF_S2_current(sigs, &(res->orients[OR_NCA]), &(res->orients[OR_NCA]), MODE_REAL);*/
				
				GAF_S2_new(sigs, As, Bs, S2, 6, MODE_REAL);
				//int GAF_S2_new(long double sig[3], struct Orient ** A, struct Orient ** B, double * S2[], int length,  int mode) {
				count += 6;
				Sums += S2NHs;
				Sumf += S2NHf;
				
				
				
				//printf("\t%f, %f, %f, %f, %f, %f\n", S2NHs, S2NCSAxs, S2NCSAys, S2NCSAxys, S2CNs, S2CaNs);
				//printf("\t%f, %f, %f, %f, %f, %f\n", S2NHf, S2NCSAxf, S2NCSAyf, S2NCSAxyf, S2CNf, S2CaNf);
				/*printf("\t[%f, %f, %f, %f]\n\t\t%f\n", res->orients[OR_NH].phi, res->orients[OR_NH].theta, res->orients[OR_NH].phi, res->orients[OR_NH].theta, S2NHs);
				
				printf("\t[%f, %f, %f, %f]\n\t\t%f\n", res->orients[OR_NCSAxx].phi, res->orients[OR_NCSAxx].theta, res->orients[OR_NCSAxx].phi, res->orients[OR_NCSAxx].theta, S2NCSAxs);
				
				
				printf("\t[%f, %f, %f, %f]\n\t\t%f\n", res->orients[OR_NCSAyy].phi, res->orients[OR_NCSAyy].theta, res->orients[OR_NCSAyy].phi, res->orients[OR_NCSAyy].theta, S2NCSAys);
				
				printf("\t[%f, %f, %f, %f]\n\t\t%f\n", res->orients[OR_NCSAxx].phi, res->orients[OR_NCSAxx].theta, res->orients[OR_NCSAyy].phi, res->orients[OR_NCSAyy].theta, S2NCSAxys);
				
				printf("\t[%f, %f, %f, %f]\n\t\t%f\n", res->orients[OR_CN].phi, res->orients[OR_CN].theta, res->orients[OR_CN].phi, res->orients[OR_CN].theta, S2CNs);
				
				printf("\t[%f, %f, %f, %f]\n\t\t%f\n", res->orients[OR_NCA].phi, res->orients[OR_NCA].theta, res->orients[OR_NCA].phi, res->orients[OR_NCA].theta, S2CaNs);
				printf("\n\n");*/

			}
		}
	}
	printf("Count: %d\nSums: %f\nSumf: %f\n", count, Sums, Sumf);
	
    /*for (s = 0.1; s <= 0.3; s += 0.1) {
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
                                    //R2GAF = GAF_15NR2(&(m.residues[residue]), &(m.residues[residue].relaxation[rel]), taus, tauf, sigs, sigf);
                                    //R2GAF_old = GAF_15NR2_old(&(m.residues[residue]), &(m.residues[residue].relaxation[rel]), taus, tauf, sigs, sigf);

                                    //if (R2GAF != R2GAF_old)
                                    //    printf("Oops\n");
                                    //R2GAF = GAF_15NR2(&(m.residues[residue]), &(m.residues[residue].relaxation[rel]), taus, tauf, sigs, sigf);



                                    //printf("R1\t\tnew: %f\n\t\told: %f\n", R1GAF, R1GAF_old);
                                //    printf("R2\t\tnew: %f\n\t\told: %f\n", R2GAF, R2GAF_old);


                               /* }
                            }
                        }
                    }
                }
            }
        }
    }*/





    return 1;
}
