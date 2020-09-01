
/* This has not been directly tested against the MATLAB model as I don't really
 * have an equivalent to it set up (NHCO_3GAFSquaredEa is likely the closest
 * but I don't have the setup around it for the Legendre polynomials) */
 /* WARNING: NCSAxy IS OMITTED */
double GAF_15NR1_old(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double * sigf) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	double field = relax->field * 1000000; // conversion to Hz

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

	double R1CSAx, R1CSAy, R1CSAxy, R1CSA, R1NH, R1NHr, R1CN, R1CaN;
	long double J1 = 0;

	// X
	J1 = J0_EMF(omega_15N, taus, S2NCSAxs, tauf, S2NCSAxf);
	R1CSAx = (1/15.) * d2x * J1; // from Bremi1997

	// Y
	J1 = J0_EMF(omega_15N, taus, S2NCSAys, tauf, S2NCSAyf);
	R1CSAy = (1/15.) * d2y * J1;

	// XY
	J1 = J0_EMF(omega_15N, taus, S2NCSAxys, tauf, S2NCSAxyf);
	R1CSAxy = (2/15.) * d2xy * J1;
	/* CSA xy is the major contributor?! */
	R1CSA = R1CSAx + R1CSAy; //- R1CSAxy;

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
double GAF_15NR2_old(struct Residue *res, struct Relaxation* relax, long double taus, long double tauf, long double * sigs, long double * sigf) {
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

	R2CSA = R2CSAx + R2CSAy;// - R2CSAxy;

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
