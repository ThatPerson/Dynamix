//
// Created by ben on 16/03/2021.
//

#include <stdlib.h>
#include <complex.h>
#include "model.h"

/** Calculates spectral density function for given frequency according to extended model free analysis.
 *  Implements J(w) = (1 - S2f) tauf / (1 + (w tauf)^2) + S2f (1 - S2s) taus / (1 + (w taus)^2),
 *  as in eq 2, Clore 1990 (assuming S2 = S2s*S2f and that there is no overall tumbling).
 */
Decimal J0(Decimal omega, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f, Decimal S2uf) {
	return S2uf * (( \
			((((1 - (Decimal) S2f)) * (Decimal) tauf)) \
			/ (1 + ((Decimal) omega * (Decimal) tauf * (Decimal) omega * (Decimal) tauf)) \
		)\
		+\
		(\
			(((Decimal) S2f) * (1 - (Decimal) S2s) * (Decimal) taus)\
			/ (1 + ((Decimal) omega * (Decimal) taus * (Decimal) omega * (Decimal) taus))\
		));
}

/** Cross correlated spectral density. Taken from Manuscript,
 * Juv(w) = (1-S2f) tf / (1 + (wtf)^2) + (1/P2) S2f (P2 - S2s) (ts / 1+(wts)^2)
 * eq 2
 */
Decimal J0_CC(Decimal omega, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f, Decimal S2uf) {
	return S2uf * (( \
			((((1 - (Decimal) S2f)) * (Decimal) tauf)) \
			/ (1 + ((Decimal) omega * (Decimal) tauf * (Decimal) omega * (Decimal) tauf)) \
		)\
		-\
		2. * (\
			(((Decimal) S2f) * (-0.5 - (Decimal) S2s) * (Decimal) taus)\
			/ (1 + ((Decimal) omega * (Decimal) taus * (Decimal) omega * (Decimal) taus))\
		));
}


Decimal Dipolar_R1(Decimal omega_obs, Decimal omega_neigh, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f, Decimal S2uf, Decimal D) {
	// for SMF and SMFT models, set S2s = 0 and taus = 0.
	Decimal q = (0.1) * sq(D) * ( \
		(J0(omega_neigh - omega_obs, taus, S2s, tauf, S2f, S2uf)) \
		 + 3 * (J0(omega_obs, taus, S2s, tauf, S2f, S2uf)) \
		  + 6 * (J0(omega_neigh + omega_obs, taus, S2s, tauf, S2f, S2uf)) \
	  );
	return (Decimal) q;
}


Decimal Dipolar_R2(Decimal omega_obs, Decimal omega_neigh, Decimal w1, Decimal wr, Decimal taus, Decimal S2s, Decimal tauf, Decimal S2f, Decimal S2uf, Decimal D) {
	// for SMF and SMFT models, set S2s = 0 and taus = 0.
	return (Decimal) (\
  		(1/20.) * sq(D) * (\
  			(2/3.) * J0(2 * M_PI * (w1 + 2 * wr), taus, S2s, tauf, S2f, S2uf) +\
  			(2/3.) * J0(2 * M_PI * (w1 - 2 * wr), taus, S2s, tauf, S2f, S2uf) +\
  			(4/3.) * J0(2 * M_PI * (w1 + wr), taus, S2s, tauf, S2f, S2uf) +\
  			(4/3.) * J0(2 * M_PI * (w1 - wr), taus, S2s, tauf, S2f, S2uf) +\
  			(3.)   * J0(omega_obs, taus, S2s, tauf, S2f, S2uf) +\
  			(1.)   * J0(omega_neigh - omega_obs, taus, S2s, tauf, S2f, S2uf) +\
  			(6.)   * J0(omega_neigh, taus, S2s, tauf, S2f, S2uf) +\
  			(6.)   * J0(omega_neigh + omega_obs, taus, S2s, tauf, S2f, S2uf)\
  		)\
  	);
}

/**
 * Calculates the R2 contribution of CSA of nuclei with frequency omega.
 * D2 is the squared D22x/D22y/D22xy variable (dependent on CSA)
 * Kurbanov2011, eqs 18-20.
 */
Decimal CSA_R2(Decimal omega, \
						 Decimal w1, \
						 Decimal wr, \
						 Decimal taus, \
						 Decimal S2s, \
						 Decimal tauf, \
						 Decimal S2f, \
						 Decimal S2uf, \
						 Decimal D2, \
						 Decimal  (*J_SD)(\
						 	Decimal,\
						 	Decimal, \
						 	Decimal, \
						 	Decimal, \
						 	Decimal, \
							Decimal
						)\
						) {
	return (Decimal) ( \
  		(1/45.) * (Decimal) D2 * ( \
  			(2/3.) * J_SD(2 * M_PI * (w1 - 2 * wr), taus, S2s, tauf, S2f, S2uf) +\
  			(2/3.) * J_SD(2 * M_PI * (w1 + 2 * wr), taus, S2s, tauf, S2f, S2uf) +\
  			(4/3.) * J_SD(2 * M_PI * (w1 - wr), taus, S2s, tauf, S2f, S2uf) +\
  			(4/3.) * J_SD(2 * M_PI * (w1 + wr), taus, S2s, tauf, S2f, S2uf) + \
  			(3.)   * J_SD(omega, taus, S2s, tauf, S2f, S2uf)\
  		)\
  	);
}

Decimal Calc_15NR1(struct Residue *res, struct Relaxation *relax, struct BCParameters *pars, struct Model *m) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	Decimal field = relax->field * 1000000; // conversion to Hz

	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a Decimal128, and
	 * software division) I've predefined it. Bond length taken as 1.02 A */
	unsigned int model = m->model;

	Decimal taus, tauf;
	if (pars->Eas != -1 || pars->Eaf != -1) {
		taus = temp_tau(pars->taus, pars->Eas, relax->T);
		tauf = temp_tau(pars->tauf, pars->Eaf, relax->T);
	} else {
		taus = pars->taus;
		tauf = pars->tauf;
	}

	Decimal omega_1H = T_DOWN * 2 * M_PI * field;
	Decimal omega_15N = T_DOWN * 2 * M_PI * field / 9.869683408806043;
	Decimal omega_13C = T_DOWN * 2 * M_PI * field / 3.976489314034722;
	Decimal d2x, d2y, d2xy, d2tot;

	Decimal *csa;
	csa = res->csaN;

	/* N CSA relaxation contribution */
	Decimal R1CSAx, R1CSAy, R1CSAxy, R1CSA, R1NH, R1NHr, R1CN, R1CaN;
	Decimal J1;

	if (model == MOD_GAF || model == MOD_GAFT || model == MOD_EGAF || model == MOD_EGAFT || model == MOD_BGF || model == MOD_BGFT) {
		d2x = (Decimal) sq(((csa[2] - csa[0]) * 0.000001) * omega_15N * T_UP);
		d2y = (Decimal) sq(((csa[1] - csa[0]) * 0.000001) * omega_15N * T_UP);
		d2xy = (Decimal) sq(0.000001 * omega_15N * T_UP) * (csa[2] - csa[0]) * (csa[1] - csa[0]);

		J1 = J0(omega_15N, taus, pars->S2NCSAxs, tauf, pars->S2NCSAxf, pars->S2uf);
		R1CSAx = (1 / 15.) * d2x * J1; // from Bremi1997
		J1 = J0(omega_15N, taus, pars->S2NCSAys, tauf, pars->S2NCSAyf, pars->S2uf);
		R1CSAy = (1 / 15.) * d2y * J1;
		J1 = J0_CC(omega_15N, taus, pars->S2NCSAxys, tauf, pars->S2NCSAxyf, pars->S2uf);
		R1CSAxy = (1 / 15.) * d2xy * J1;
		/** Eq A30, Bremi1997 */
		R1CSA = R1CSAx + R1CSAy + 2. * R1CSAxy;
	} else {
		d2tot= (sq(csa[2]) + sq(csa[1]) + sq(csa[0]));
		d2tot+=(-(csa[2] * csa[1] + csa[1] * csa[0] + csa[2] * csa[0]));
		d2tot*=sq(0.000001 * omega_15N * T_UP);
		R1CSA = (2/15.) * d2tot * J0(omega_15N, taus, pars->S2NCSAxs, tauf, pars->S2NCSAxf, pars->S2uf);
	}
	/* N Dipolar Interactions Contributions */
	R1NH = Dipolar_R1(omega_15N, omega_1H, taus, pars->S2NHs, tauf, pars->S2NHf, pars->S2uf, D_NH);
	R1NHr = Dipolar_R1(omega_15N, omega_1H, taus, pars->S2NHrs, tauf, pars->S2NHrf, pars->S2uf, D_NHr);
	R1CN = Dipolar_R1(omega_15N, omega_13C, taus, pars->S2CNs, tauf, pars->S2CNf, pars->S2uf, D_CN);
	R1CaN = Dipolar_R1(omega_15N, omega_13C, taus, pars->S2CaNs, tauf, pars->S2CaNf, pars->S2uf, D_NCA);


	return (Decimal) (R1CSA + R1NH + R1NHr + R1CN + R1CaN) * T_DOWN;
}

Decimal Calc_15NR2(struct Residue *res, struct Relaxation* relax, struct BCParameters *pars, struct Model *m) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	Decimal field = relax->field * 1000000; // conversion to Hz
	Decimal w1 = relax->w1, wr = relax->wr;
	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a Decimal128, and
	 * software division) I've predefined it. Bond length taken as 1.02 A */
	unsigned int model = m->model;

	Decimal taus, tauf;
	if (pars->Eas != -1 || pars->Eaf != -1) {
		taus = temp_tau(pars->taus, pars->Eas, relax->T);
		tauf = temp_tau(pars->tauf, pars->Eaf, relax->T);
	} else {
		taus = pars->taus;
		tauf = pars->tauf;
	}

	Decimal omega_1H = 2 * M_PI * field;
	Decimal omega_13C, omega_15N;
	Decimal d2x, d2y, d2xy, d2tot;

	omega_13C = 2 * M_PI * field / 3.976489314034722;
	omega_15N = 2 * M_PI * field / 9.869683408806043;

	Decimal *csa;
	csa = res->csaN;

	d2x = (Decimal) sq(((csa[2] - csa[0]) * 0.000001) * omega_15N);
	d2y = (Decimal) sq(((csa[1] - csa[0]) * 0.000001) * omega_15N);
	d2xy= (Decimal) sq(0.000001 * omega_15N) * (csa[2] - csa[0]) * (csa[1] - csa[0]);

	/* N CSA relaxation contribution */
	Decimal R2CSAx, R2CSAy, R2CSAxy, R2CSA, R2NH, R2NHr, R2CN, R2CaN;

	w1 *= T_DOWN;
	wr *= T_DOWN;
	omega_1H *= T_DOWN;
	omega_13C *= T_DOWN;
	omega_15N *= T_DOWN;

	Decimal J0sum = 0;


	if (model == MOD_GAF || model == MOD_GAFT || model == MOD_EGAF || model == MOD_EGAFT || model == MOD_BGF || model == MOD_BGFT) {
		R2CSAx = CSA_R2(omega_15N, w1, wr, taus, pars->S2NCSAxs, tauf, pars->S2NCSAxf, pars->S2uf, d2x, J0);
		R2CSAy = CSA_R2(omega_15N, w1, wr, taus, pars->S2NCSAys, tauf, pars->S2NCSAyf, pars->S2uf, d2y, J0);
		R2CSAxy = CSA_R2(omega_15N, w1, wr, taus, pars->S2NCSAxys, tauf, pars->S2NCSAxyf, pars->S2uf, d2xy, J0_CC);
		/**
		 * Equation for spectral density for xy taken from Manuscript
		 */
		R2CSA = R2CSAx + R2CSAy + 2. * R2CSAxy;
	} else {
		J0sum += (2/3.) * J0(2 * M_PI * (w1 - 2 * wr), taus, pars->S2NCSAxs, tauf, pars->S2NCSAxf, pars->S2uf);
		J0sum += (2/3.) * J0(2 * M_PI * (w1 + 2 * wr), taus, pars->S2NCSAxs, tauf, pars->S2NCSAxf, pars->S2uf);
		J0sum += (4/3.) * J0(2 * M_PI * (w1 - wr), taus, pars->S2NCSAxs, tauf, pars->S2NCSAxf, pars->S2uf);
		J0sum += (4/3.) * J0(2 * M_PI * (w1 + wr), taus, pars->S2NCSAxs, tauf, pars->S2NCSAxf, pars->S2uf);

		d2tot= (sq(csa[2]) + sq(csa[1]) + sq(csa[0]));
		d2tot+=(-(csa[2] * csa[1] + csa[1] * csa[0] + csa[2] * csa[0]));
		d2tot*=sq(0.000001 * omega_15N * T_UP);
		R2CSA = (1/45.) * d2tot * (J0sum + 3 * J0(omega_15N, taus, pars->S2NCSAxs, tauf, pars->S2NCSAxf, pars->S2uf));

	}

	/* N Dipolar Interactions Contributions */
	R2NH = Dipolar_R2(omega_15N, omega_1H, w1, wr, taus, pars->S2NHs, tauf, pars->S2NHf, pars->S2uf, D_NH);
	R2NHr = Dipolar_R2(omega_15N, omega_1H, w1, wr, taus, pars->S2NHrs, tauf, pars->S2NHrf, pars->S2uf, D_NHr);
	R2CN = Dipolar_R2(omega_15N, omega_13C, w1, wr, taus, pars->S2CNs, tauf, pars->S2CNf, pars->S2uf, D_CN);
	R2CaN = Dipolar_R2(omega_15N, omega_13C, w1, wr, taus, pars->S2CaNs, tauf, pars->S2CaNf, pars->S2uf, D_NCA);
	return (R2CSA + R2NH + R2NHr + R2CN + R2CaN)*T_DOWN;
}

Decimal Calc_13CR1(struct Residue *res, struct Relaxation* relax, struct BCParameters *pars, struct Model *m) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	Decimal field = relax->field * 1000000; // conversion to Hz

	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a Decimal128, and
	 * software division) I've predefined it. Bond length taken as 1.02 A */
	unsigned int model = m->model;

	Decimal taus, tauf;
	if (pars->Eas != -1 || pars->Eaf != -1) {
		taus = temp_tau(pars->taus, pars->Eas, relax->T);
		tauf = temp_tau(pars->tauf, pars->Eaf, relax->T);
	} else {
		taus = pars->taus;
		tauf = pars->tauf;
	}

	Decimal omega_1H = 2 * M_PI * field;
	Decimal omega_13C, omega_15N, wCOCa;
	Decimal d2x, d2y, d2xy, d2tot;

	omega_13C = 2 * M_PI * field / 3.976489314034722;
	omega_15N = 2 * M_PI * field / 9.869683408806043;
	wCOCa = 120 * omega_13C * 0.000001;

	omega_1H *= T_DOWN;
	omega_13C *= T_DOWN;
	omega_15N *= T_DOWN;
	wCOCa *= T_DOWN;

	Decimal *csa;
	csa = res->csaC;

/* N CSA relaxation contribution */
	Decimal R1CSAx, R1CSAy, R1CSAxy, R1CSA, R1CH, R1CHr, R1CN, R1CCAc, R1CCAp;
	Decimal J1;
	if (model == MOD_GAF || model == MOD_GAFT || model == MOD_EGAF || model == MOD_EGAFT || model == MOD_BGF || model == MOD_BGFT) {
		d2x = (Decimal) sq(((csa[2] - csa[0]) * 0.000001) * omega_13C * T_UP);
		d2y = (Decimal) sq(((csa[1] - csa[0]) * 0.000001) * omega_13C * T_UP);
		d2xy= (Decimal) sq(0.000001 * omega_13C * T_UP) * (csa[2] - csa[0]) * (csa[1] - csa[0]);
		J1 = J0(omega_13C, taus, pars->S2CSAxs, tauf, pars->S2CSAxf, pars->S2uf);
		R1CSAx = (1/15.) * d2x * J1; // from Bremi1997
		J1 = J0(omega_13C, taus, pars->S2CSAys, tauf, pars->S2CSAyf, pars->S2uf);
		R1CSAy = (1/15.) * d2y * J1;
		J1 = J0_CC(omega_13C, taus, pars->S2CSAxys, tauf, pars->S2CSAxyf, pars->S2uf);
		R1CSAxy = (1/15.) * d2xy * J1;
		R1CSA = R1CSAx + R1CSAy + 2.*R1CSAxy;
	} else {
		d2tot= (sq(csa[2]) + sq(csa[1]) + sq(csa[0]));
		d2tot+=(-(csa[2] * csa[1] + csa[1] * csa[0] + csa[2] * csa[0]));
		d2tot*=sq(0.000001 * omega_13C * T_UP);
		R1CSA = (2/15.) * d2tot * J0(omega_13C, taus, pars->S2CSAxs, tauf, pars->S2CSAxf, pars->S2uf);
	}



	/* N Dipolar Interactions Contributions */
	R1CH = Dipolar_R1(omega_13C, omega_1H, taus, pars->S2CHs, tauf, pars->S2CHf, pars->S2uf, D_CH);
	R1CHr = Dipolar_R1(omega_13C, omega_1H, taus, pars->S2CHrs, tauf, pars->S2CHrf, pars->S2uf, D_CHr);
	R1CN = Dipolar_R1(omega_13C, omega_15N, taus, pars->S2CNs, tauf, pars->S2CNf, pars->S2uf, D_CN);
	R1CCAp = Dipolar_R1(omega_13C, omega_13C - wCOCa, taus, pars->S2CCAps, tauf, pars->S2CCApf, pars->S2uf, D_CCAp);
	R1CCAc = Dipolar_R1(omega_13C, omega_13C - wCOCa, taus, pars->S2CCAcs, tauf, pars->S2CCAcf, pars->S2uf, D_CCAc);

	return (Decimal) (R1CSA + R1CH + R1CHr + R1CN + R1CCAp + R1CCAc)*T_DOWN;
}

Decimal Calc_13CR2(struct Residue *res, struct Relaxation* relax, struct BCParameters *pars, struct Model *m) {
	/* Takes in residue and relaxation data, and outputs an R1 for given tau and S2. */
	Decimal field = relax->field * 1000000; // conversion to Hz
	Decimal w1 = relax->w1, wr = relax->wr;
	/* In the original MATLAB code the dipolar coupling constant was calculated on the fly.
	 * Here, because Planck's constant is 10^-34 (which would require a Decimal128, and
	 * software division) I've predefined it. Bond length taken as 1.02 A */
	unsigned int model = m->model;
	Decimal taus, tauf;
	if (pars->Eas != -1 || pars->Eaf != -1) {
		taus = temp_tau(pars->taus, pars->Eas, relax->T);
		tauf = temp_tau(pars->tauf, pars->Eaf, relax->T);
	} else {
		taus = pars->taus;
		tauf = pars->tauf;
	}

	Decimal omega_1H = 2 * M_PI * field;
	Decimal omega_13C, omega_15N, wCOCa;
	Decimal d2x, d2y, d2xy, d2tot;

	omega_13C = 2 * M_PI * field / 3.976489314034722;
	omega_15N = 2 * M_PI * field / 9.869683408806043;
	wCOCa = 120 * omega_13C * 0.000001;

	w1 *= T_DOWN;
	wr *= T_DOWN;
	omega_1H *= T_DOWN;
	omega_13C *= T_DOWN;
	omega_15N *= T_DOWN;
	wCOCa *= T_DOWN;

	Decimal *csa;
	csa = res->csaC;

	/* CSA relaxation contribution */
	Decimal R2CSAx, R2CSAy, R2CSAxy, R2CSA, R2CH, R2CHr, R2CN, R2CCAc, R2CCAp;
	Decimal J0sum = 0;
	if (model == MOD_GAF || model == MOD_GAFT || model == MOD_EGAF || model == MOD_EGAFT  || model == MOD_BGF || model == MOD_BGFT) {
		d2x = (Decimal) sq(((csa[2] - csa[0]) * 0.000001) * omega_13C * T_UP);
		d2y = (Decimal) sq(((csa[1] - csa[0]) * 0.000001) * omega_13C * T_UP);
		d2xy= (Decimal) sq(0.000001 * omega_13C * T_UP) * (csa[2] - csa[0]) * (csa[1] - csa[0]);

		R2CSAx = CSA_R2(omega_13C, w1, wr, taus, pars->S2CSAxs, tauf, pars->S2CSAxf, pars->S2uf, d2x, J0);
		R2CSAy = CSA_R2(omega_13C, w1, wr, taus, pars->S2CSAys, tauf, pars->S2CSAyf, pars->S2uf, d2y, J0);
		R2CSAxy = CSA_R2(omega_13C, w1, wr, taus, pars->S2CSAxys, tauf, pars->S2CSAxyf, pars->S2uf, d2xy, J0_CC);
		R2CSA = R2CSAx + R2CSAy + 2.*R2CSAxy;
	} else {
		J0sum += (2/3.) * J0(2 * M_PI * (w1 - 2 * wr), taus, pars->S2CSAxs, tauf, pars->S2CSAxf, pars->S2uf);
		J0sum += (2/3.) * J0(2 * M_PI * (w1 + 2 * wr), taus, pars->S2CSAxs, tauf, pars->S2CSAxf, pars->S2uf);
		J0sum += (4/3.) * J0(2 * M_PI * (w1 - wr), taus, pars->S2CSAxs, tauf, pars->S2CSAxf, pars->S2uf);
		J0sum += (4/3.) * J0(2 * M_PI * (w1 + wr), taus, pars->S2CSAxs, tauf, pars->S2CSAxf, pars->S2uf);

		d2tot= (sq(csa[2]) + sq(csa[1]) + sq(csa[0]));
		d2tot+=(-(csa[2] * csa[1] + csa[1] * csa[0] + csa[2] * csa[0]));
		d2tot*=sq(0.000001 * omega_13C * T_UP);
		R2CSA = (1/45.) * d2tot * (J0sum + 3 * J0(omega_13C, taus, pars->S2CSAxs, tauf, pars->S2CSAxf, pars->S2uf));

	}
	/* N Dipolar Interactions Contributions */
	R2CH = Dipolar_R2(omega_13C, omega_1H, w1, wr, taus, pars->S2CHs, tauf, pars->S2CHf, pars->S2uf, D_CH);
	R2CHr = Dipolar_R2(omega_13C, omega_1H, w1, wr, taus, pars->S2CHrs, tauf, pars->S2CHrf, pars->S2uf, D_CHr);
	R2CN = Dipolar_R2(omega_13C, omega_15N, w1, wr, taus, pars->S2CNs, tauf, pars->S2CNf, pars->S2uf, D_CN);
	R2CCAp = Dipolar_R2(omega_13C, omega_13C - wCOCa, w1, wr, taus, pars->S2CCAps, tauf, pars->S2CCApf, pars->S2uf, D_CCAp);
	R2CCAc = Dipolar_R2(omega_13C, omega_13C - wCOCa, w1, wr, taus, pars->S2CCAcs, tauf, pars->S2CCAcf, pars->S2uf, D_CCAc);

	return (Decimal) ((R2CSA + R2CH + R2CHr + R2CN + R2CCAp + R2CCAc)*(Decimal)T_DOWN);
}



/**
 * Calculates the order parameter relating to two orientations A and B undergoing anisotropic axial fluctuations of magnitude sig.
 * This function parallelises the calculation such that the loop over l, m, mp, k, kp only needs to occur once per motion - also
 * meaning that the calculation of the exponential only needs to be done once. Checks for conjugation are made to reduce the need
 * to calculate complex powers only to the cases where this is necessary.\n
 * Lienin1998, eq 6
 * @param sig
 *  Array containing [alpha, beta, gamma] deflection angles
 * @param A
 *  Array containing pointers to orientation vectors for 'A' for each S2 being calculated.
 * @param B
 *  Array containing pointers to orientation vectors for 'B' for each S2 being calculated.
 * @param S2
 *  Array containing pointers to S2 values - this is where output will be put.
 * @param length
 *  Number of S2 values being calculated (eg length of S2)
 * @param mode
 *  MODE_REAL or MODE_IMAG depending on which result is being returned.
 * @return Returns 1 in all cases.
 */
int GAF_S2(Decimal sig[3], struct Orient ** A, struct Orient ** B, Decimal * S2[], int length,  unsigned int mode) {
	int l, m, mp, k, kp, i;
	Complex * Amp = (Complex *) malloc(sizeof(Complex) * (long unsigned int) length);
	Complex temp, ttemp;


	for (i = 0; i < length; i++) {
		Amp[i] = 0;
	}



	Decimal sqsigs[3] = {sq(sig[0]), sq(sig[1]), sq(sig[2])};

	/* sq_s is an array where sq_s[i+2] = i^2 */
	int sq_s[5] = {4, 1, 0, 1, 4};

	/* sq_d is an array where sq_d[i+2][j+2] = (i^2 + j^2) / 2
	 * This speeds precalculates loop values */
	Decimal sq_d[5][5] = {
			{4.0, 2.5, 2.0, 2.5, 4.0},
			{2.5, 1.0, 0.5, 1.0, 2.5},
			{2.0, 0.5, 0.0, 0.5, 2.0},
			{2.5, 1.0, 0.5, 1.0, 2.5},
			{4.0, 2.5, 2.0, 2.5, 4.0}
	};

	Decimal lexp, kexp, mexp;

	/*
	 * The GAF order parameter equation is composed of four main parts. These are
	 *  1. Imaginary power of k-kp
	 *  2. Exponential of angular deflections with (m^2 + mp^2) etc
	 *  3. Product of Wigner matrices
	 *  4. Product of Spherical Harmonics
	 * As this is a nested sum it can be quite slow. In order to speed it up, a number of optimizations have been made.
	 *
	 *  1. If we are considering only real order parameters, then given spherical harmonics with
	 *	 no imaginary component it can be shown that any imaginary part of the i^(k-kp) is ignored.
	 *	 As (k-kp) is an integer, i^(k-kp) goes (0) 1, (1) i, (2) -1, (3) -i, (4) 1 and so on.
	 *	 So creal(i^(k)) is 0 for any odd k, and so these cases may be ignored.
	 *	 If k is divisible by 4, then this factor is 1, and if not, it is -1.
	 *	 This avoids the need for calling cpow() on each iteration, which is a hard operation.
	 *	 In the case that the spherical harmonics have an imaginary part, however, this must be done.
	 *  2. The component for each of alpha, beta, and gamma is done at their section in the loop. This means
	 *	 that the beta component is only calculated 5 times, the gamma component 125 times, and the alpha
	 *	 3125 times.
	 *	 The sigma values are presquared in sqsigs. While squaring is easy, doing this so many times can slow
	 *	 down the iteration. Additionally, the summed squares of the indices is precalculated (sq_s, sq_d).
	 *	 It has been considered to precalculate the sigma exponentials and then take these to the power of the
	 *	 sum of squares. This doesn't lead to a significant speed bump, but does make the code more complicated.
	 *  3. The product of Wigner matrices over k, kp, l, m, mp and exponentiation is done once and this is then
	 *	 used for all S2 parameters being calculated.
	 *  4. No optimization is done to the products of spherical harmonics.
	 */

	for (l = -2; l <= 2; l++) {
		lexp = -(sqsigs[1] * sq_s[l+2]);
		for (m = -2; m <= 2; m++) {
			for (mp = -2; mp <= 2; mp++) {
				mexp = -(sqsigs[2] * sq_d[m+2][mp+2]);
				for (k = -2; k <= 2; k++) {
					for (kp = -2; kp <= 2; kp++) {
						kexp = -(sqsigs[0] * sq_d[k+2][kp+2]);
						temp = 1;
						temp *= exp(lexp + kexp + mexp);
						temp *= Dwig[k+2][l+2] * Dwig[kp+2][l+2] * Dwig[m+2][k+2] * Dwig[mp+2][kp+2];
						for (i = 0; i < length; i++) {
							ttemp = temp;
							if (cimag(A[i]->Y2[m+2]) < 0.00001 && cimag(B[i]->Y2[mp+2]) < 0.00001 && mode==MODE_REAL) {
								/* In this case then Y2 and Y2c ar || model == MOD_BGF || model == MOD_BGFTe real numbers and therefore the last step of this is real.
								* If mode is MODE_REAL, then, we may safely ignore any case in which
								*  cpow(-1*I, k-kp)
								* is imaginary (as this will give only an imaginary component to the Amp).
								* This is the case for any (k-kp)%2 != 0.
								* In this case, if (k-kp)%4 == 0 the cpow function gives 1, else -1.
								*/
								if ((k-kp)%2 != 0)
									continue;
								if ((k - kp) % 4 == 0)
									ttemp *= 1;
								else
									ttemp *= -1.;
							} else {
								/* Other wise we have to go the long route... */
								Complex multip = cpow(-1 * I, k-kp);
								ttemp *= multip;
							}
							ttemp *= A[i]->Y2[m+2] * B[i]->Y2c[mp+2];
							Amp[i] += ttemp;
						}
					}
				}
			} 
		}
	}
	//Amp *= (4 * M_PI / 5.);
	for (i = 0;i < length; i++) {
		Amp[i] = Amp[i] * (4 * M_PI / 5.);
		switch (mode) {
			case MODE_REAL: *S2[i] = creal((Complex) Amp[i]); break;
			case MODE_IMAG: *S2[i] = cimag((Complex) Amp[i]); break;
			default: break;
		}
	}
	free(Amp);
	return 1;
}


/**
 * Calculates the order parameter for anisotropic model free motion. Essentially, the order_params define an ellipsoid.
 * The order parameter for a given theta, phi is the radius of the ellipsoid in that direction. eg;
 *   Let Sa, Sb, Sc be the squared order parameters along the following directions;
 *	 Sa - CaCa axis (eg the axis with gamma deflections about it)
 *	 Sb - ~N-H axis (eg the axis with alpha deflections about it)
 *	 Sc - perpendicular to plane axis (eg beta).
 *   Then the equation for an ellipsoid is;
 *	(sin(theta)^2 sin(phi)^2 / Sa^2 + sin(theta)^2 cos(phi)^2 / Sb^2 + cos(theta)^2 / Sc^2 = 1 / r^2
 *   Evaluating this for r for a given theta, phi gives us the amplitude of order parameter in that orientation.
 * @param order_params
 *  Array containing [Sa, Sb, Sc] squared order parameters
 * @param A
 *  Array containing pointers to orientation vectors for 'A' for each S2 being calculated.
 * @param S2
 *  Array containing pointers to S2 values - this is where output will be put.
 * @param length
 *  Number of S2 values being calculated (eg length of S2)
 * @return Returns 1 in all cases.
 */
int AIMF_S2(Decimal order_params[3], struct Orient ** A, Decimal * S2[], int length) {
	int i;
	Decimal ct, st, cp, sp;
	Decimal S2a = pow(order_params[0], 2.), S2b = pow(order_params[1], 2.), S2c = pow(order_params[2], 2.);
	Decimal SaSbSc = S2a * S2b * S2c;
	Decimal denom;
	for (i = 0; i < length; i++) {
		// A[i].theta, A[i].phi
		// order_params[0], [1], [2] (x, y, z)
		ct = cos(A[i]->theta);
		st = sin(A[i]->theta);
		cp = cos(A[i]->phi);
		sp = sin(A[i]->phi);

		denom = 0;
		denom += S2b * S2c * pow(st * sp, 2.);
		denom += S2a * S2c * pow(st * cp, 2.);
		denom += S2a * S2b * pow(ct, 2.);

		*S2[i] = sqrt(SaSbSc / denom);

	}
	return 1;
}



