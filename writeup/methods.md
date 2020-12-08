---
bibliography: /home/ben/library.bib
title: "Methods."
author:
	- Ben P. Tatman
geometry: margin=2cm
fontsize: 12pt
csl: "/home/ben/nature.csl"
header-includes:
    - \usepackage{setspace}
    - \doublespacing
    - \usepackage{lineno}
    - \linenumbers
---

Dynamics Fitting
================

The various models were fit to the relaxation data using an in-house fitting model implemented using the C programming language. This model can be found at http://www.website.com. Expressions for the rates and spectral densities for the models can be found below. Minimization was performed using the Nelder-Mead algorithm[@Nelder.1965], using N repeats with random starting points. Table 1 shows the variation in the random starting points for each model, along with the bounds used in the fitting to ensure physical solutions were found[^1].

Table 1: Parameters used in the dynamic model fitting. SMF refers to Simple Model Free, EMF refers to Extended Model Free, GAF refers to Gaussian Axial Fluctuation Model, 3D and 6D refer to the dimensionality of the axial fluctuations. A prefix of T indicates temperature dependence. \* indicates these are bound via the non-temperature dependent forms.


| Parameter | Description | Units | Models | Start Range | Bound Range |
|-|-|-|-|-|-|
|$\tau$|Correlation time|s|SMF|$0-10^{-8}$|$\geq 0$|
|$\tau_0$|Correlation preexponential|s|SMFT|$0-10^{-15}$|$\geq 0$|
|$\tau^{\text{s}}$|Slow motion correlation time|s|EMF/GAF|$0-10^{-8}$|$10^{-5} \geq \tau \geq 0$|
|$\tau^{\text{f}}$|Fast motion correlation time|s|EMF/GAF|$0-10^{-11}$|$10^{-5} \geq \tau \geq 0$|
|$\tau^{\text{s}}_0$|Slow preexponential factor|s|EMFT/GAFT|$0-10^{-15}$|\*|
|$\tau^{\text{f}}_0$|Fast preexponential factor|s|EMFT/GAFT|$0-10^{-20}$|\*|
|$S^2$|Order parameter|\-|SMF|$0.5 - 1.0$|$1\geq S^2 \geq 0$|
|$S^2_{\text{s}}$|Slow order parameter|\-|EMF/EMFT|$S^{2}_{NH} - 1.0$|$1\geq S^2 \geq 0$|
|$S^2_{\text{f}}$|Fast order parameter|\-|EMF/EMFT/3D-GAF/3D-GAFT|$S^{2}_{NH} - 1.0$|$1\geq S^2 \geq 0$|
|$\sigma_{\alpha/\beta/\gamma}^{\text{s}}$|Angular deflections for slow motions|radians|3D/6D-GAF/GAFT|$0.00-0.25$|$0.52 \geq \sigma \geq 0.00$|
|$\sigma_{\alpha/\beta/\gamma}^{\text{f}}$|Angular deflections for fast motions|radians|6D-GAF/GAFT|$0.00-0.25$|$0.52 \geq \sigma \geq 0.00$|
|$\text{Ea}$|Correlation time activation energy|J mol$^{-1}$|SMFT|$0-60000$|\*|
|$\text{Ea}^{\text{s}}$|Slow correlation time activation energy|J mol$^{-1}$|EMFT/GAFT|$0-60000$|\*|
|$\text{Ea}^{\text{f}}$|Fast correlation time activation energy|J mol$^{-1}$|EMFT/GAFT|$0-60000$|\*|


The best fit parameters were found for all models by minimizing the $\chi^{2}$ target function:

$$\chi^{2} = \ \sum_{i}^{}{W_{i}\frac{{(X_{i,exp} - \ X_{i,\ calc})}^{2}}{\sigma_{i,\ exp}^{2}}}$$

where $X_{i}$ are relaxation rates and dipolar order parameter measurements and $\sigma_{i}$ appropriate experimental errors[^2]. $W_{i}$ is a weighting which is taken to be $1$ for all relaxation measurements and $100$ for dipolar order parameters to reflect the increased model reliance on these order parameters [@Lamley.2015]. The rigid limit NH distance was assumed to be N $\text{Ang}$.

Errors were estimated using Monte Carlo error analysis using N iterations. Relaxation rates were back-calculated from the model fit parameters, and random noise within experimental error was added using the Box-Muller method. The model was then fit to these new relaxation rates. This procedure was repeated N times per residue, and the error is given as twice the standard deviation of all of these repeats[^3].


Simple Model Free
-----------------

The Simple Model Free model was implemented as outlined in Lipari *et al*. 1982[@Lipari.1982]. In brief, the spectral density is taken as[^4]:

$$J\left( \omega\ (rad\ s^{- 1}) \right) = \ \frac{\left( 1 - \ S^{2} \right)\tau}{(1 + \left( \omega\tau \right)^{2})}\ $$

Where $\omega$ is the frequency of interest, and the other parameters are as defined in Table 1. The contribution to the $R_{1}$ relaxation rate from a dipolar interaction was obtained using[@Kurbanov.2011][^5]:

$$R_{1,\ dipolar}\left( \Omega,\omega \right) = \frac{1}{10}\omega_{D}^{2}(J\left( \omega - \ \Omega \right) + 3 \times J\left( \omega \right) + 6 \times J\left( \omega + \ \Omega \right))$$

Where $\Omega$ is the Larmor frequency for the nuclei being observed, $\omega$ is the Larmor frequency for the dipolar coupled nuclei, and the dipolar coupling is given by $\omega_{D}$. The CSA contribution is[^6]:

$$R_{1,\ CSA}\left( \Omega \right) = \ \frac{2}{15}\ \left( \delta_{11}^{2} + \ \delta_{22}^{2} + \ \delta_{33}^{2} - \ \delta_{11}\delta_{22} - \delta_{11}\delta_{33} - \ \delta_{22}\delta_{33} \right)\Omega^{2} \times J(\Omega)\ \ $$

In which $\delta_{11}^{2}$ are components of the CSA, parameterized using the isotropic chemical shift of the given nuclei. Defining $J_{k}$ as follows (where $\omega_{r}$ and $\omega_{1}$ are the spinning frequency and spin lock frequency in Hz, respectively)[^7].

$$J_{k}\left( \omega_{r},\omega_{1} \right) = \frac{2}{3}J(2\pi\left( \omega_{1} - 2\omega_{r} \right) + \ \frac{2}{3}J(2\pi\left( \omega_{1} + 2\omega_{r} \right) + \ \frac{4}{3}J(2\pi\left( \omega_{1} - \omega_{r} \right) + \ \frac{4}{3}J(2\pi\left( \omega_{1} + \omega_{r} \right)$$

The dipolar contribution to the $R_{1\rho}$ rate is taken as[^8];

$$R_{1\rho,\ dipolar}\left( \Omega,\omega \right) = \frac{1}{20}\omega_{D}^{2}(J_{k} + 3 \times J\left( \Omega \right) + \ J\left( \omega - \ \Omega \right) + 6 \times J\left( \omega \right) + 6 \times J\left( \omega + \ \Omega \right))$$

And correspondingly[^9]

$$R_{1\rho,\ CSA}\left( \Omega \right) = \ \frac{1}{45}\ \left( \delta_{11}^{2} + \ \delta_{22}^{2} + \ \delta_{33}^{2} - \ \delta_{11}\delta_{22} - \delta_{11}\delta_{33} - \ \delta_{22}\delta_{33} \right)\Omega^{2} \times (J_{k} + 3 \times J\left( \Omega \right))\ \ $$

For ^15^N, contributions from N-H, N-Hr, N-C and N-Ca dipolar interactions were considered. For ^13^C, C-H, C-Hr, C-N and C-C dipolar interactions were considered. The effective distances used for calculation of dipolar couplings are given in Table 2.


  **Bond**      **Length / Å**
  ------------- ---------------------
  N-H           1.02 [@Ferrage.2006]
  N-H (rest)    1.80 [@Ferrage.2006]
  N-C           1.33 [@Ferrage.2006]
  N-C$\alpha$   1.46 [@Ferrage.2006]
  C-H           2.04 [@Engh.1991]
  C-H (rest)    1.82 [@Ferrage.2006]
  C-N           1.33 [@Ferrage.2006]
  C-C           1.53 [@Engh.1991]

Extended Model Free
-------------------

The Extended Model Free spectral density was taken as [@Clore.1990][^10]:

$$J\left( \omega\ (rad\ s^{- 1}) \right)\  = \ \frac{\left( 1 - \ S_{\text{fast}}^{2} \right)\tau_{\text{fast}}}{(1 + \left( \omega\tau_{\text{fast}} \right)^{2})} + \ \frac{S_{\text{fast}}^{2}\left( 1 - S_{\text{slow}}^{2} \right)\tau_{\text{slow}}}{\left( 1 + \left( \omega\tau_{\text{slow}} \right)^{2} \right)}\ $$

The form of the relaxation rates was as in Simple Model Free.


Gaussian Axial Fluctuations
---------------------------

The GAF spectral density was taken as[^11]

$$J\left( \omega \right) = \ \frac{\left( 1 - \ S_{\text{fast}}^{2} \right)\tau_{\text{fast}}}{(1 + \left( \omega\tau_{\text{fast}} \right)^{2})} + \ \frac{1}{P_{2}(\cos\left( \theta_{\text{\mu
u}} \right))}\frac{S_{\text{fast}}^{2}\left( P_{2}(\cos\left( \theta_{\text{\mu
u}} \right)) - S_{\text{slow}}^{2} \right)\tau_{\text{slow}}}{\left( 1 + \left( \omega\tau_{\text{slow}} \right)^{2} \right)}\ $$

Where $P_{2}\left( \cos\left( \theta_{\text{\mu
u}} \right) \right)$ equals 1 for autocorrelated motions and -1/2 for cross correlated motions. The GAF order parameters were obtained as by Lienin et al, 1998[@Lienin.1998] via[^12]


$$S_{\text{\mu
u}}^{2} = \ \frac{4\pi}{5}\sum_{l,\ k,k^{'},m,m^{'} = - 2}^{2}{\left( - i \right)^{k - k^{'}}\exp\left\lbrack - \frac{\sigma_{\alpha}^{2}\left( k^{2} + k^{'2} \right)}{2} - \ \sigma_{\beta}^{2}l^{2} - \ \frac{\sigma_{\gamma}^{2}\left( m^{2} + m^{'2} \right)}{2} \right\rbrack \times \ d_{\text{kl}}^{\left( 2 \right)}\left( \frac{\pi}{2} \right)d_{k^{'}l}^{\left( 2 \right)}\left( \frac{\pi}{2} \right)d_{\text{mk}}^{\left( 2 \right)}\left( \frac{\pi}{2} \right)d_{m^{'}k^{'}}^{\left( 2 \right)}\left( \frac{\pi}{2} \right)Y_{2m}\left( \mathbf{e}_{\mu}^{\text{pp}} \right)Y_{2m^{'}}^{*}(\mathbf{e}_{\nu}^{\text{pp}})}$$

In which $d_{mm^{'}}^{\left( 2 \right)}\left( \frac{\pi}{2} \right)$ are reduced Wigner matrix elements evaluated at $\frac{\pi}{2}$, $Y_{2m}\left( \theta(t),\ \phi(t) \right)$ are second order spherical harmonics, and $\mathbf{e}_{\alpha}^{\text{pp}}$, $\mathbf{e}_{\beta}^{\text{pp}}$, $\mathbf{e}_{\gamma}^{\text{pp}}$ are principal axes rigidly attached to the peptide plane. The deflection angles, $\sigma_{\alpha}$, $\sigma_{\beta}$, $\sigma_{\gamma}$ refer to Gaussian rotations about these axes. Here, we have taken the same definition of these principal axes as in Lienin et al, 1998
[@Lienin.1998]. Relaxation rates were obtained. Dipolar contributions to the relaxation rates were taken as in EMF, while CSA contributions were calculated as [@Bremi.1997][^13]:


$$R_{1,\ CSA}^{X}\left( \Omega \right) = \ \frac{1}{15} \times \left( \delta_{\text{xx}} - \delta_{\text{zz}} \right)^{2} \times \Omega^{2} \times J(\Omega)\ $$

$$R_{1,\ CSA}^{Y}\left( \Omega \right) = \ \frac{1}{15} \times \left( \delta_{\text{yy}} - \delta_{\text{zz}} \right)^{2} \times \Omega^{2} \times J(\Omega)$$

$$R_{1,\ CSA}^{\text{XY}}\left( \Omega \right) = \ \frac{1}{15} \times (\delta_{\text{xx}} - \ \delta_{\text{zz}})(\delta_{\text{yy}} - \delta_{\text{zz}}) \times \Omega^{2} \times J(\Omega)$$

$$R_{1,\ CSA}^{}\left( \Omega \right) = \ R_{1,\ CSA}^{X}\left( \Omega \right) + R_{1,\ CSA}^{Y}\left( \Omega \right) + 2 \times R_{1,\ CSA}^{\text{XY}}\left( \Omega \right)$$

For $R_{1\rho}$, CSA contributions to the relaxation rate were taken as[^14]:

$$R_{1\rho,\ CSA}^{X}\left( \Omega \right) = \ \frac{1}{15} \times \left( \delta_{\text{xx}} - \delta_{\text{zz}} \right)^{2} \times \Omega^{2} \times (J_{k} + 3 \times J\left( \Omega \right))\ $$

$$R_{1\rho,\ CSA}^{Y}\left( \Omega \right) = \ \frac{1}{15} \times \left( \delta_{\text{yy}} - \delta_{\text{zz}} \right)^{2} \times \Omega^{2} \times (J_{k} + 3 \times J\left( \Omega \right))$$

$$R_{1\rho,\ CSA}^{\text{XY}}\left( \Omega \right) = \ \frac{1}{15} \times (\delta_{\text{xx}} - \ \delta_{\text{zz}})(\delta_{\text{yy}} - \delta_{\text{zz}}) \times \Omega^{2} \times (J_{k} + 3 \times J\left( \Omega \right))$$

$$R_{1\rho,\ CSA}^{}\left( \Omega \right) = \ R_{1\rho,\ CSA}^{X}\left( \Omega \right) + R_{1\rho,\ CSA}^{Y}\left( \Omega \right) + 2 \times R_{1\rho,\ CSA}^{\text{XY}}\left( \Omega \right)$$


Temperature Dependence
----------------------

For the temperature dependent models we assumed that the timescales of motion would be time dependent according to the Arrhenius equation[^15]:

$$\tau\left( T \right) = \ \tau^{0}\exp\left( \frac{\text{Ea}}{R T} \right)$$


Fitted Correlation Functions
----------------------------

Using the parameters obtained from the dynamics model fitting described above we determined the form of the correlation functions for these motions. For Simple Model Free analysis, the correlation functions were plotted according to [@Lipari.1982][^16]:

$$C\left( t \right) = S^{2} + \left( 1 - S^{2} \right)\exp\left( \frac{- t}{\tau} \right)$$

For Extended Model Free, the following form was used [@Clore.1990][^17]:

$$C\left( t \right) = S^{2} + \left( 1 - S_{f}^{2} \right)\exp{\left( \frac{- t}{\tau_{f}} \right) + S_{f}^{2}\left( 1 - S_{s}^{2} \right)\exp\left( \frac{- t}{\tau_{s}} \right)}$$


AMBER MD
========

A molecular dynamics trajectory for GB1 in a supercell was computed using AMBER MD[@Case.2005]. The coordinates of the X-ray structure of GB1 (PDB:2gi9) were taken as a starting conformation. A unit cell containing 4 GB1 molecules was created using the AMBER MD *UnitCell* utility, which was propagated by *PropPDB* to produce a 3x3x3 supercell containing 108 GB1 proteins. To this, 108 PO~4~^3-^ counter ions were added using *AddToBox*, along with 108 MPD and 648 IPA cocrystallising ligands. Finally, 12852 water molecules were added to the box. This was then charge balanced by addition of sodium ions, to give an overall box size of 75.591 Å x 107.152 Å x 150.822 Å. The ff14SB [@Maier.2015] forcefield was used for the GB1 proteins, while TIP3P was used for water and GAFF for the cocrystals.

The system was heated to 300 K and allowed to equilibrate for 25 ps with a timestep of 0.5 fs. The system was then simulated for a full 1 $\mu$s run with a 2 fs timestep and a cut off of 11 Å for non-bonded interactions. The temperature was maintained at 300 K using Langevin thermostat, and the SHAKE algorithm[@Ryckaert.1977] was applied to all bond lengths involving a hydrogen atom. Anisotropic pressure scaling was used with periodic boundary conditions.

Correlation functions for NH, CH and CN bonds were extracted from the final trajectories for each of the 108 GB1 molecules simulated in the following manner. *Pytraj*[@Nguyen.2016] was used align the supercells to minimize the effects of overall tumbling. Atomic coordinates were then extracted, and for each time step the orientation of the relevant atomic bond vector was obtained. The correlation function was then calculated using [@Bremi.1997]


$$C_{\mu\nu}^{\text{int}}\left( k\Delta t \right) = \ \frac{1}{(N - k)}\sum_{i = 1}^{N - k}\frac{3\left( \mathbf{e}_{\mu,i + k} \cdot \mathbf{e}_{\nu,i} \right)^{2} - 1}{2}$$

Bibliography
============



[^1]: (see chisq.c lines 34-211, main.c 135-204. I've omitted orientation variation for now.). crosen.c includes Nelder-Mead method.

[^2]: See chisq.c lines 246-248 for relaxation data fitting, and lines 254-288 for order parameter fitting

[^3]: Errors.c lines 26-35 for the Box-Muller method, 48-60 for standard deviation calculation. Errors calculated in calc\_errors(). Back calculation of rates is line 131, adding error is done 132. Simplex performed 141. Statistics and output done main.c lines 457-460.

[^4]: Smf.c, lines 15-18

[^5]: Smf.c, SMF\_Dipolar\_R1 33-40

[^6]: Smf.c eg 123-125

[^7]: This is J0sum, eg 177-181 in smf.c and 203-207 in emf.c

[^8]: Smf.c 58-65, emf.c 81-88.

[^9]: Smf.c 205-207 and emf.c 231-233.

[^10]: Emf.c 14-24.

[^11]: The case for autocorrelated is just J0\_EMF, emf.c 14-24. The cross correlate is emf.c 30-40.

[^12]: Gaf.c 103-177. Note I've messed it around a bit. The individual terms of the exponential are calculated once per loop, and I've considered that the Y2m components have no real components and so if we're selecting for real order parameters we can ignore any cases where (-i)\^(k-k') is imaginary.

[^13]: Gaf.c, 248-256, 416-422.

[^14]: Gaf.c, 36-60 implements each summation, then 333-336 and 501-503.

[^15]: eg chisq.c 46, 82-83, 128-129, 172-173.

[^16]: Correlation.c, 171-172

[^17]: Correlation.c, 134-136


