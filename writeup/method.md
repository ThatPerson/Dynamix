# Methods

Fitting of relaxation data to the above models was performed using an in house program developed in C. Expressions for the spectral densities and rates are detailed below. Model fitting was performed using the Nelder-Mead method[cite], starting from random starting points, repeated N times. The Box-Muller transform was used to obtain normally distributed random numbers. The best-fit parameters for each model was obtained by minimizing the $\chi^{2}$ target function:

$$ \chi^{2} = \sum_i \frac{(\chi_{i, calc} - \chi_{i, exp})^{2}}{\sigma_{i, exp}^{2}} $$

In which $\chi_{i}$ are relaxation rates and dipolar coupling measurements and $\sigma_i$ are experimental errors. Errors were fitted as in previous publications[@Lamley2015a]; in brief, errors were calculated using Monte Carlo error analysis with N iterations. Random noise was added to relaxation rates which were then refit to the model, with the error being taken as twice the standard deviation of these model fits.

## Simple Model Free

Simple model free was implemented as previously[@Lipari1982,@Lamley2015a]. The spectral density function was taken as:

$$ J_{\text{SMF}}(\omega) = \frac{(1 - S^{2}) \tau_{\text{eff}}}{1 + (\omega \tau_{\text{eff}})^2} $$

Where $S^{2}$ is the squared order parameter and $\tau_{\text{eff}}$ the effective time scale of motion.

Expressions for the relaxation rates are given below (in which $X = ^{15}N/ ^{13}C$ depending on nuclei observed, and $Y = ^{1}H$):

$$ R_{1, dipolar, {X}Y} = \frac{1}{10} \left(\frac{\mu_0}{2\pi} \frac{\gamma_{X} \gamma_{Y}}{\hbar r_{XY}^{3}}\right)^{2} (J_0(\omega_Y - \omega_X) + 3 J_1(\omega_X) + 6 J_2(\omega_Y + \omega_X)) $$

$$ R_{1, CSA} = \frac{2}{15} \omega_{X}^{2} (\sigma_{11}^2 + \sigma_{22}^{2} + \sigma_{33}^{2} - \sigma_{11}\sigma_{22} - \sigma_{11}\sigma_{33} - \sigma_{22}\sigma_{33}) J_1(\omega_X) $$

$$ R_{1\rho, dipolar, XY} = \frac{1}{20} \left(\frac{\mu_0}{2\pi} \frac{\gamma_{X} \gamma_{Y}}{\hbar r_{XY}^{3}}\right)^{2} (4 J_{0}(\omega_1) + 3 J_{1}(\omega_X) + J_{0}(\omega_Y - \omega_X) + 6 J_{1}(\omega_Y) + 6 J_2(\omega_{Y} + \omega_{X})) $$

$$ R_{1\rho, CSA} = \frac{1}{45} \omega_{X}^{2} (\sigma_{11}^2 + \sigma_{22}^{2} + \sigma_{33}^{2} - \sigma_{11}\sigma_{22} - \sigma_{11}\sigma_{33} - \sigma_{22}\sigma_{33}) (\frac{1}{3}( 2 J_0(2 \pi (\omega_1 - 2 \omega_r)) + 2 J_0(2 \pi (\omega_1 + 2 \omega_r)) + 4 J_0(2 \pi (\omega_1 - \omega_r)) + 4 J_0(2 \pi (\omega_1 + \omega_r))) + 3 J_{0}(\omega_X)))$$

In which $\sigma_{11} > \sigma_{22} > \sigma_{33}$ being the CSA components for $^{15}N$ or $^{13}C$. 

## Extended Model Free

The spectral density function for the extended model free model was[@Clore1990]:

$$ J_{\text{EMF}}(\omega) = \frac{(1 - S^{2}_{\text{fast}}) \tau_{\text{fast}}}{1 + (\omega \tau_{\text{fast}})^{2}} + \frac{S^{2}_{\text{fast}} (1 - S^{2}_{\text{slow}}) \tau_{\text{slow}}}{1 + (\omega \tau_{\text{slow}})^{2}} $$

Expressions for relaxation rates were taken as in simple model free.

## Gaussian Axial Fluctuations

The spectral density function for the models including gaussian axial fluctuations was taken as:

$$ J_{\mu\nu, \text{GAF}}(\omega) = \frac{(1 - S^{2}_{\text{fast}}) \tau_{\text{fast}}}{1 + (\omega \tau_{\text{fast}})^{2}} + \frac{S^{2}_{\text{fast}} (P_{2}(\cos(\theta_{\mu\nu})) - S^{2}_{\text{slow}}) \tau_{\text{slow}}}{P_{2}(\cos(\theta_{\mu\nu})) ( 1 + (\omega \tau_{\text{slow}})^{2})} $$

Where $\nu$ and $\nu$ refer to the dipolar vector of pricipal axis of the CSA tensor, and $P_2(\cos(\theta_{\mu\nu}))$ is the second order Legendre polynomial with $\theta_{\mu\nu}$ being the angle between principal axes $\mu$ and $\nu$. This is equal to $1$ for autocorrelated processes (giving the EMF form of the spectral density), and $-\frac{1}{2}$ for cross correlated processes. 

The order parameters under gaussian axial fluctuations were given by the following expression from Lienin 1998[@Lienin1998].

$$ S^{2}_{\mu\nu} = \frac{4\pi}{5} \sum_{l,k,k',m,m' = -2}^{2} (-i)^{k-k'} \exp(-\frac{\sigma_{\alpha}^{2} (k^2 + k'^2)}{2} - \sigma_{\beta}^{2} l^2 - \frac{\sigma_{\gamma}^{2} (m^{2} + m'^{2})}{2}) \times d_{kl}^{(2)}(\pi/2) d_{k'l}^{(2)}(\pi/2) d_{mk}^{(2)}(\pi/2) d_{m'k'}^{(2)}(\pi/2) Y_{2m}(e^{pp}_{\mu}) Y^*_{2m')(e^{pp}_{\nu}) $$

In this work we have taken the angle pairs $e^{pp}_{\mu}$ as in Lienin et al.[@Lienin1998]. In order to validate this approach we have also utilised the model allowing the axes attached to the peptide plane to rotate (not shown), however this was found to have no significant effect.

## Temperature Dependence

We model temperature dependence using the Arrhenius equation on the timescales.

$$ \tau = \tau_0 \exp(E_{a} / RT) $$

