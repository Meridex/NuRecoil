# The architecture of `NuRecoil` package


## CEvNS-like calculations

The output of such calculations is the event rate in a given energy bin (w/ or w/o detector effects). The formula is given as

$$
\frac{dR}{dE_{\rm det}} = N_T\int \frac{d\Phi}{dE_\nu}\frac{d\sigma}{dE_R}(E_\nu,E_R)f_{\rm res}(f_Q(E_R)E_R,E_{\rm det})dE_\nu dE_R
$$

- $\frac{d\Phi}{dE_\nu}$ is the neutrino flux
- $\frac{d\sigma}{dE_R}(E_\nu, E_R)$ is the differential cross section (wrt. the recoil energy) as function of the incoming neutrino energy and the recoil energy
- $f_Q(E_R)$ is the quenching factor at given recoil energy such that the actual energy that can be measured is given as $E_R^{ee}=f_Q\times E_R$
- $f_{\rm res}(f_QE_R,E_{\rm det})$ is the resolution function which is assumed to be Gaussian: $f_{\rm res} = \frac{1}{\sqrt{2\pi}\Delta E}e^{-\frac{(E_R-E_{\rm det})^2}{2\Delta E^2}}$ which has been correctly normalized. $\Delta E$ is the resolution depending on the experiments

## The neutrino flux

We should implement our own neutrino flux calculations, but can also rely on public available library, i.e. CONFLUX. In either way, we should keep the same API and provide the possibility for the user to use either implementations.

The neutrino flux will be given as function of the neutrino energy $E_\nu$ in $\rm MeV$ as well as the reactor thermal power $P$ ($\rm GW$), and the distance to the reactor core $L$ ($\rm m$) and given in the unit of $\#/{\rm MeV}/{\rm cm^2}$

- Phenomenological Fitting formula
- Interpolation with provided discrete data points
- Calling CONFLUX library

## The differential cross section

We will provide the CEvNS differential cross section in the SM, but also provide the possibility of any BSM modifications, i.e. different weak mixing angles, NSI, BSM neutrino EM properties etc. For this purpose, the differential cross section calculation should have uniform API for easier function call. The Standard CEvNS calculation only consider the nuclear recoil, but it should be possible to include electron recoil as well as the Migdal effect.

### SM CEvNS

$$
\frac{d\sigma}{dE_R} = \frac{G_F^2M}{\pi}\left(1- \frac{E_R}{E_\nu} + \frac{1}{2}\left(\frac{E_R}{E_\nu}\right)^2 - \frac{ME_R}{2E_\nu^2}\right)\left\vert Q_W(|\vec{q}|)\right\vert^2
$$

- $G_F$ is the Fermi constant
- $M$ is the target nucleus mass
- $E_\nu$ is the incoming neutrino energy
- $E_R$ is the recoil energy.
- $Q_W$ is the weak charge of the nucleus
- $|\vec{q}|\simeq \sqrt{2ME_R}$

The weak charge of the nucleus is given in the SM as
$$
Q_W(|\vec{q}|) = g_V^nNF_N(|\vec{q}|)+g_V^pZF_Z(|\vec{q}|)
$$
with $F_{N,Z}(|\vec{q}|)$ the neutron/proton form factors (representing the neutron/proton distributions within the corresponding nucleus) which have many different choices.

$g_V^{n,p}$ are the vector current interactions of the neutron and proton through mediator $Z$ and are given as
$$
\begin{align*}
g_V^n &= 2g_V^d+g_V^u \stackrel{\rm LO}{=} -\frac{1}{2}\\
g_V^p &= 2g_V^u+g_V^d \stackrel{\rm LO}{=} \frac{1}{2} - 2s_W^2
\end{align*}
$$
where the corresponding vector current couplings of the quarks are given by their quantum numbers
$$
g_V^f \stackrel{\rm LO}{=} T_3 - 2Qs_W^2
$$

### Other BSM Scenarios



## Quenching Factor

There are many different choices for the quenching factor. Linhard is the most commonly used one which depends on several phenomenological parameters. Other quenching factor can also be used.
