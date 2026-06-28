# About the Quenching Factor

For nuclear recoil signals, the quenching factor maps the true nuclear recoil
energy $E_R$ to an electron-equivalent detector energy:

$$
	\begin{align}
	E_{\mathrm{ee}} = Q(E_R)E_R
	\end{align}
$$

where $Q(E_R)$ is dimensionless.  In semiconductor ionization detectors $Q$ is
usually the ionization efficiency.  In scintillators $Q$ is usually the light-yield
ratio for nuclear recoils relative to electron recoils.  In noble-liquid TPCs the
response is better described by separate light and charge yields, so a single
scalar $Q(E_R)$ is usually only an approximation.

The input energy unit used in `NuRecoil` is MeV, but most quenching-factor
papers write formulas in keV.  In the formulas below, unless explicitly stated,
$E_R$ is in $\mathrm{keV_{nr}}$.


## Standard Lindhard Quenching Factor

The standard implementation-ready Lindhard formula is

$$
	\begin{align}
	Q_{\mathrm{L}}(E_R; k) &= \frac{k g(\epsilon)}{1+k g(\epsilon)} \\
	g(\epsilon) &= 3\epsilon^{0.15}+0.7\epsilon^{0.6}+\epsilon \\
	\epsilon &= 11.5\left(\frac{E_R}{\mathrm{keV}}\right)Z^{-7/3}.
	\end{align}
$$

The conventional semi-empirical value is

$$
	\begin{align}
	k_{\mathrm{L}} = 0.133 Z^{2/3}A^{-1/2}.
	\end{align}
$$

For Ge, this gives $k_{\mathrm{L}}\simeq 0.157$.  Many analyses instead fit
$k$ from calibration data.  Useful values are:

| Material | Suggested $k$ | Source | Notes |
| --- | ---: | --- | --- |
| Ge | 0.157 | Lindhard semi-empirical formula | $Z=32$, $A\simeq72.6$ |
| Ge | $0.1789^{+0.0014}_{-0.0010}$ | Scholz et al. 2016 | with adiabatic correction |
| Ge | $0.162\pm0.004$ | Bonhomme et al. 2022 | no adiabatic correction required above $\sim0.6\,\mathrm{keV_{nr}}$ |
| Si | $0.146$ | Lindhard semi-empirical formula | $Z=14$, $A\simeq28.1$ |

Implementation notes:

- Convert the package energy from MeV to keV before evaluating $\epsilon$.
- Return $Q=0$ below the detector-specific threshold only if the detector model
  requires such a cutoff; the Lindhard formula itself is smooth down to zero.
- For isotope-specific work use the isotope mass number $A$; otherwise natural
  abundance average $A$ is sufficient at the precision of this model.

References:

- J. Lindhard, V. Nielsen, M. Scharff, and P. V. Thomsen, "Integral equations
  governing radiation effects", Mat. Fys. Medd. Dan. Vid. Selsk. 33, no. 10
  (1963).
- J. Lindhard, M. Scharff, and H. E. Schiott, "Range concepts and heavy ion
  ranges", Mat. Fys. Medd. Dan. Vid. Selsk. 33, no. 14 (1963).
- J. D. Lewin and P. F. Smith, "Review of mathematics, numerical factors, and
  corrections for dark matter experiments based on elastic nuclear recoil",
  Astropart. Phys. 6, 87 (1996).
- B. J. Scholz et al., "Measurement of the low-energy quenching factor in
  germanium using an $^{88}$Y/Be photoneutron source", Phys. Rev. D 94, 122003
  (2016), [arXiv:1608.03588](https://arxiv.org/abs/1608.03588).
- A. Bonhomme et al., "Direct measurement of the ionization quenching factor of
  nuclear recoils in germanium in the keV energy range", Eur. Phys. J. C 82, 815
  (2022), [arXiv:2202.03754](https://arxiv.org/abs/2202.03754).


## Modified Lindhard with an Adiabatic Correction

Scholz et al. use a multiplicative adiabatic correction factor to suppress the
Lindhard ionization response for slow nuclear recoils:

$$
	\begin{align}
	Q_{\mathrm{Scholz}}(E_R;k,\xi)
		&= Q_{\mathrm{L}}(E_R;k) F_{\mathrm{AC}}(E_R;\xi), \\
	F_{\mathrm{AC}}(E_R;\xi)
		&= 1-\exp\left[-\frac{E_R}{\xi}\right].
	\end{align}
$$

Here $E_R$ and $\xi$ are both in $\mathrm{keV_{nr}}$.  Their best-fit Ge
parameters at $\sim77\,\mathrm{K}$ are

$$
	\begin{align}
	k &= 0.1789^{+0.0014}_{-0.0010}, \\
	\xi &= 0.16^{+0.10}_{-0.13}\,\mathrm{keV_{nr}}.
	\end{align}
$$

The corresponding recoil-to-electron-equivalent conversion for a multi-scatter
event in the Scholz analysis is

$$
	\begin{align}
	E_{\mathrm{ee}}
		= \sum_i E_{R,i} Q_{\mathrm{L}}(E_{R,i};k)
		  \eta(x_i;\delta,\tau)F_{\mathrm{AC}}(E_{R,i};\xi),
	\end{align}
$$

where $\eta(x_i;\delta,\tau)$ is a detector charge-collection correction.  For a
generic quenching model in `NuRecoil`, implement only the product
$Q_{\mathrm{L}}F_{\mathrm{AC}}$ unless a detector-specific dead-layer model is
also being implemented.

Suggested model class:

- `AdiabaticLindhardQuenching(k=0.1789, xi_keV=0.16)` for Ge.
- Set `xi_keV=0` or disable this factor to recover ordinary Lindhard.

Reference:

- B. J. Scholz et al., Phys. Rev. D 94, 122003 (2016),
  [arXiv:1608.03588](https://arxiv.org/abs/1608.03588), Eq. (5), Table I.


## Lindhard with Binding-Energy Cutoff

Sarkis, Aguilar-Arevalo, and D'Olivo reintroduce an effective atomic binding
energy into the Lindhard integral-equation treatment.  Define the reduced
variables

$$
	\begin{align}
	\epsilon_R &= c_Z E_R, \\
	u &= c_Z U, \\
	c_Z &= 11.5 Z^{-7/3}\,\mathrm{keV}^{-1},
	\end{align}
$$

where $U$ is an effective binding-energy parameter in keV.  Their quenching
factor is

$$
	\begin{align}
	f_n = \frac{\bar{\eta}}{\epsilon_R}
		= \frac{\epsilon+u-\bar{\nu}(\epsilon)}{\epsilon+u}.
	\end{align}
$$

In the constant-$u$ model, the QF vanishes when

$$
	\begin{align}
	\epsilon_R \le 2u
	\quad\Longleftrightarrow\quad
	E_R \le 2U.
	\end{align}
$$

### Semi-Analytic Ansatz

A simple implementation path is their ansatz

$$
	\begin{align}
	\bar{\nu}(\epsilon)
		&= \bar{\nu}_{\mathrm{L}}(\epsilon)+C_0\epsilon^{1/2}+C_1+u, \\
	\bar{\nu}_{\mathrm{L}}(\epsilon)
		&= \frac{\epsilon}{1+k g(\epsilon)}.
	\end{align}
$$

The QF is then

$$
	\begin{align}
	Q_{\mathrm{BE}}(E_R)
		= \max\left[
			0,\,
			\frac{\epsilon_R-\bar{\nu}(\epsilon)}{\epsilon_R}
		\right],
		\qquad
		\epsilon = \epsilon_R-u.
	\end{align}
$$

Use this ansatz only for $\epsilon\ge u$; set $Q_{\mathrm{BE}}=0$ for
$\epsilon_R\le2u$.  The fit parameters reported by Sarkis et al. are:

| Material | $C_0$ | $C_1$ | $U$ [keV] | Notes |
| --- | ---: | ---: | ---: | --- |
| Si | $(9.1\pm4.4)\times10^{-3}$ | $(3.33\pm1.2)\times10^{-5}$ | $0.15\pm0.06$ | cutoff near $300\,\mathrm{eV}$ |
| Ge | $(3.0\pm1.3)\times10^{-4}$ | $(0.62\pm0.12)\times10^{-5}$ | $0.02\pm0.01$ | cutoff near $40\,\mathrm{eV}$ |

### Numerical Integral-Equation Solution

The preferred Sarkis et al. model solves

$$
	\begin{align}
	k\epsilon^{1/2}\bar{\nu}'(\epsilon)
	-\frac{1}{2}k\epsilon^{3/2}\bar{\nu}''(\epsilon)
	=
	\int_{\epsilon u}^{\epsilon^2}
	\frac{f(t^{1/2})}{2t^{3/2}}
	\left[
		\bar{\nu}\left(\epsilon-\frac{t}{\epsilon}\right)
		+\bar{\nu}\left(\frac{t}{\epsilon}-u\right)
		-\bar{\nu}(\epsilon)
	\right]dt.
	\end{align}
$$

The boundary prescription is

$$
	\begin{align}
	\bar{\nu}(\epsilon) &=
	\begin{cases}
	\epsilon+u, & \epsilon<u, \\
	\epsilon+u-\lambda(\epsilon), & \epsilon\ge u,
	\end{cases} \\
	\lambda(u)&=0,\\
	\alpha_1 &\equiv \lambda'(u^+)=1+\frac{1}{2u}\alpha_2,\\
	-\frac{2}{u} &\le \alpha_2 \equiv \lambda''(u^+) \le 0.
	\end{align}
$$

The physical solution must satisfy

$$
	\begin{align}
	0\le\bar{\nu}'(\epsilon)\le1,
	\qquad
	\lim_{\epsilon\to\infty}\bar{\nu}''(\epsilon)=0^-.
	\end{align}
$$

Fit parameters for the numerical solution:

| Material | $k$ | $U$ [keV] | Suggested validity |
| --- | ---: | ---: | --- |
| Si | $0.161^{+0.029}_{-0.020}$ | $0.15^{+0.10}_{-0.05}$ | use above $\sim0.5\,\mathrm{keV}$ |
| Ge | $0.162^{+0.017}_{-0.024}$ | $0.02^{+0.015}_{-0.010}$ | use above $\sim50\,\mathrm{eV}$ |

The paper reports supplemental tables for the numerical solutions $f_n(E_R)$
for Si and Ge.  If this model is implemented, the most stable approach is to
ship those tables as data and interpolate them, rather than solving the
integro-differential equation at runtime.

Binding-energy scale data from the same paper:

| Material | Shell / process | Energy [eV] |
| --- | --- | ---: |
| Si | 2p shell above [Ne] core | 100 |
| Si | average electron-hole creation | 3.7 |
| Si | dislocation / Frenkel-pair creation | 36 |
| Ge | 3d shell above [Ar] core | 30 |
| Ge | average electron-hole creation | 3.0 |
| Ge | dislocation / Frenkel-pair creation | 23 |

Reference:

- Y. Sarkis, A. Aguilar-Arevalo, and J. C. D'Olivo, "Study of the ionization
  efficiency for nuclear recoils in pure crystals", Phys. Rev. D 101, 102001
  (2020), [arXiv:2001.06503](https://arxiv.org/abs/2001.06503), Eqs. (1), (2),
  (9), (11), (13), (15)-(20), Tables I, III, IV.


## Experimental Data and Interpolation

For a concrete detector, the most robust implementation is a tabulated QF:

$$
	\begin{align}
	Q(E_R)=\mathrm{Interp}\left[\{E_i,Q_i,\sigma_i\}\right].
	\end{align}
$$

Recommended implementation choices:

- Use linear interpolation by default.
- Reject extrapolation by default, or clamp to the nearest endpoint only when the
  caller explicitly asks for `extrapolate="clamp"`.
- Keep statistical and systematic errors in the data file, even if the first code
  implementation only interpolates central values.
- For repeated measurements at nearly the same $E_R$, either keep all points or
  compute an inverse-variance weighted average:

$$
	\begin{align}
	\bar{Q} =
	\frac{\sum_i Q_i/\sigma_i^2}{\sum_i 1/\sigma_i^2},
	\qquad
	\sigma_{\bar{Q}} =
	\left(\sum_i 1/\sigma_i^2\right)^{-1/2}.
	\end{align}
$$

### Ge Ionization QF Data from Bonhomme et al.

The following values are directly transcribed from Bonhomme et al. 2022, Tables
5 and 6.  They are per-point constant-Q fits to coincidence distributions.  The
column `dq` includes the geometrical uncertainty quoted in the table; correlated
systematic uncertainties are discussed separately in their Table 4.

| $E_R$ [keVnr] | $\delta E_R$ | $Q$ | $\delta Q_{\rm stat}$ | $\delta Q$ | Beam [keV] |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.97 | 0.05 | 0.181 | 0.010 | 0.013 | 250 |
| 1.17 | 0.07 | 0.192 | 0.015 | 0.019 | 250 |
| 1.17 | 0.07 | 0.182 | 0.011 | 0.016 | 250 |
| 1.54 | 0.07 | 0.190 | 0.007 | 0.011 | 250 |
| 1.54 | 0.07 | 0.183 | 0.010 | 0.013 | 250 |
| 1.58 | 0.07 | 0.193 | 0.007 | 0.010 | 250 |
| 1.58 | 0.07 | 0.191 | 0.008 | 0.011 | 250 |
| 1.92 | 0.04 | 0.190 | 0.020 | 0.020 | 250 |
| 1.92 | 0.04 | 0.182 | 0.005 | 0.006 | 250 |
| 1.99 | 0.04 | 0.187 | 0.004 | 0.006 | 250 |
| 1.99 | 0.04 | 0.193 | 0.006 | 0.007 | 250 |
| 0.79 | 0.11 | 0.171 | 0.008 | 0.026 | 640 |
| 0.80 | 0.11 | 0.178 | 0.008 | 0.027 | 640 |
| 0.81 | 0.12 | 0.148 | 0.012 | 0.027 | 640 |
| 1.05 | 0.13 | 0.175 | 0.002 | 0.023 | 640 |
| 1.06 | 0.13 | 0.187 | 0.005 | 0.025 | 640 |
| 1.09 | 0.15 | 0.165 | 0.000 | 0.024 | 640 |
| 1.60 | 0.13 | 0.172 | 0.001 | 0.015 | 640 |
| 1.64 | 0.13 | 0.164 | 0.002 | 0.014 | 640 |
| 2.49 | 0.11 | 0.185 | 0.001 | 0.009 | 640 |
| 2.49 | 0.11 | 0.189 | 0.001 | 0.009 | 640 |
| 3.02 | 0.18 | 0.194 | 0.001 | 0.012 | 640 |
| 3.02 | 0.18 | 0.196 | 0.001 | 0.012 | 640 |
| 3.99 | 0.17 | 0.203 | 0.001 | 0.009 | 640 |
| 3.99 | 0.17 | 0.206 | 0.001 | 0.009 | 640 |
| 4.08 | 0.17 | 0.205 | 0.001 | 0.009 | 640 |
| 4.08 | 0.17 | 0.207 | 0.001 | 0.009 | 640 |
| 4.97 | 0.11 | 0.218 | 0.001 | 0.005 | 640 |
| 4.97 | 0.11 | 0.219 | 0.001 | 0.005 | 640 |
| 5.15 | 0.10 | 0.215 | 0.001 | 0.005 | 640 |
| 5.15 | 0.10 | 0.216 | 0.001 | 0.005 | 640 |
| 0.60 | 0.08 | 0.159 | 0.021 | 0.031 | 500 |
| 0.61 | 0.09 | 0.166 | 0.006 | 0.024 | 500 |
| 0.61 | 0.10 | 0.159 | 0.013 | 0.028 | 500 |
| 0.80 | 0.10 | 0.182 | 0.015 | 0.028 | 500 |
| 0.80 | 0.10 | 0.177 | 0.007 | 0.024 | 500 |
| 0.82 | 0.12 | 0.159 | 0.011 | 0.025 | 500 |
| 1.22 | 0.10 | 0.175 | 0.002 | 0.015 | 500 |
| 1.25 | 0.11 | 0.167 | 0.001 | 0.014 | 500 |
| 1.73 | 0.09 | 0.187 | 0.002 | 0.010 | 500 |
| 1.73 | 0.09 | 0.190 | 0.002 | 0.010 | 500 |
| 1.90 | 0.09 | 0.186 | 0.002 | 0.009 | 500 |
| 2.30 | 0.14 | 0.195 | 0.002 | 0.012 | 500 |
| 2.35 | 0.13 | 0.202 | 0.001 | 0.012 | 500 |
| 2.35 | 0.13 | 0.203 | 0.002 | 0.012 | 500 |
| 3.04 | 0.13 | 0.199 | 0.001 | 0.009 | 500 |
| 3.04 | 0.13 | 0.200 | 0.001 | 0.009 | 500 |
| 3.11 | 0.13 | 0.204 | 0.001 | 0.009 | 500 |
| 3.11 | 0.13 | 0.202 | 0.001 | 0.009 | 500 |
| 3.78 | 0.08 | 0.216 | 0.001 | 0.005 | 500 |
| 3.78 | 0.08 | 0.212 | 0.001 | 0.005 | 500 |
| 3.92 | 0.09 | 0.214 | 0.001 | 0.005 | 500 |
| 3.92 | 0.09 | 0.215 | 0.001 | 0.005 | 500 |
| 0.96 | 0.14 | 0.181 | 0.010 | 0.027 | 800 |
| 0.98 | 0.14 | 0.195 | 0.002 | 0.027 | 800 |
| 0.99 | 0.16 | 0.142 | 0.001 | 0.023 | 800 |
| 1.29 | 0.16 | 0.192 | 0.005 | 0.025 | 800 |
| 1.30 | 0.17 | 0.199 | 0.004 | 0.025 | 800 |
| 1.33 | 0.18 | 0.163 | 0.006 | 0.023 | 800 |
| 1.97 | 0.17 | 0.180 | 0.002 | 0.015 | 800 |
| 2.01 | 0.17 | 0.164 | 0.002 | 0.014 | 800 |
| 2.79 | 0.14 | 0.192 | 0.003 | 0.010 | 800 |
| 2.79 | 0.14 | 0.198 | 0.003 | 0.010 | 800 |
| 3.79 | 0.22 | 0.221 | 0.002 | 0.013 | 800 |
| 3.79 | 0.22 | 0.220 | 0.002 | 0.013 | 800 |
| 4.89 | 0.21 | 0.211 | 0.003 | 0.009 | 800 |
| 4.89 | 0.21 | 0.215 | 0.002 | 0.009 | 800 |
| 5.01 | 0.21 | 0.223 | 0.002 | 0.010 | 800 |
| 5.01 | 0.21 | 0.222 | 0.002 | 0.010 | 800 |
| 6.10 | 0.14 | 0.224 | 0.002 | 0.005 | 800 |
| 6.10 | 0.14 | 0.227 | 0.002 | 0.006 | 800 |
| 6.32 | 0.14 | 0.236 | 0.002 | 0.006 | 800 |
| 6.32 | 0.14 | 0.233 | 0.002 | 0.006 | 800 |

Bonhomme et al. also report a combined Lindhard fit

$$
	\begin{align}
	k_{\mathrm{Ge}} = 0.162\pm0.004
	\end{align}
$$

for the range $0.6-6.3\,\mathrm{keV_{nr}}$.

### Other Data Sources to Add Later

| Material | Data availability | Where to extract |
| --- | --- | --- |
| Ge | plotted historical data and data-set ranges | Sarkis et al. 2020, Fig. 5 and Table II |
| Si | plotted data and fit parameters | Sarkis et al. 2020, Fig. 4 and Tables II-IV |
| NaI(Tl) | experimental Na/I QF points | Tretyak 2010 Fig. 13; Cintas et al. 2024 Fig. 11 and supplementary material if available |
| CsI(Tl) | experimental Cs/I QF points | Tretyak 2010 Fig. 11 |
| CaWO$_4$ | O/Ca/W QF data | Tretyak 2010 Fig. 9; CRESST publications cited there |
| LXe/LAr | prefer light/charge yield tables | NEST tables/code or Sorensen-Dahl model below |

When the paper only provides a figure, use a plot digitizer and store the
digitized points in a data file with columns:

```text
material, recoil_species, energy_keVnr, q, q_err_low, q_err_high, source, figure
```

References:

- A. Bonhomme et al., Eur. Phys. J. C 82, 815 (2022),
  [arXiv:2202.03754](https://arxiv.org/abs/2202.03754), Tables 4-6, Fig. 14.
- Y. Sarkis et al., Phys. Rev. D 101, 102001 (2020),
  [arXiv:2001.06503](https://arxiv.org/abs/2001.06503), Table II, Figs. 4-5.
- D. Cintas et al., "A measurement of the sodium and iodine scintillation
  quenching factors across multiple NaI(Tl) detectors to identify systematics",
  [arXiv:2402.12480](https://arxiv.org/abs/2402.12480).


## Birks Law and the Tretyak Model for Scintillators

For scintillators, Tretyak's implementation of Birks law is the most direct
code-ready model:

$$
	\begin{align}
	\frac{dL}{dr} = S\frac{dE/dr}{1+k_B dE/dr}.
	\end{align}
$$

The integrated light yield is

$$
	\begin{align}
	L(E)=\int_0^E \frac{S\,dE'}{1+k_B(dE/dr)(E')}.
	\end{align}
$$

The ion quenching factor is

$$
	\begin{align}
	Q_i(E)=
	\frac{
		\int_0^E \frac{dE'}{1+k_B(dE/dr)_i(E')}
	}{
		\int_0^E \frac{dE'}{1+k_B(dE/dr)_e(E')}
	}.
	\end{align}
$$

The relative light yield sometimes reported in scintillator papers is

$$
	\begin{align}
	R_i(E) =
	\frac{L_i(E)/E}{L_e(E_0)/E_0}
	=
	Q_i(E)\frac{L_e(E)/E}{L_e(E_0)/E_0}.
	\end{align}
$$

If electron light yield is proportional in the relevant range, then
$R_i(E)\simeq Q_i(E)$.

Implementation inputs:

- `energy_keV`: recoil energy grid.
- `stopping_power_ion`: total stopping power for the recoil ion in
  $\mathrm{MeV\,cm^2/g}$ from SRIM/TRIM.
- `stopping_power_electron`: electron stopping power in the same material and
  units from ESTAR.
- `kB`: Birks factor in $\mathrm{g/(MeV\,cm^2)}$.

Implementation details:

- Interpolate stopping powers on energy, preferably in log-log space.
- Numerically integrate over $E'$ using the same energy unit as the stopping
  tables.
- Tretyak usually uses total ion stopping power.  For CaWO$_4$, he also shows an
  electronic-stopping-only fit; expose `stopping_component="total"|"electronic"`
  if both tables are available.

Useful $k_B$ values from Tretyak:

| Material | Recoil / data used | $k_B$ [$\mathrm{g/(MeV\,cm^2)}$] | Source location |
| --- | --- | ---: | --- |
| CsI(Tl) | Cs/I recoils | $3.2\times10^{-3}$ | Fig. 11 |
| CsI(Na) | Cs/I recoils | $5.5\times10^{-3}$ | Fig. 12 |
| NaI(Tl) | Na/I recoils, older data | $3.8\times10^{-3}$ | Fig. 13a |
| NaI(Tl) | low-energy Na data | $6.5\times10^{-3}$ | Fig. 13b |
| NaI(Tl) | alpha data | $1.25\times10^{-3}$ | Fig. 14 |
| CaWO$_4$ | O/Ca/W, total stopping | $8.0\times10^{-3}$ | Fig. 9a |
| CaWO$_4$ | O/Ca/W, electronic stopping | $9.8\times10^{-3}$ | Fig. 9b |
| ZnWO$_4$ | alpha fit, recoils predicted | $9.0\times10^{-3}$ | Fig. 7 |
| LXe | Xe recoils, data set [60] | $3.5\times10^{-4}$ | Fig. 16a |
| LXe | Xe recoils, data set [61] | $1.7\times10^{-3}$ | Fig. 16b |
| LXe | Xe recoils, data set [62] | $5.5\times10^{-4}$ | Fig. 16c |
| LAr | Ar recoils | $1.25\times10^{-3}$ | Fig. 17 |

The spread in $k_B$ is real and reflects different scintillator samples,
temperatures, dopant concentrations, and signal integration windows.  For
NaI(Tl) and CsI(Tl), prefer detector-specific measured QF tables over a universal
Birks parameter when such data exist.

References:

- J. B. Birks, "Scintillations from Organic Crystals: Specific Fluorescence and
  Relative Response to Different Radiations", Proc. Phys. Soc. A 64, 874 (1951).
- J. B. Birks, The Theory and Practice of Scintillation Counting, Pergamon Press
  (1964).
- V. I. Tretyak, "Semi-empirical calculation of quenching factors for ions in
  scintillators", Astropart. Phys. 33, 40 (2010),
  [arXiv:0911.3041](https://arxiv.org/abs/0911.3041), Eqs. (2), (4)-(8), Figs.
  9, 11-17.


## Stopping-Power Fraction Estimates

For a rough non-scintillator estimate one can use the electronic fraction of
the stopping power:

$$
	\begin{align}
	Q_{\mathrm{SP}}(E_R)
	=
	\frac{1}{E_R}
	\int_0^{E_R}
	\frac{S_e(E')}{S_e(E')+S_n(E')}dE',
	\end{align}
$$

where $S_e$ is electronic stopping and $S_n$ is nuclear stopping.  This is not a
complete detector-response model; it ignores charge collection, recombination,
scintillation non-proportionality, lattice effects, and thresholds.

Implementation input table:

```text
energy_keV, electronic_stopping_MeV_cm2_per_g, nuclear_stopping_MeV_cm2_per_g
```

Use this model only for exploratory estimates or as input to the Birks model.

References:

- J. F. Ziegler, J. P. Biersack, and U. Littmark, The Stopping and Range of Ions
  in Matter, Pergamon Press (1985).
- J. F. Ziegler, M. D. Ziegler, and J. P. Biersack, "SRIM - The stopping and
  range of ions in matter (2010)", Nucl. Instrum. Meth. B 268, 1818 (2010).


## Noble-Liquid Light and Charge Yield Model

For LXe and LAr, implement light and charge yields instead of a scalar QF when
possible.  A compact model following Sorensen and Dahl is:

$$
	\begin{align}
	N_q &= n_e+n_\gamma = \frac{E_R f_n(E_R)}{W}, \\
	f_n(E_R) &= Q_{\mathrm{L}}(E_R;k),
	\end{align}
$$

where $W\simeq 13.8\,\mathrm{eV}$ for LXe in Sorensen and Dahl, and for Xe
$k=0.166$ from the Lindhard semi-empirical formula.  If
$\alpha=N_{\mathrm{ex}}/N_i$, then

$$
	\begin{align}
	N_i &= \frac{N_q}{1+\alpha}, \\
	N_{\mathrm{ex}} &= \alpha N_i.
	\end{align}
$$

For nuclear recoils in LXe, Sorensen and Dahl find
$N_{\mathrm{ex}}/N_i\simeq1$ with $\sim15\%$ uncertainty.  Electron escape is
modeled with the Thomas-Imel box model:

$$
	\begin{align}
	\frac{n_e}{N_i} &= \frac{1}{\xi}\ln(1+\xi),\\
	\xi &= N_i\frac{\beta}{4},
	\end{align}
$$

where $\beta\equiv\alpha_{\mathrm{TI}}/(a^2v)$ is the fitted Thomas-Imel
parameter.  The photon count is then

$$
	\begin{align}
	n_\gamma = N_q-n_e.
	\end{align}
$$

Useful LXe parameter choices from Sorensen and Dahl Table I:

| Fit choice | $N_{\mathrm{ex}}/N_i$ | $\beta=\alpha_{\mathrm{TI}}/(a^2v)$ |
| --- | ---: | ---: |
| Manzur data | 0.86 | 0.028 |
| Manzur data, $E_R>7\,\mathrm{keV}$ | 1.05 | 0.025 |
| Manzur + corrected low-E data | 1.04 | 0.030 |
| Aprile data | 1.13 | 0.042 |

Suggested default for a simple LXe yield model:

$$
	\begin{align}
	W &= 13.8\,\mathrm{eV},\\
	k &= 0.166,\\
	N_{\mathrm{ex}}/N_i &= 1.0,\\
	\beta &= 0.030.
	\end{align}
$$

The resulting yields are

$$
	\begin{align}
	L_y(E_R) &= \frac{n_\gamma}{E_R},\\
	Q_y(E_R) &= \frac{n_e}{E_R}.
	\end{align}
$$

with $E_R$ in keV if yields are desired in photons/keV and electrons/keV.

For production LXe/LAr detector simulations, prefer NEST or detector-specific
yield tables because electric-field dependence, recombination fluctuations, and
detector-specific calibration enter the observed S1/S2 response.

References:

- P. Sorensen and C. E. Dahl, "Nuclear recoil energy scale in liquid xenon with
  application to the direct detection of dark matter", Phys. Rev. D 83, 063501
  (2011), [arXiv:1101.6080](https://arxiv.org/abs/1101.6080), Eqs. (2), (3),
  (5)-(7), Table I.
- J. Thomas and D. A. Imel, "Recombination of electron-ion pairs in liquid argon
  and liquid xenon", Phys. Rev. A 36, 614 (1987).
- M. Szydagis et al., "NEST: A Comprehensive Model for Scintillation Yield in
  Liquid Xenon", JINST 6, P10002 (2011),
  [arXiv:1106.1613](https://arxiv.org/abs/1106.1613).
- M. Szydagis et al., "Enhancement of NEST capabilities for simulating low-energy
  recoils in liquid xenon", JINST 8, C10003 (2013),
  [arXiv:1307.6601](https://arxiv.org/abs/1307.6601).


## Implementation Priority in `NuRecoil`

The most useful implementation order is:

1. `TabulatedQuenching`: linear interpolation of measured QF tables, initially
   with the Bonhomme Ge data above.
2. `AdiabaticLindhardQuenching`: $Q_{\mathrm{L}}(E_R;k)[1-\exp(-E_R/\xi)]$.
3. `BindingEnergyLindhardQuenching`: implement the Sarkis ansatz first; later
   replace or supplement it with interpolated numerical-solution tables.
4. `BirksQuenching`: requires stopping-power tables and material-specific $k_B$.
5. `NobleLiquidYield`: return $(L_y,Q_y)$ rather than a scalar QF.
