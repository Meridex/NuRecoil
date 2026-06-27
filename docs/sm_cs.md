# About the differential cross section in the SM

## The differential Cross Section for CEvNS

From [arXiv: 2203.07361] and complement from [JHEP 05 (20022) 037]

$$
	\begin{align}
	\frac{d\sigma}{dT} = \frac{G_{F}^{2}M}{4\pi}\left( 1-\frac{MT}{2E_{\nu}^{2}}-\frac{T}{E_{\nu}} + \frac{1}{2}\left( \frac{T}{E_{\nu}} \right)^{2}\right) Q_{w}^{2}[F_{w}(q^{2})]^{2} + \frac{G_{F}^{2}M}{4\pi}\left( 1+\frac{MT}{2E_{\nu}^{2}}-\frac{T}{E_{\nu}} \right) F_{A}(q^{2})
	\end{align}
$$
with
$$
	\begin{align}
	Q_{w} &\equiv Z(1-4s_{W}^{2}) - N \\
	F_{w}(0) &= 1
	\end{align}
$$
$T=E_{R} = \frac{q^{2}}{2M}=E_{\nu}-E_{\nu}'$ is the nuclear recoil energy, $T\in\left[ 0, \frac{2E_{\nu}^{2}}{M+2E_{\nu}} \right]$.

In most cases, the axial-vector contribution is suppressed compared with the vector contribution, then we have
$$
	\begin{align}
	\frac{d\sigma}{dT} = \frac{G_{F}^{2}M}{4\pi}\left( 1-\frac{MT}{2E_{\nu}^{2}}-\frac{T}{E_{\nu}} + \frac{1}{2}\left( \frac{T}{E_{\nu}} \right)^{2} \right) Q_{w}^{2}[F_{w}(q^{2})]^{2}
	\end{align}
$$
Asumming $T\ll E_{\nu}$ for CEvNS (note that this assumption may failed for some other cases), we reach
$$
	\begin{align}
	\frac{d\sigma}{dT} = \frac{G_{F}^{2}M}{4\pi}\left( 1- \frac{MT}{2E_{\nu^{2}}} \right) Q_{w}^{2}[F_{w}(q^{2})]^{2}
	\end{align}
$$
There exists another convention as:
$$
	\begin{align}
	\frac{d\sigma}{dT} &= \frac{G_{F}^{2}M}{\pi}\left( 1- \frac{MT}{2E_{\nu}^{2}} \right) Q_{w}'^{2}(F_{w}(q^{2}))^{2} \\
	\text{with }Q_{w}' &= g_{p}^{V}Z + g_{n}^{V}N, \\
	g_{n}^{V} &= -\frac{1}{2} \\
	g_{p}^{V} &= \frac{1}{2}-2s_{W}^{2}
	\end{align}
$$
For the weak mixing angle
- From [JHEP 05 (2022) 037] and [arXiv: hep-ph/0409169] and [arXiv: 1712.09146], $s_{W}^{2} (q^{2}\to 0) = 0.23868$
- From [arXiv: 2411.03122] and [Phys. Rev. D 110 (2024) 030001] $s_{W}^{2}(q^{2}\to 0)=023873\pm 0.00005$

For the couplings with radiation corrections, we have from [arXiv: 2411.03122] and [Eur. Phys. J. C 83(7):683 (2023)]
$$
	\begin{align}
	g_{V}^{n} &= -\frac{1}{2} + \mathrm{r.c.} = -0.5117 \\
	g_{V}^{p}(\nu_{\ell}) &= \frac{1}{2} - 2s_{W}^{2} + \mathrm{r.c.} = \begin{cases}
	0.0382 && \text{for }\nu_{e} \\
	0.0300 && \text{for }\nu_{\mu} \\
	0.0256 && \text{for }\nu_{\tau}
	\end{cases}
	\end{align}
$$


## The differential cross section for EvES

For the scattering with the electron, there are more contributions. [arXiv: 1907.03379] gives a comprehensive discussion about the neutrino-electron scattering. Here, I follow [arXiv: 2411.03122]

$$
	\begin{align}
	\frac{d\sigma^{\mathrm{EvES,SM}}}{dT} = Z_{\mathrm{eff}}^{\mathcal{A}}(E_{R}) \frac{G_{F}^{2}m_{e}}{2\pi}\left[ (g_{V}^{\nu_{\ell}}+g_{A}^{\nu_{\ell}})^{2} + (g_{V}^{\nu_{\ell}}-g_{A}^{\nu_{\ell}})^{2}\left( 1-\frac{T}{E_{\nu}} \right)^{2} - ((g_{V}^{\nu_{\ell}})^{2}-(g_{A}^{\nu_{\ell}})^{2}) \frac{m_{e}T}{E_{\nu}^{2}} \right]
	\end{align}
$$
The couplings are given as [arXiv: 2207.05036] where the weak mixing angle at zero momentum transfer is given as $s_{W}^{2} = 0.23857$
$$
	\begin{align}
    g_{V}^{\nu_{e}} &= 2s_{W}^{2} + \frac{1}{2} + \mathrm{r.c.} = 0.9521,  &\qquad g_{A}^{\nu_{e}} &= \frac{1}{2} + \mathrm{r.c.} = 0.4938 \\
	g_{V}^{\nu_{\mu}} &= 	2s_{w}^{2} - \frac{1}{2} + \mathrm{r.c.} = -0.0397  &\qquad g_{A}^{\nu_{\mu}} &= -\frac{1}{2} + \mathrm{r.c.} = -0.5062 \\
	g_{V}^{\nu_{\tau}} &= 2s_{W}^{2} - \frac{1}{2} + \mathrm{r.c.} = -0.0353  &\qquad g_{A}^{\nu_{\tau}} &= -\frac{1}{2} + \mathrm{r.c.} = -0.5062
	\end{align}
$$
Note that for anti-neutrinos, the axial couplings change sign!


### About the atomic ionization

$Z_{\mathrm{eff}}^{\mathcal{A}}(T)$ quantifies the effective number of electrons which can be ionized at $T$:

#### For Ge

From [JHEP 09 (2022) 164]

$$
	\begin{align}
	Z_{\mathrm{eff}}^{\mathrm{Ge}}(T) = \begin{cases}
	32, & T>11.103\,\mathrm{keV} \\
	30, & 11.103\,\mathrm{keV}\geq T>1.4146\,\mathrm{keV} \\
	28, & 1.4146\,\mathrm{keV}\geq T>1.2481\,\mathrm{keV} \\
	26, & 1.2481\,\mathrm{keV}\geq T>1.217\,\mathrm{keV} \\
	22, & 1.217\,\mathrm{keV}\geq T>0.1801\,\mathrm{keV} \\
	20, & 0.1801\,\mathrm{keV}\geq T>0.1249\,\mathrm{keV} \\
	18, & 0.1249\,\mathrm{keV}\geq T>0.1208\,\mathrm{keV} \\
	14, & 0.1208\,\mathrm{keV}\geq T>0.0298\,\mathrm{keV} \\
	10, & 0.0298\,\mathrm{keV}\geq T>0.0292\,\mathrm{keV} \\
	4, & 0.0292\,\mathrm{keV}\geq T
	\end{cases}
	\end{align}
$$

This is obtained from the Electron binding energies given at https://xdb.lbl.gov, for Ge it reads (given in eV)

| Element                                         | $K\ 1s$ | $L_{1}\ 2s$ | $L_{2}\ 2p_{1 / 2}$ | $L_{3}\ 2p_{3 / 2}$ | $M_{1}\ 3s$ | $M_{2}\ 3p_{1 / 2}$ | $M_{3}\ 3p_{3 / 2}$ | $M_{4}\ 3d_{3 / 2}$ | $M_{5}\ 3d_{5 / 2}$ | $N_{1}\ 4s$ | $N_{2}\ 4p_{1 / 2}$ |
| ----------------------------------------------- | ------- | ----------- | ------------------- | ------------------- | ----------- | ------------------- | ------------------- | ------------------- | ------------------- | ----------- | ------------------- |
| Ge-32                                           | 11103   | 1414.6      | 1248.1              | 1217.0              | 180.1       | 124.9               | 120.8               | 29.8                | 29.2                | -           | -                   |
| # of electron                                   | 2       | 2           | 2                   | 4                   | 2           | 2                   | 4                   | 4                   | 6                   | 2           | 2                   |
| # of electron accumulated from inner to outside | 2       | 4           | 6                   | 10                  | 12          | 14                  | 18                  | 22                  | 28                  | 30          | 32                  |
| # of electron accumulated from outside to inner | 32      | 30          | 28                  | 26                  | 22          | 20                  | 18                  | 14                  | 10                  | 4           | 2                   |

For the out most electrons (belongs to N orbit), they are more of valence band, binding energy is not a good description.

#### For Cs and I

From [JHEP 09 (2022) 164]
$$
	\begin{align}
	Z_{\mathrm{eff}}^{\mathrm{Cs}}(T) = \begin{cases}
	55, & T > 35.99\,\mathrm{keV} \\
	53, & 35.99\,\mathrm{keV}\geq T>5.71\,\mathrm{keV} \\
	51, & 5.71\,\mathrm{keV}\geq T>5.36\,\mathrm{keV} \\
	49, & 5.36\,\mathrm{keV}\geq T>5.01\,\mathrm{keV} \\
	45, & 5.01\,\mathrm{keV}\geq T>1.21\,\mathrm{keV} \\
	43, & 1.21\,\mathrm{keV}\geq T>1.07\,\mathrm{keV} \\
	41, & 1.07\,\mathrm{keV}\geq T>1\,\mathrm{keV} \\
	37, & 1\,\mathrm{keV}\geq T>0.74\,\mathrm{keV} \\
	33, & 0.74\,\mathrm{keV}\geq T>0.73\,\mathrm{keV} \\
	27, & 0.73\,\mathrm{keV}\geq T>0.23\,\mathrm{keV} \\
	25, & 0.23\,\mathrm{keV}\geq T>0.17\,\mathrm{keV} \\
	23, & 0.17\,\mathrm{keV}\geq T>0.16\,\mathrm{keV} \\
	19, & 0.16\,\mathrm{keV}\geq T
	\end{cases}
	\end{align}
$$

This is obtained from the Electron binding energies given at https://xdb.lbl.gov, for Cs it reads (given in eV)

| Element                                         | $K\ 1s$ | $L_{1}\ 2s$ | $L_{2}\ 2p_{1 / 2}$ | $L_{3}\ 2p_{3 / 2}$ | $M_{1}\ 3s$ | $M_{2}\ 3p_{1 / 2}$ | $M_{3}\ 3p_{3 / 2}$ | $M_{4}\ 3d_{3 / 2}$ | $M_{5}\ 3d_{5 / 2}$ | $N_{1}\ 4s$ | $N_{2}\ 4p_{1 / 2}$ | $N_{3}\ 4p_{3 / 2}$ | $N_{4}\ 4d_{3 / 2}$ | $N_{5}\ 4d_{5 / 2}$ | $N_{6}\ 4f_{5 / 2}$ | $N_{7}\ 4f_{7 / 2}$ | $O_{1}\ 5s$ | $O_{2}\ 5p_{1 / 2}$ | $O_{3}\ 5p_{3 / 2}$ | $P_{1}\ 6s$ |
| ----------------------------------------------- | ------- | ----------- | ------------------- | ------------------- | ----------- | ------------------- | ------------------- | ------------------- | ------------------- | ----------- | ------------------- | ------------------- | ------------------- | ------------------- | ------------------- | ------------------- | ----------- | ------------------- | ------------------- | ----------- |
| Cs 55                                           | 35985   | 5714        | 5359                | 5012                | 1211        | 1071                | 1003                | 740.5               | 726.6               | 232.3       | 172.4               | 161.3               | 79.8                | 77.5                | -                   | -                   | 23.3        | 13.4                | 12.1                | -           |
| # of electron                                   | 2       | 2           | 2                   | 4                   | 2           | 2                   | 4                   | 4                   | 6                   | 2           | 2                   | 4                   | 4                   | 6                   | 0 out of 6          | 0 out of 8          | 2           | 2                   | 4                   | 1 out of 2  |
| # of electron accumulated from inner to outside | 2       | 4           | 6                   | 10                  | 12          | 14                  | 18                  | 22                  | 28                  | 30          | 32                  | 36                  | 40                  | 46                  | -                   | -                   | 48          | 50                  | 54                  | 55          |
| # of electron accumulated from outside to inner | 55      | 53          | 51                  | 49                  | 45          | 43                  | 41                  | 37                  | 33                  | 27          | 25                  | 23                  | 19                  | 15                  |                     |                     | 9           | 7                   | 5                   | 1           |



$$
	\begin{align}
	Z_{\mathrm{eff}}^{\mathrm{I}}(T) = \begin{cases}
	53, & T>33.17\,\mathrm{keV} \\
	51, & 33.17\,\mathrm{keV}\geq T>5.19\,\mathrm{keV} \\
	49, & 5.19\,\mathrm{keV}\geq T>4.86\,\mathrm{keV} \\
	47, & 4.86\,\mathrm{keV}\geq T>4.56\,\mathrm{keV} \\
	43, & 4.56\,\mathrm{keV}\geq T>1.07\,\mathrm{keV} \\
	41, & 1.07\,\mathrm{keV}\geq T>0.93\,\mathrm{keV} \\
	39, & 0.93\,\mathrm{keV}\geq T>0.88\,\mathrm{keV} \\
	35, & 0.88\,\mathrm{keV}\geq T>0.63\,\mathrm{keV} \\
	31, & 0.63\,\mathrm{keV}\geq T>0.62\,\mathrm{keV} \\
	25, & 0.62\,\mathrm{keV}\geq T>0.19\,\mathrm{keV} \\
	23, & 0.19\,\mathrm{keV}\geq T>0.124\,\mathrm{keV} \\
	21, & 0.124\,\mathrm{keV}\geq T>0.123\,\mathrm{keV} \\
	17, & 0.123\,\mathrm{keV}\geq T
	\end{cases}
	\end{align}
$$

This is obtained from the Electron binding energies given at https://xdb.lbl.gov, for I it reads (given in eV)

| Element                                         | $K\ 1s$ | $L_{1}\ 2s$ | $L_{2}\ 2p_{1 / 2}$ | $L_{3}\ 2p_{3 / 2}$ | $M_{1}\ 3s$ | $M_{2}\ 3p_{1 / 2}$ | $M_{3}\ 3p_{3 / 2}$ | $M_{4}\ 3d_{3 / 2}$ | $M_{5}\ 3d_{5 / 2}$ | $N_{1}\ 4s$ | $N_{2}\ 4p_{1 / 2}$ | $N_{3}\ 4p_{3 / 2}$ | $N_{4}\ 4d_{3 / 2}$ | $N_{5}\ 4d_{5 / 2}$ | $N_{6}\ 4f_{5 / 2}$ | $N_{7}\ 4f_{7 / 2}$ | $O_{1}\ 5s$ | $O_{2}\ 5p_{1 / 2}$ | $O_{3}\ 5p_{3 / 2}$ |
| ----------------------------------------------- | ------- | ----------- | ------------------- | ------------------- | ----------- | ------------------- | ------------------- | ------------------- | ------------------- | ----------- | ------------------- | ------------------- | ------------------- | ------------------- | ------------------- | ------------------- | ----------- | ------------------- | ------------------- |
| I 53                                            | 33169   | 5188        | 4852                | 4557                | 1072        | 931                 | 875                 | 630.8               | 619.3               | 186         | 123                 | 123                 | 50.6                | 48.9                | -                   | -                   | -           | -                   | -                   |
| # of electron                                   | 2       | 2           | 2                   | 4                   | 2           | 2                   | 4                   | 4                   | 6                   | 2           | 2                   | 4                   | 4                   | 6                   | 0 out of 6          | 0 out of 8          | 2           | 2                   | 3 out of 4          |
| # of electron accumulated from inner to outside | 2       | 4           | 6                   | 10                  | 12          | 14                  | 18                  | 22                  | 28                  | 30          | 32                  | 36                  | 40                  | 46                  | -                   | -                   | 48          | 50                  | 53                  |
| # of electron accumulated from outside to inner | 53      | 51          | 49                  | 47                  | 43          | 41                  | 39                  | 35                  | 31                  | 25          | 23                  | 21                  | 17                  | 13                  | -                   | -                   | 7           | 5                   | 3                   |
