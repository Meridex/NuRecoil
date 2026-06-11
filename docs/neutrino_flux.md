# About the Neutrino Flux from Reactor

$$
\frac{d\Phi}{dE_\nu} = \frac{1}{4\pi L^{2}} \frac{P_{\mathrm{th}}}{\sum_{i}\left( \frac{f_{i}}{F}e_{i} \right)} \sum_{i} \frac{f_{i}}{F} \frac{dN_{i}}{dE_{\nu}}
$$
- $L$ is the distance to the reactor core [m]
- $P_{\mathrm{th}}$ the thermal power of the reactor [PW -> MeV/s]
- $f_{i}$ is the number of fissions from actinide $i$ [fission]
	- $F=\sum_{i}f_{i}$ [fission]
- $e_{i}$ is the effective thermal energy per fission contributed by each actinide $i$ [MeV/fission]
- $\frac{dN_{i}}{dE_{\nu}}$ the cumulative antineutrino spectrum of $i$ per fission [\#/fission/MeV]

With the above definitions of the unit, we obtain the final unit for $\frac{d\Phi}{dE_{\nu}}$ is

$$
	\rm [m^{-2}][MeV\cdot s^{-1}][MeV\cdot fission^{-1}]^{-1}[\#\cdot fission^{-1}\cdot MeV^{-1}] = \# / m^{2} / MeV / s
$$

The user will provide the distance $L$, and the thermal power $P_{\mathrm{th}}$. In principle, user should also provide the number of fission $f_{i}$ or the fission fraction $\frac{f_{i}}{F}$ according to the actual reactor situation. But we could offer some default choices for the fission fraction for common reactor (Also keep in mind that, the fraction will change with time as the fuel burn out)

The choices for $\frac{f_{i}}{F}$ for $i=\rm ^{235}U, ^{238}U, ^{239}Pu, ^{241}Pu$
- For Kuo-Sheng Nuclear Power Station [arXiv: 2402.06416]
	- U235: 55%
	- U238: 7%
	- Pu239: 32%
	- Pu241: 6%
- For CONUS [arXiv:2401.07684]
	- U235: 49.1%
	- U238: 7.4%
	- Pu239: 36.1%
	- Pu241: 7.4%
- Typical Commercial Reactor [arXiv: 2310.113070]
	- U235: 58%
	- U238: 8%
	- Pu239: 29%
	- Pu241: 5%
- Daya Bay [arXiv: 2210.01068]
	- U235: 56.4%
	- U238: 7.6%
	- Pu239: 30.4%
	- Pu241: 5.6%

Another quantity we need to implement is the effective thermal energy released per fission for different actinides $e_{i}$, (including both thermal fission and neutron capture)

- From [Journal of Nuclear Energy 23 (1969) 5117, Table A.2, last column]
	- U235: $201.7\pm 0.6\,\rm MeV / fission$
	- U238: $205.0\pm 0.9\,\rm MeV / fission$
	- Pu239: $210.0\pm 0.9\,\rm MeV / fission$
	- Pu241: $212.4\pm 1.0\,\rm MeV / fission$
- From [arXiv: hep-ph/0410100, Table 4]
	- U235: $201.92\pm 0.46\,\rm MeV / fission$
	- U238: $205.52\pm 0.96\,\rm MeV / fission$
	- Pu239: $209.99\pm 0.60\,\rm MeV / fission$
	- Pu241: $213.60\pm 0.65\,\rm MeV / fission$
- From [Phys.Rev.C 88 (2013) 0114605, Table IX]
	- U235: $202.36\pm 0.26\,\rm MeV / fission$
	- U238: $205.99\pm 0.52\, \rm MeV / fission$
	- Pu239: $211.12\pm 0.34\, \rm MeV / fission$
	- Pu241: $214.26\pm 0.33\,\rm MeV / fission$
-


The neutrino spectrum per fission for different actinides

- From [arXiv: 1101.2663, TABLE VI] valid from 2 MeV to 8 MeV
	- $\frac{dN_{k}}{dE_{\nu}}\equiv S_{k}(E_{\nu}) = \exp\left( \sum_{p=1}^{6}\alpha_{pk}E_{\nu}^{p-1} \right)$
	- $k=U_{235}$:
		- $\alpha_{1}=3.217$
		- $\alpha_{2}=-3.111$
		- $\alpha_{3}=1.395$
		- $\alpha_{4}=-3.690\times 10^{-1}$
		- $\alpha_{5}=4.445\times 10^{-2}$
		- $\alpha_{6}=-2.053\times 10^{-3}$
	- $k=U_{238}$:
		- $\alpha_{1}=4.833\times 10^{-1}$
		- $\alpha_{2}=1.927\times 10^{-1}$
		- $\alpha_{3}=-1.283\times 10^{-1}$
		- $\alpha_{4}=-6.762\times 10^{-3}$
		- $\alpha_{5}=2.233\times 10^{-3}$
		- $\alpha_{6}=-1.536 \times 10^{-4}$
	- $k=Pu_{239}$:
		- $\alpha_{1}=6.413$
		- $\alpha_{2}=-7.432$
		- $\alpha_{3}=3.535$
		- $\alpha_{4}=-8.820\times 10^{-1}$
		- $\alpha_{5}=1.025\times 10^{-1}$
		- $\alpha_{6}=-4.550\times 10^{-3}$
	- $k=Pu_{241}$:
		- $\alpha_{1}=3.251$
		- $\alpha_{2}=-3.204$
		- $\alpha_{3}=1.428$
		- $\alpha_{4}=-3.675 \times 10^{-1}$
		- $\alpha_{5}=4.254\times 10^{-2}$
		- $\alpha_{6}=-1.896 \times 10^{-3}$
- From [Phys.Rev.D 39 (1989) 3378, Table I] From 2 MeV to 8 MeV
	- $\frac{dN_{\nu}}{dE_{\nu}} = \exp(a_{0}+a_{1}E_{\nu}+a_{2}E_{\nu}^{2})$
	- When applied to energy from 8 MeV to 12 MeV, the spectra are overestimated by a factor of 2 - 3
	- U235
		- $a_{0}=0.870$
		- $a_{1}=-0.160$
		- $a_{2}=-0.0910$
	- U238
		- $a_{0}=0.976$
		- $a_{1}=-0.162$
		- $a_{2}=-0.0790$
	- Pu239
		- $a_{0}=0.896$
		- $a_{1}=-0.239$
		- $a_{2}=-0.0981$
	- Pu241
		- $a_{0}=0.793$
		- $a_{1}=-0.080$
		- $a_{2}=-0.1085$
- From [arXiv: 1106.0687, Table III], just for reference, does not have result for U238, cannot be used
	- $\frac{dN_{\nu}}{dE_{\nu}}=\exp\left( \sum_{i=1}^{6}\alpha_{i}E_{\nu}^{i-1} \right)$
	- U235:
		- $\alpha_{1}=4.367$
		- $\alpha_{2}=-4.577$
		- $\alpha_{3}=2.100$
		- $\alpha_{4}=-5.294\times 10^{-1}$
		- $\alpha_{5}=6.186\times 10^{-2}$
		- $\alpha_{6}=-2.777\times 10^{-3}$
	- Pu239:
		- $\alpha_{1}=4.757$
		- $\alpha_{2}=-5.392$
		- $\alpha_{3}=2.563$
		- $\alpha_{4}=-6.596\times 10^{-1}$
		- $\alpha_{5}=7.820\times 10^{-2}$
		- $\alpha_{6}=-3.536\times 10^{-3}$
	- Pu241:
		- $\alpha_{1}=2.990$
		- $\alpha_{2}=-2.882$
		- $\alpha_{3}=1.278$
		- $\alpha_{4}=-3.343\times 10^{-1}$
		- $\alpha_{5}=3.905\times 10^{-2}$
		- $\alpha_{6}=-1.754\times 10^{-3}$





Instead of using the above fitting formula, we can also use data table from references (here, I just list the table number, for detailed data, please check the references)
- From [arXiv: 1101.2663, TABLE III, IV, V], mueller2011, from 2 MeV to 8 MeV, for U235, U238, Pu239, Pu241 respectively
- From [Phys.Rev.D 39 (1989) 3378, Table II], vogel1989, Below 2 MeV, for U235, U238, Pu239, Pu241 respectively
- From [arXiv: 2304.14992, Supplementary materials], CEA2023, from 0 to 12.5 MeV, for U235, U238, Pu239, Pu241 respectively
- From [Phys.Rev.Lett. 123 (2019) 022502, Supplementary materials], estienne2019, from 0.0 to 10.0 MeV, for U235, U238, Pu239, Pu241 respectively
- From [arXiv: hep-ph/9904384, Fig.1], kopeikin1999 from 0 to 1.5 MeV, total spectrum
- From [Physics of Atomic Nuclei 75 (2012) 2, 143, Table 3], kopeikin2012, from 0.01 MeV to 9 MeV, total spectrum
- From [https://inspirehep.net/literature/2890702, pp 81-82] not used, coved by above list


## Todo

- [ ] Add the neutrino flux calculation using CONFLUX
