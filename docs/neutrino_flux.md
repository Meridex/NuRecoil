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


The neutrino spectrum per fission for different actinides

