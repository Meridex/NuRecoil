# About the Nuclear Form Factor

c.f. [arXiv: 2203.07361]

In general, the form factor is separately given for proton and neutron, and the total weak form factor is given as
$$
	\begin{align}
	F_{w}(q^{2}) = \frac{1}{Q_{w}}\left[ ZQ_{w}^{p}F_{p}(q^{2}) + NQ_{w}^{n}F_{n}(q^{2}) \right] 
	\end{align}
$$
in aligned with the differential cross section
$$
	\begin{align}
	\frac{d\sigma}{dT} = \frac{G_{F}^{2}M}{4\pi}\left( 1-\frac{MT}{2E_{\nu}^{2}} \right) Q_{w}^{2}(F_{w}(q^{2}))^{2}
	\end{align}
$$

The form factor for either proton or neutron is related to the proton and neutron densities
$$
	\begin{align}
	F_{n}(q^{2}) &= \frac{4\pi}{N} \int dr r^{2} \frac{\sin(qr)}{qr} \rho_{n}(r) \\
	F_{p}(q^{2}) &= \frac{4\pi}{Z} \int dr r^{2} \frac{\sin(qr)}{qr} \rho_{p}(r)
	\end{align}
$$
where $\rho_{n,p}$ are neutron and proton density distributions normalized to the neutron and proton numbers. In a concrete estimation of the form factor, we need to obtain the information of the neutron/proton distributions. However, there are phenomenological form factors that are based on empirical fits to elastic electron scattering data and similar parameterizations are assumed for the neutron form factor.

See [arXiv: 1902.07398] for a detail studies about uncertainties from form factors for CEvNS.

## Helm Form Factor

[Phys. Rev. 104 (1956) 1466]

The nucleon distribution is given by the convolution of a uniform density with radius $R_{0}$ and a Gaussian profile with width $s$ (the surface thickness), this results in the following form factor:
$$
	\begin{align}
	F_{\mathrm{Helm}}(q^{2}) = \frac{3j_{1}(qR_{0})}{qR_{0}}e^{-q^{2}s^{2}/2}
	\end{align}
$$
where $j_{1}(x)$ is the spherical Bessel function of order one and can also be given as
$$
	\begin{align}
	j_{1}(x) = \frac{\sin x}{x^{2}} - \frac{\cos x}{x}
	\end{align}
$$

[TODO] How to obtain $R_{0}$ and $s$ for target nucleon

- From [JHEP 02 (2020) 123, arXiv: 1911.00762], for CSI, $R_{0}=4.83\,\mathrm{fm}$ and $s=0.9\,\mathrm{fm}$
- From [arXiv: hep-ph/0608035], $R_{1}=\sqrt{ c^{2}+\frac{7}{3}\pi^{2}a^{2}-5s^{2} }$, $c\simeq 1.23 A^{1/3}-0.60\,\mathrm{fm}$, $s\simeq 0.9\,\mathrm{fm}$, $a\simeq 0.52\,\mathrm{fm}$
	- See also, [Nuclear and Particle Physics Proceedings 273-275 (2016) 414-418]
	- See also, [1101.3049], where many other types of form factors are listed
	- Proposed by Lewin and Smith [Astropart. Phys. 6 (1996) 87] and fitting to data from Fricke [Atomic Data and Nuclear Data Tables 60 (1995) 177-285]

## Klein-Nystrand Form Factor

[Phys. Rev. C 60 (1999) 014903, arXiv: hep-ph/9902259]

This approach relies on  a surface-diffuse distribution that results from folding a short-range Yukawa potential with range $a_{k}$ over a hard sphere distribution with radius $R_{A}$, and the form factor becomes:
$$
	\begin{align}
	F_{\mathrm{KN}}(q^{2}) = \frac{3j_{1}(qR_{A})}{qR_{A}}\left[ \frac{1}{1+q^{2}a_{k}^{2}} \right] 
	\end{align}
$$

[TODO] How to obtain $R_{A}$ and $a_{k}$

- [arXiv: 2104.01811] $R_{A}=A^{1/3}r_{0}$,  $a=0.7\,\mathrm{fm}$, $r_{0}=1.3\,\mathrm{fm}$