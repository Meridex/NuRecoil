# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.3
#   kernelspec:
#     display_name: NURECOIL
#     language: python
#     name: python3
# ---

# %% [markdown]
# # SM CEvNS — Standard Model Coherent Elastic Neutrino-Nucleus Scattering
#
# This notebook uses the **NuRecoil** package to compute the SM CEvNS differential
# cross section, total cross section, and detector event rate.
#
# Units follow `misc/plans/units_convention.md`:
# - Energy: **MeV** (internal)
# - Cross section: $\text{cm}^2$
# - Flux: $\#\,/\,\text{MeV}\,/\,\text{cm}^2$

# %% [markdown]
# ## 1. Imports

# %%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.integrate import quad

# NuRecoil 模块
import nurecoil.constants as C
from nurecoil.nucleus import Nucleus, ISOTOPE_TABLE
from nurecoil.cross_section.form_factors import HelmFormFactor
from nurecoil.cross_section.sm_cevns import SMCEvNS
from nurecoil.detector.resolution import CDEXResolution
from nurecoil.flux.phenomenological import PhenomenologicalFlux
from nurecoil.rate import compute_cevns_spectrum, target_count_from_mass

plt.rcParams.update({"figure.dpi": 120, "font.size": 11})

# %% [markdown]
# ## 2. Physical Constants
#
# Values from `nurecoil.constants` (PDG 2024). All internal energies are in **MeV**.
#
# Note: `sin2_theta_W` is the low-$q^2$ running value $\approx 0.23868$,
# appropriate for CEvNS kinematics (not the MS-bar value at $M_Z$).

# %%
print(f"G_F          = {C.G_F:.7e}  MeV⁻²")
print(f"G_F_GeV2     = {C.G_F_GeV2:.7e}  GeV⁻²  (PDG reference)")
print(f"ℏc           = {C.hbar_c:.7f}  MeV·fm")
print(f"sin²θ_W      = {C.sin2_theta_W}  (low-q² running value)")
print(f"u_to_MeV     = {C.u_to_MeV}  MeV")
print(f"cm_per_fm    = {C.cm_per_fm:.1e}")
print(f"keV_per_MeV  = {C.keV_per_MeV:.0f}")

# %% [markdown]
# ## 3. Weak Mixing Angle and Couplings
#
# SM tree-level vector couplings (low-$q^2$):
#
# $$
# g_V^p = \frac{1}{2} - 2\sin^2\!\theta_W, \qquad g_V^n = -\frac{1}{2}
# $$
#
# With radiative corrections (arXiv:2411.03122):
#
# $$
# g_V^n = -0.5117, \qquad g_V^p = \begin{cases}
# 0.0382 & \nu_e \\ 0.0300 & \nu_\mu \\ 0.0256 & \nu_\tau
# \end{cases}
# $$
#
# The weak charge including nuclear form factors is:
#
# $$
# Q_W\!\left(|\vec{q}|\right) = g_V^n\, N\, F\!\left(|\vec{q}|\right)
#                              + g_V^p\, Z\, F\!\left(|\vec{q}|\right)
# $$

# %%
s2w   = C.sin2_theta_W
g_V_p = 0.5 - 2.0 * s2w   # proton weak vector coupling (tree-level)
g_V_n = -0.5               # neutron weak vector coupling (tree-level)

# Radiative-corrected values
g_V_p_rc = {"nu_e": 0.0382, "nu_mu": 0.0300, "nu_tau": 0.0256}
g_V_n_rc = -0.5117

targets = {
    "Ge76":  Nucleus(Z=32, A=76),
    "Si28":  Nucleus(Z=14, A=28),
    "Xe132": Nucleus(Z=54, A=132),
    "Cs133": Nucleus(Z=55, A=133),
    "Ar40":  Nucleus(Z=18, A=40),
}

print(f"{'Nucleus':<8}  {'Z':>4}  {'N':>4}  {'Q_W tree':>10}  {'Q_W ν_e RC':>12}  {'Q_W ν_μ RC':>12}")
print("-" * 58)
for name, nuc in targets.items():
    Q_tree = g_V_n * nuc.N + g_V_p * nuc.Z
    Q_nue  = g_V_n_rc * nuc.N + g_V_p_rc["nu_e"]  * nuc.Z
    Q_numu = g_V_n_rc * nuc.N + g_V_p_rc["nu_mu"] * nuc.Z
    print(f"{name:<8}  {nuc.Z:>4}  {nuc.N:>4}  {Q_tree:>10.3f}  {Q_nue:>12.3f}  {Q_numu:>12.3f}")

# %% [markdown]
# ## 4. CEvNS Differential Cross Section
#
# $$
# \frac{d\sigma}{dE_R} = \frac{G_F^2\, M}{\pi}
# \left(1 - \frac{E_R}{E_\nu}
#      + \frac{1}{2}\!\left(\frac{E_R}{E_\nu}\right)^{\!2}
#      - \frac{M E_R}{2E_\nu^2}\right)
# \left|Q_W\!\left(|\vec{q}|\right)\right|^2
# $$
#
# Momentum transfer: $|\vec{q}| \approx \sqrt{2M E_R}$ [MeV]. The Helm form factor
# modifies $Q_W$ at finite $|\vec{q}|$.
#
# The plot below shows $d\sigma/dE_R$ for $^{76}$Ge at several neutrino energies.

# %%
ge76 = Nucleus(Z=32, A=76)
xs_ge76 = SMCEvNS(ge76, HelmFormFactor())

fig, ax = plt.subplots(figsize=(7, 4.5))

E_nu_values = [2.0, 4.0, 6.0, 8.0]   # MeV
colors = plt.cm.plasma(np.linspace(0.15, 0.85, len(E_nu_values)))

for E_nu, col in zip(E_nu_values, colors):
    E_R_max = xs_ge76.E_R_max(E_nu)
    E_R = np.linspace(0.0, E_R_max * 0.999, 300)   # MeV
    dsig = xs_ge76(E_nu, E_R)                        # cm² / MeV

    # Plot in keV for readability
    ax.plot(E_R * C.keV_per_MeV, dsig, color=col, label=f"$E_\\nu = {E_nu}$ MeV")

ax.set_xlabel("$E_R$ [keV]")
ax.set_ylabel(r"$d\sigma/dE_R\, {\rm [cm^{2} / MeV]}$")
ax.set_title("$^{76}$Ge  SM CEvNS Differential Cross-section")
ax.set_yscale("log")
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 5. Effect of Radiative Corrections and Neutrino Flavor
#
# The `flavor` parameter selects between tree-level couplings (`"tree"`) and
# radiatively corrected couplings (`"nu_e"`, `"nu_mu"`, `"nu_tau"`).
# The difference is a few percent in the total cross section.

# %%
flavors = {
    "tree":   ("Tree-level",  "royalblue",  "-"),
    "nu_e":   (r"$\nu_e$ RC", "tomato",     "--"),
    "nu_mu":  (r"$\nu_\mu$ RC","seagreen",  "-."),
    "nu_tau": (r"$\nu_\tau$ RC","darkorange",":" ),
}

E_nu_plot = 5.0   # MeV
E_R_max_ref = SMCEvNS(ge76).E_R_max(E_nu_plot)
E_R = np.linspace(0.0, E_R_max_ref * 0.999, 400)

fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))

# Left: absolute dσ/dE_R
ax = axes[0]
for flavor, (label, col, ls) in flavors.items():
    xs = SMCEvNS(ge76, HelmFormFactor(), flavor=flavor)
    ax.plot(E_R * C.keV_per_MeV, xs(E_nu_plot, E_R),
            color=col, ls=ls, lw=2, label=label)
ax.set_xlabel("$E_R$ [keV]")
ax.set_ylabel(r"$d\sigma/dE_R\, {\rm [cm^{2} / MeV]}$")
ax.set_title(f"$^{{76}}$Ge, $E_\\nu = {E_nu_plot}$ MeV")
ax.set_yscale("log")
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Right: ratio to tree-level
ax = axes[1]
xs_tree_vals = SMCEvNS(ge76, HelmFormFactor(), flavor="tree")(E_nu_plot, E_R)
for flavor, (label, col, ls) in flavors.items():
    if flavor == "tree":
        continue
    xs = SMCEvNS(ge76, HelmFormFactor(), flavor=flavor)
    ratio = xs(E_nu_plot, E_R) / np.where(xs_tree_vals > 0, xs_tree_vals, np.nan)
    ax.plot(E_R * C.keV_per_MeV, ratio, color=col, ls=ls, lw=2, label=label)
ax.axhline(1.0, color="royalblue", ls="-", lw=1.5, label="Tree-level")
ax.set_xlabel("$E_R$ [keV]")
ax.set_ylabel("Ratio to tree-level")
ax.set_title("Radiative Correction Effect")
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# %% [markdown]
# ## 6. Total Cross Section via $E_R$ Integration
#
# $$
# \sigma(E_\nu) = \int_0^{E_R^{\max}} \frac{d\sigma}{dE_R}\, dE_R,
# \qquad
# E_R^{\max} = \frac{2E_\nu^2}{M + 2E_\nu}
# $$

# %%
def total_xsec(xs_model, E_nu_arr):
    """Numerically integrate dσ/dE_R from 0 to E_R_max for each E_nu."""
    sigma = np.empty_like(E_nu_arr)
    for i, E_nu in enumerate(E_nu_arr):
        E_R_max = xs_model.E_R_max(float(E_nu))
        val, _ = quad(
            lambda E_R: float(xs_model(float(E_nu), E_R)),
            0.0, E_R_max,
            limit=80,
        )
        sigma[i] = val
    return sigma   # cm²

E_nu_grid = np.linspace(0.5, 10.0, 80)   # MeV
sigma_ge76 = total_xsec(xs_ge76, E_nu_grid)

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(E_nu_grid, sigma_ge76 / 1e-40, color="royalblue", lw=2)
ax.set_xlabel("$E_\\nu$ [MeV]")
ax.set_ylabel(r"$\sigma$ [$10^{-40}\,\rm cm^{2}$]")
ax.set_title("$^{76}$Ge  SM CEvNS Total Cross-section")
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

print(fr"σ(E_ν=5 MeV) ≈ {np.interp(5.0, E_nu_grid, sigma_ge76):.3e} cm²")

# %% [markdown]
# ## 7. CEvNS Event Rate for Different Target Nuclei
#
# Convolve the reactor antineutrino flux (Huber–Mueller parametrisation) with the
# SM CEvNS cross section; integrate over the full phase space to get the rate per
# target atom per second:
#
# $$
# R = \int_{E_\nu} \int_{E_R}
#     \frac{d\Phi}{dE_\nu}\,
#     \frac{d\sigma}{dE_R}\;
#     dE_R\, dE_\nu
# $$
#
# Parameters: $P = 3\ \text{GW}$, $L = 30\ \text{m}$, $N_T = 1$ (per target atom).

# %%
from scipy.integrate import dblquad

flux_model = PhenomenologicalFlux()

P_GW  = 3.0   # GW
L_m   = 30.0  # m

E_nu_min, E_nu_max = 1.8, 8.0   # MeV  (reactor ν̄_e window)

def rate_per_atom(nucleus):
    """Integrated rate [events / atom] (no detector smearing)."""
    xs_model = SMCEvNS(nucleus, HelmFormFactor())

    def integrand(E_R, E_nu):
        phi  = float(flux_model(E_nu, P_GW, L_m))
        dsig = float(xs_model(E_nu, E_R))
        return phi * dsig

    val, _ = dblquad(
        integrand,
        E_nu_min, E_nu_max,
        0.0,
        lambda E_nu: xs_model.E_R_max(E_nu),
        epsrel=1e-3,
    )
    return val   # events / atom / s  (flux already per s via GW conversion)

print("Computing rates (this may take ~30 s) …")
rates = {}
for name, nuc in targets.items():
    rates[name] = rate_per_atom(nuc)
    print(f"  {name:<8}  {rates[name]:.3e} events/atom/s")
print("Done.")

# %% [markdown]
# ## 8. Visualisation: Differential and Total Cross Sections
#
# Compare the differential cross-section shape for Ge, Si, and Xe at fixed
# $E_\nu = 5\ \text{MeV}$, and the total cross section $\sigma(E_\nu)$.

# %%
E_nu_plot = 5.0   # MeV

plot_targets = {
    "Ge76":  (Nucleus(Z=32, A=76),  "royalblue"),
    "Si28":  (Nucleus(Z=14, A=28),  "tomato"),
    "Xe132": (Nucleus(Z=54, A=132), "seagreen"),
}

fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

# --- Left: dσ/dE_R ---
ax = axes[0]
for name, (nuc, col) in plot_targets.items():
    xs_m = SMCEvNS(nuc, HelmFormFactor())
    E_R_max = xs_m.E_R_max(E_nu_plot)
    E_R = np.linspace(0.0, E_R_max * 0.999, 400)
    dsig = xs_m(E_nu_plot, E_R)
    ax.plot(E_R * C.keV_per_MeV, dsig, color=col, lw=2, label=f"${name}$")

ax.set_xlabel("$E_R$ [keV]")
ax.set_ylabel(r"$d\sigma/dE_R\, {\rm [cm^{2} / MeV]}$")
ax.set_title(f"Differential Cross-section ($E_\\nu = {E_nu_plot}$ MeV)")
ax.set_yscale("log")
ax.legend()
ax.grid(True, alpha=0.3)

# --- Right: σ(E_nu) ---
ax = axes[1]
for name, (nuc, col) in plot_targets.items():
    xs_m = SMCEvNS(nuc, HelmFormFactor())
    sig  = total_xsec(xs_m, E_nu_grid)
    ax.plot(E_nu_grid, sig / 1e-40, color=col, lw=2, label=f"${name}$")

ax.set_xlabel("$E_\\nu$ [MeV]")
ax.set_ylabel(r"$\sigma$ [$10^{-40}\,\rm cm^{2}$]")
ax.set_title("Total Cross-section vs Neutrino Energy")
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# %% [markdown]
# ## 9. Event-Rate Comparison Across Target Nuclei
#
# The CEvNS cross section scales approximately as $N^2$ (coherent neutron contribution).
# The plots below visualise the rate for each nucleus and verify the $\sigma \propto Q_W^2$
# scaling directly.

# %%
names  = list(rates.keys())
R_vals = np.array([rates[n] for n in names])
N_vals = np.array([targets[n].N for n in names])
Q_vals = np.array([g_V_n * targets[n].N + g_V_p * targets[n].Z for n in names])

fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

bar_colors = ["royalblue", "tomato", "seagreen", "darkorange", "purple"]

# --- Left: event rate bar chart ---
ax = axes[0]
bars = ax.bar(names, R_vals, color=bar_colors)
ax.set_ylabel("Event Rate [events / atom / s]")
ax.set_title(f"CEvNS Rate  ($P={P_GW}$ GW, $L={L_m}$ m)")
ax.set_yscale("log")
ax.bar_label(bars, fmt="%.1e", padding=3, fontsize=8)
ax.grid(True, axis="y", alpha=0.3)

# --- Right: rate vs Q_W^2 (N^2 scaling check) ---
ax = axes[1]
Q2 = Q_vals ** 2
ax.scatter(Q2, R_vals, c=bar_colors, s=100, zorder=5)
for n, x, y in zip(names, Q2, R_vals):
    ax.annotate(n, (x, y), textcoords="offset points", xytext=(5, 3), fontsize=9)

# Fit a proportionality line
slope = np.polyfit(Q2, R_vals, 1)[0]
x_fit = np.linspace(0, Q2.max() * 1.05, 100)
ax.plot(x_fit, slope * x_fit, "k--", alpha=0.4, label=f"$\\propto Q_W^2$")
ax.set_xlabel("$Q_W^2$ ($q\\to 0$ limit)")
ax.set_ylabel("Event Rate [events / atom / s]")
ax.set_title("$\\sigma \\propto Q_W^2$ Check")
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# %% [markdown]
# ## 10. Binned CEvNS Detected-Energy Spectrum
#
# The high-level rate API assembles the reactor flux, SM CEvNS cross section,
# Lindhard quenching, detector resolution, and target normalization into the
# detected-energy bin spectrum
#
# $$
# N_T \int_{\Delta E_{\rm det}} dE_{\rm det}
#     \int dE_\nu \int dE_R\,
#     \frac{d\Phi}{dE_\nu}
#     \frac{d\sigma}{dE_R}
#     f_{\rm res}\!\left(f_Q(E_R)E_R, E_{\rm det}\right).
# $$
#
# Here `N_T` is the number of target atoms.  For a single-isotope target the
# package uses the approximation $N_T = m_{\rm target}/A \times N_A$.

# %%
target_mass_g = 1000.0       # 1 kg Ge76
exposure_days = 1.0
exposure_s = exposure_days * 24.0 * 3600.0

E_det_edges_keV = np.linspace(0.0, 0.30, 9)      # detected energy bins [keVee]
E_det_edges_MeV = E_det_edges_keV / C.keV_per_MeV
E_det_centers_keV = 0.5 * (E_det_edges_keV[:-1] + E_det_edges_keV[1:])
bin_widths_keV = np.diff(E_det_edges_keV)

N_T_ge76 = target_count_from_mass(target_mass_g, ge76)
cevns_bin_rates = compute_cevns_spectrum(
    E_det_edges_MeV,
    flux_model,
    nucleus=ge76,
    target_mass_g=target_mass_g,
    P=P_GW,
    L=L_m,
    E_nu_range=(E_nu_min, E_nu_max),
    resolution=CDEXResolution(),
)
cevns_bin_counts = compute_cevns_spectrum(
    E_det_edges_MeV,
    flux_model,
    nucleus=ge76,
    target_mass_g=target_mass_g,
    P=P_GW,
    L=L_m,
    E_nu_range=(E_nu_min, E_nu_max),
    resolution=CDEXResolution(),
    exposure_s=exposure_s,
)

print(f"Target mass: {target_mass_g / 1000:.1f} kg Ge76")
print(f"N_T: {N_T_ge76:.3e} target atoms")
print(f"Exposure: {exposure_days:.1f} day")
print(f"Total CEvNS rate in plotted window: {cevns_bin_rates.sum():.3e} events/s")
print(f"Total CEvNS counts in plotted window: {cevns_bin_counts.sum():.3e} events")

fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

ax = axes[0]
ax.step(E_det_edges_keV[:-1], cevns_bin_rates, where="post", color="royalblue", lw=2)
ax.set_xlabel("$E_{\\rm det}$ [keVee]")
ax.set_ylabel("Rate [events / s / bin]")
ax.set_title("Binned CEvNS Detected Spectrum")
ax.set_yscale("log")
ax.grid(True, alpha=0.3)

ax = axes[1]
ax.bar(E_det_centers_keV, cevns_bin_counts, width=bin_widths_keV,
       color="royalblue", alpha=0.8, edgecolor="black", linewidth=0.4)
ax.set_xlabel("$E_{\\rm det}$ [keVee]")
ax.set_ylabel(f"Counts / bin / {exposure_days:.0f} day")
ax.set_title("Exposure-Scaled CEvNS Counts")
ax.set_yscale("log")
ax.grid(True, axis="y", alpha=0.3)

plt.tight_layout()
plt.show()
