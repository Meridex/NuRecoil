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
# # SM EvES — Standard Model Elastic neutrino-Electron Scattering
#
# This notebook uses the **NuRecoil** package to compute the SM EvES differential
# cross section and its dependence on neutrino flavor, neutrino vs anti-neutrino,
# and the atomic ionization threshold effect ($Z_{\rm eff}$).
#
# Units:
# - Energy: **MeV** (internal); plots shown in **keV** where appropriate
# - Cross section: $\text{cm}^2 / \text{MeV}$

# %% [markdown]
# ## 1. Imports

# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

import nurecoil.constants as C
from nurecoil.nucleus import Nucleus
from nurecoil.cross_section.sm_eves import SMEvES

plt.rcParams.update({"figure.dpi": 120, "font.size": 11})

# %% [markdown]
# ## 2. EvES Cross-section Formula
#
# $$
# \frac{d\sigma^{\rm EvES}}{dT} =
#     Z_{\rm eff}^{\mathcal{A}}(T)\,
#     \frac{G_F^2\, m_e}{2\pi}
#     \Bigl[
#         (g_V + g_A)^2
#         + (g_V - g_A)^2 \left(1 - \frac{T}{E_\nu}\right)^{\!2}
#         - (g_V^2 - g_A^2)\,\frac{m_e\, T}{E_\nu^2}
#     \Bigr]
# $$
#
# Kinematic upper bound: $T_{\max} = 2E_\nu^2 / (m_e + 2E_\nu)$.
#
# Flavor-dependent couplings (with radiative corrections, arXiv:2207.05036):
#
# | Flavor | $g_V$ | $g_A$ |
# |--------|--------|--------|
# | $\nu_e$ | $+0.9521$ | $+0.4938$ |
# | $\nu_\mu$ | $-0.0397$ | $-0.5062$ |
# | $\nu_\tau$ | $-0.0353$ | $-0.5062$ |
#
# For anti-neutrinos: $g_A \to -g_A$.

# %%
print(f"m_e       = {C.m_e} MeV")
print(f"G_F       = {C.G_F:.7e} MeV⁻²")

# Show T_max for a few neutrino energies
print("\nKinematic upper bound T_max = 2E_ν² / (m_e + 2E_ν):")
print(f"{'E_ν [MeV]':>12}  {'T_max [MeV]':>14}  {'T_max [keV]':>14}")
print("-" * 45)
for E_nu in [1.0, 2.0, 5.0, 10.0]:
    T_max = 2 * E_nu**2 / (C.m_e + 2 * E_nu)
    print(f"{E_nu:>12.1f}  {T_max:>14.4f}  {T_max * C.keV_per_MeV:>14.2f}")

# %% [markdown]
# ## 3. Differential Cross Section — Flavor Comparison
#
# The three neutrino flavors have very different vector couplings:
# $\nu_e$ couples through both CC and NC ($g_V \approx +0.95$),
# while $\nu_\mu$ and $\nu_\tau$ couple through NC only ($g_V \approx -0.04$).

# %%
ge76 = Nucleus(Z=32, A=76)

flavors = {
    r"$\nu_e$":    ("nu_e",   "royalblue",  "-"),
    r"$\nu_\mu$":  ("nu_mu",  "tomato",     "--"),
    r"$\nu_\tau$": ("nu_tau", "seagreen",   "-."),
}

E_nu = 5.0   # MeV
T_max = SMEvES(ge76, flavor="nu_e").E_R_max(E_nu)
T_arr = np.linspace(0.0, T_max * 0.999, 500)   # MeV

fig, ax = plt.subplots(figsize=(7, 4.5))
for label, (flavor, col, ls) in flavors.items():
    xs = SMEvES(ge76, flavor=flavor)
    ax.plot(T_arr * C.keV_per_MeV, xs(E_nu, T_arr),
            color=col, ls=ls, lw=2, label=label)

ax.set_xlabel("$T$ [keV]")
ax.set_ylabel(r"$d\sigma/dT\, {\rm [cm^{2}/MeV]}$")
ax.set_title(f"SM EvES Differential Cross-section  ($E_\\nu = {E_nu}$ MeV, Ge)")
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 4. Neutrino vs Anti-neutrino
#
# For anti-neutrinos the axial coupling changes sign ($g_A \to -g_A$), which
# swaps the roles of the $(g_V + g_A)^2$ and $(g_V - g_A)^2$ terms.
# The effect is most visible at large $T/E_\nu$ where the $(1 - T/E_\nu)^2$
# suppression differentiates the two terms.

# %%
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))

for i, flavor in enumerate(["nu_e", "nu_mu"]):
    ax = axes[i]
    xs_nu    = SMEvES(ge76, flavor=flavor, antineutrino=False)
    xs_nubar = SMEvES(ge76, flavor=flavor, antineutrino=True)

    T_max_here = xs_nu.E_R_max(E_nu)
    T_plot = np.linspace(0.0, T_max_here * 0.999, 500)

    dsig_nu    = xs_nu(E_nu, T_plot)
    dsig_nubar = xs_nubar(E_nu, T_plot)

    ax.plot(T_plot * C.keV_per_MeV, dsig_nu,    lw=2, label=f"$\\{flavor}$")
    ax.plot(T_plot * C.keV_per_MeV, dsig_nubar, lw=2, ls="--",
            label=f"$\\bar{{\\{flavor}}}$")

    ax.set_xlabel("$T$ [keV]")
    ax.set_ylabel(r"$d\sigma/dT\, {\rm [cm^{2}/MeV]}$")
    ax.set_title(f"{flavor.replace('_', ' ')}  vs anti  ($E_\\nu = {E_nu}$ MeV)")
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# %% [markdown]
# ## 5. Atomic Ionization Threshold — $Z_{\rm eff}(T)$
#
# $Z_{\rm eff}^{\mathcal{A}}(T)$ counts the number of electrons that can be ionized
# at recoil energy $T$.  It is a step function derived from atomic binding energies
# (JHEP 09 (2022) 164; data from https://xdb.lbl.gov).
#
# Currently implemented for **Ge** (Z=32), **Cs** (Z=55), and **I** (Z=53).

# %%
# Visualise Z_eff by plotting dσ/dT with a fine T grid so the steps are visible
targets_zeff = {
    "Ge76":  (Nucleus(Z=32, A=76),  "royalblue"),
    "Cs133": (Nucleus(Z=55, A=133), "tomato"),
    "I127":  (Nucleus(Z=53, A=127), "seagreen"),
}

E_nu_zeff = 50.0   # MeV — high enough that T_max > all binding energy thresholds

fig, ax = plt.subplots(figsize=(9, 5))
for name, (nuc, col) in targets_zeff.items():
    xs = SMEvES(nuc, flavor="nu_e")
    T_max_here = xs.E_R_max(E_nu_zeff)
    T_plot = np.linspace(1e-5, T_max_here * 0.999, 5000)
    dsig = xs(E_nu_zeff, T_plot)
    ax.plot(T_plot * C.keV_per_MeV, dsig, lw=1.5, label=name, color=col)

ax.set_xlabel("$T$ [keV]")
ax.set_ylabel(r"$d\sigma/dT\, {\rm [cm^{2}/MeV]}$")
ax.set_title(f"EvES $Z_{{\\rm eff}}$ Step Structure  ($E_\\nu = {E_nu_zeff}$ MeV, $\\nu_e$)")
ax.set_xscale("log")
ax.legend()
ax.grid(True, alpha=0.3, which="both")
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 6. Total EvES Cross Section vs Neutrino Energy

# %%
def total_eves_xsec(xs_model, E_nu_arr):
    """Numerically integrate dσ/dT from 0 to T_max for each E_nu."""
    sigma = np.empty_like(E_nu_arr)
    for i, E_nu in enumerate(E_nu_arr):
        T_max = xs_model.E_R_max(float(E_nu))
        val, _ = quad(
            lambda T: float(xs_model(float(E_nu), T)),
            0.0, T_max,
            limit=100,
        )
        sigma[i] = val
    return sigma   # cm²

E_nu_grid = np.linspace(0.5, 10.0, 60)   # MeV

fig, ax = plt.subplots(figsize=(7, 4.5))
for label, (flavor, col, ls) in flavors.items():
    xs = SMEvES(ge76, flavor=flavor)
    sig = total_eves_xsec(xs, E_nu_grid)
    ax.plot(E_nu_grid, sig / 1e-44, color=col, ls=ls, lw=2, label=label)

ax.set_xlabel("$E_\\nu$ [MeV]")
ax.set_ylabel(r"$\sigma$ [$10^{-44}\,\rm cm^{2}$]")
ax.set_title("SM EvES Total Cross-section (Ge76)")
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 7. CEvNS vs EvES — Relative Contribution
#
# For completeness, compare the total CEvNS and EvES cross sections at the same
# target (Ge76) to see their relative magnitudes.  CEvNS dominates at MeV energies
# due to coherent $N^2$ enhancement.

# %%
from nurecoil.cross_section.sm_cevns import SMCEvNS
from nurecoil.cross_section.form_factors import HelmFormFactor

xs_cevns = SMCEvNS(ge76, HelmFormFactor(), flavor="nu_e")
xs_eves  = SMEvES(ge76, flavor="nu_e")

def total_cevns_xsec(xs_model, E_nu_arr):
    sigma = np.empty_like(E_nu_arr)
    for i, E_nu in enumerate(E_nu_arr):
        T_max = xs_model.E_R_max(float(E_nu))
        val, _ = quad(
            lambda E_R: float(xs_model(float(E_nu), E_R)),
            0.0, T_max, limit=80,
        )
        sigma[i] = val
    return sigma

sig_cevns = total_cevns_xsec(xs_cevns, E_nu_grid)
sig_eves  = total_eves_xsec(xs_eves,  E_nu_grid)

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(E_nu_grid, sig_cevns / 1e-40, lw=2, color="royalblue", label="CEvNS ($^{76}$Ge nucleus)")
ax.plot(E_nu_grid, sig_eves  / 1e-40, lw=2, color="tomato",    label=r"EvES ($\nu_e$, 32 electrons)")
ax.set_xlabel("$E_\\nu$ [MeV]")
ax.set_ylabel(r"$\sigma$ [$10^{-40}\,\rm cm^{2}$]")
ax.set_title("CEvNS vs EvES — Total Cross-section ($^{76}$Ge, $\\nu_e$)")
ax.set_yscale("log")
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

print("At E_ν = 5 MeV:")
print(f"  CEvNS: {np.interp(5.0, E_nu_grid, sig_cevns):.3e} cm²")
print(f"  EvES:  {np.interp(5.0, E_nu_grid, sig_eves):.3e} cm²")
print(f"  Ratio CEvNS/EvES: {np.interp(5.0, E_nu_grid, sig_cevns) / np.interp(5.0, E_nu_grid, sig_eves):.1f}")
