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
#     name: nurecoil
# ---

# %% [markdown]
# # Nuclear Form Factors
#
# This notebook demonstrates the form-factor models available in
# `nurecoil.cross_section.form_factors` and compares them for typical
# CEvNS target nuclei.
#
# | Model | Class | Reference |
# |-------|-------|-----------|
# | Helm (Lewin-Smith) | `HelmFormFactor` | Lewin & Smith (1996); arXiv:hep-ph/0608035 |
# | Klein-Nystrand | `KleinNystrandFormFactor` | Phys. Rev. C 60 (1999) 014903 |
# | Gaussian | `GaussianFormFactor` | — (fast approximation) |
#
# All form factors satisfy $F(0) = 1$ and are dimensionless.
# Momentum transfer $q = \sqrt{2 M E_R}$ is in **MeV** (natural units).

# %% [markdown]
# ## 1. Imports

# %%
import numpy as np
import matplotlib.pyplot as plt

from nurecoil.nucleus import Nucleus, ISOTOPE_TABLE
from nurecoil.cross_section.form_factors import (
    HelmFormFactor,
    KleinNystrandFormFactor,
    GaussianFormFactor,
)
import nurecoil.constants as C

plt.rcParams.update({"figure.dpi": 120, "font.size": 11})

# %% [markdown]
# ## 2. Target Nuclei

# %%
GE76  = Nucleus(Z=32, A=76)   # Germanium  (CDEX, TEXONO)
CS133 = Nucleus(Z=55, A=133)  # Caesium    (COHERENT CsI)
SI28  = Nucleus(Z=14, A=28)   # Silicon    (light target)

nuclei = {"Ge-76": GE76, "Cs-133": CS133, "Si-28": SI28}

# %% [markdown]
# ## 3. Form Factor vs. Momentum Transfer
#
# We plot $F(q)$ as a function of $q$ [MeV] for each model and each nucleus.

# %%
q = np.linspace(0.0, 300.0, 500)  # MeV

ff_helm  = HelmFormFactor()
ff_kn    = KleinNystrandFormFactor()
ff_gauss = GaussianFormFactor()

fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)

for ax, (name, nuc) in zip(axes, nuclei.items()):
    ax.plot(q, ff_helm(q, nuc),  label="Helm (Lewin-Smith)", lw=2)
    ax.plot(q, ff_kn(q, nuc),   label="Klein-Nystrand",     lw=2, ls="--")
    ax.plot(q, ff_gauss(q, nuc), label="Gaussian",           lw=2, ls=":")
    ax.axhline(0, color="grey", lw=0.5, ls="--")
    ax.set_title(name)
    ax.set_xlabel(r"$q$ [MeV]")
    ax.set_xlim(0, 300)
    ax.set_ylim(-0.3, 1.05)

axes[0].set_ylabel(r"$F(q)$")
axes[1].legend(loc="upper right", fontsize=9)
fig.suptitle("Nuclear Form Factors", y=1.02)
fig.tight_layout()
plt.show()

# %% [markdown]
# ## 4. Form Factor vs. Recoil Energy
#
# In practice $q \approx \sqrt{2 M E_R}$, so it is useful to plot $F$ directly
# against the nuclear recoil energy $E_R$ [keV].

# %%
E_R_keV = np.linspace(0.1, 100.0, 500)   # keV
E_R_MeV = E_R_keV / C.keV_per_MeV        # convert to MeV (internal unit)

fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)

for ax, (name, nuc) in zip(axes, nuclei.items()):
    q_arr = np.sqrt(2.0 * nuc.M * E_R_MeV)  # MeV

    ax.plot(E_R_keV, ff_helm(q_arr, nuc),  label="Helm",          lw=2)
    ax.plot(E_R_keV, ff_kn(q_arr, nuc),   label="Klein-Nystrand", lw=2, ls="--")
    ax.plot(E_R_keV, ff_gauss(q_arr, nuc), label="Gaussian",       lw=2, ls=":")
    ax.axhline(0, color="grey", lw=0.5, ls="--")
    ax.set_title(name)
    ax.set_xlabel(r"$E_R$ [keV]")
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.3, 1.05)

axes[0].set_ylabel(r"$F(q(E_R))$")
axes[1].legend(loc="upper right", fontsize=9)
fig.suptitle(r"Nuclear Form Factors vs. Recoil Energy", y=1.02)
fig.tight_layout()
plt.show()

# %% [markdown]
# ## 5. Helm Parameter Sensitivity
#
# The Lewin-Smith Helm parametrisation has three free parameters:
# skin thickness $s$, surface diffuseness $a$, and radius coefficient $c_\mathrm{coeff}$.
# Here we vary each one around its default to illustrate the sensitivity
# for Ge-76.

# %%
q_sens = np.linspace(0.0, 250.0, 400)

fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)

# --- vary s ---
ax = axes[0]
for s in [0.7, 0.9, 1.1]:
    F = HelmFormFactor(s_fm=s)(q_sens, GE76)
    ax.plot(q_sens, F, label=rf"$s={s}$ fm")
ax.set_title(r"Vary skin thickness $s$")
ax.set_xlabel(r"$q$ [MeV]")
ax.set_ylabel(r"$F(q)$")
ax.legend(fontsize=9)

# --- vary a ---
ax = axes[1]
for a in [0.40, 0.52, 0.65]:
    F = HelmFormFactor(a_fm=a)(q_sens, GE76)
    ax.plot(q_sens, F, label=rf"$a={a}$ fm")
ax.set_title(r"Vary diffuseness $a$")
ax.set_xlabel(r"$q$ [MeV]")
ax.legend(fontsize=9)

# --- vary c_coeff ---
ax = axes[2]
for cc in [1.1, 1.23, 1.35]:
    F = HelmFormFactor(c_coeff=cc)(q_sens, GE76)
    ax.plot(q_sens, F, label=rf"$c_\mathrm{{coeff}}={cc}$ fm")
ax.set_title(r"Vary $c_\mathrm{coeff}$")
ax.set_xlabel(r"$q$ [MeV]")
ax.legend(fontsize=9)

for ax in axes:
    ax.axhline(0, color="grey", lw=0.5, ls="--")
    ax.set_xlim(0, 250)
    ax.set_ylim(-0.3, 1.05)

fig.suptitle("Helm Form Factor — Parameter Sensitivity (Ge-76)", y=1.02)
fig.tight_layout()
plt.show()

# %% [markdown]
# ## 6. Impact on the CEvNS Cross Section
#
# The form factor suppresses the cross section at high recoil energies.
# Here we compare $d\sigma/dE_R$ with the three models for Ge-76
# at a typical reactor neutrino energy $E_\nu = 4$ MeV.

# %%
from nurecoil.cross_section.sm_cevns import SMCEvNS

E_nu = 4.0  # MeV

xs_helm  = SMCEvNS(GE76, form_factor=ff_helm)
xs_kn    = SMCEvNS(GE76, form_factor=ff_kn)
xs_gauss = SMCEvNS(GE76, form_factor=ff_gauss)

E_R_max  = xs_helm.E_R_max(E_nu)
E_R_MeV2 = np.linspace(1e-6, E_R_max * 0.999, 400)
E_R_keV2 = E_R_MeV2 * C.keV_per_MeV

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(E_R_keV2, xs_helm(E_nu, E_R_MeV2),  label="Helm",          lw=2)
ax.plot(E_R_keV2, xs_kn(E_nu, E_R_MeV2),   label="Klein-Nystrand", lw=2, ls="--")
ax.plot(E_R_keV2, xs_gauss(E_nu, E_R_MeV2), label="Gaussian",       lw=2, ls=":")
ax.set_xlabel(r"$E_R$ [keV]")
ax.set_ylabel(r"$d\sigma/dE_R$ [cm$^2$/MeV]")
ax.set_title(rf"CEvNS on Ge-76, $E_\nu = {E_nu}$ MeV")
ax.legend()
ax.set_yscale("log")
fig.tight_layout()
plt.show()
