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
# # Phenomenological vs. Interpolated Flux Models
#
# This notebook compares the two families of reactor antineutrino flux models
# available in `nurecoil`:
#
# | Family | Classes | Description |
# |--------|---------|-------------|
# | **Phenomenological** | `PhenomenologicalFlux` | Exponential polynomial fit to nuclear data |
# | **Interpolated** | `IsotopeInterpolatedFlux`, `InterpolatedFlux` | Tabulated spectra from summation calculations / measurements |
#
# Both share the same fission-fraction weighting and $\bar{e}$ bookkeeping via
# `ReactorMix`, so they can be compared directly.
#
# Topics covered:
# 1. Available models and their energy ranges
# 2. Absolute flux comparison on a common grid
# 3. Ratio plots — all models relative to Estienne 2019
# 4. Effect of fission fraction preset on the ratio
# 5. Per-isotope spectral shape comparison (pheno coefficients vs. isospec tables)

# %%
import numpy as np
import matplotlib.pyplot as plt

from nurecoil.flux import (
    PhenomenologicalFlux,
    IsotopeInterpolatedFlux,
    InterpolatedFlux,
    make_flux,
    list_sources,
    list_isotopes,
)
from nurecoil.flux.phenomenological import SPECTRUM_MODELS, _SPECTRUM_ENERGY_RANGE

# %% [markdown]
# ## 1. Available models and energy ranges

# %%
print("Phenomenological models:")
for model, (lo, hi) in _SPECTRUM_ENERGY_RANGE.items():
    note = "  ← no U238" if model == "huber_2011" else ""
    print(f"  {model:<20}  {lo:.1f} – {hi:.1f} MeV{note}")

print()
print("Interpolated / tabulated sources:")
for src in list_sources():
    isos = list_isotopes(src)
    kind = "isospec  " if isos else "composite"
    flux = make_flux(src)
    print(f"  {src:<22}  {kind}   {flux.E_min:.3f} – {flux.E_max:.2f} MeV")

# %% [markdown]
# ## 2. Absolute flux comparison
#
# Common evaluation grid: 2–8 MeV, where all phenomenological models and most
# isospec sources are valid.  All use the default "typical" PWR fission fractions.

# %%
P = 3.0   # GW thermal power
L = 10.0  # m baseline

ISOSPEC_SOURCES = [s for s in list_sources() if list_isotopes(s)]
COMPOSITE_SOURCES = [s for s in list_sources() if not list_isotopes(s)]

# Style maps
src_colors = {"CEA2023": "C0", "estienne2019": "C1", "mueller2011": "C2", "vogel1989": "C3"}
src_styles  = {
    "CEA2023":      dict(lw=1.5, ls="-"),
    "estienne2019": dict(lw=1.5, ls="--"),
    "mueller2011":  dict(lw=1.5, ls="-."),
    "vogel1989":    dict(lw=1.5, ls=":"),
}
pheno_colors = {"huber_mueller": "k", "mueller_2011": "dimgrey", "huber_2011": "C6", "vogel_1985": "C7"}
pheno_styles  = {
    "huber_mueller": dict(lw=2.0, ls="-"),
    "mueller_2011":  dict(lw=1.5, ls="--"),
    "huber_2011":    dict(lw=1.5, ls="-."),
    "vogel_1985":    dict(lw=1.5, ls=":"),
}

E_cmp = np.linspace(2.0, 8.0, 400)

# Build all flux objects
isospec_fluxes = {src: make_flux(src) for src in ISOSPEC_SOURCES}

pheno_fluxes = {}
for model in SPECTRUM_MODELS:
    kw = {}
    if model == "huber_2011":
        # huber_2011 has no U238 coefficients
        kw["fission_fractions"] = {"U235": 0.600, "U238": 0.000, "Pu239": 0.270, "Pu241": 0.130}
    pheno_fluxes[model] = PhenomenologicalFlux(spectrum_model=model, **kw)

# %%
fig, ax = plt.subplots(figsize=(9, 5))

for src, fx in isospec_fluxes.items():
    phi = fx(E_cmp, P, L)
    ax.semilogy(E_cmp, phi, label=f"interp: {src}",
                color=src_colors[src], **src_styles[src])

for model, fx in pheno_fluxes.items():
    phi = fx(E_cmp, P, L)
    ax.semilogy(E_cmp, phi, label=f"pheno: {model}",
                color=pheno_colors[model], **pheno_styles[model])

ax.set_xlabel(r"$E_\nu$ (MeV)")
ax.set_ylabel(r"$d\Phi/dE_\nu$  (cm$^{-2}$ MeV$^{-1}$)")
ax.set_title(f"All models — absolute flux  (P = {P} GW, L = {L} m, preset = 'typical')")
ax.legend(fontsize=8, ncol=2)
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 3. Ratio to Estienne 2019
#
# Estienne 2019 is used as reference because it covers the full 0.05–10 MeV
# range with 101 well-spaced points from a recent summation calculation.

# %%
flux_ref = isospec_fluxes["estienne2019"]
phi_ref  = flux_ref(E_cmp, P, L)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# --- left: isospec sources vs. Estienne 2019 ---
ax = axes[0]
ax.axhline(1.0, color="C1", lw=2, ls="--", label="Estienne 2019 (ref)")
for src, fx in isospec_fluxes.items():
    if src == "estienne2019":
        continue
    phi   = fx(E_cmp, P, L)
    ratio = np.where(phi_ref > 0, phi / phi_ref, np.nan)
    ax.plot(E_cmp, ratio, label=src, color=src_colors[src], **src_styles[src])
ax.set_xlabel(r"$E_\nu$ (MeV)")
ax.set_ylabel("Ratio to Estienne 2019")
ax.set_title("Interpolated sources — ratio to Estienne 2019")
ax.set_ylim(0.5, 1.5)
ax.axhline(0.0, color="none")  # keep ylim symmetric
ax.legend(fontsize=9)

# --- right: pheno models vs. Estienne 2019 ---
ax = axes[1]
ax.axhline(1.0, color="C1", lw=2, ls="--", label="Estienne 2019 (ref)")
for model, fx in pheno_fluxes.items():
    phi   = fx(E_cmp, P, L)
    ratio = np.where(phi_ref > 0, phi / phi_ref, np.nan)
    ax.plot(E_cmp, ratio, label=model,
            color=pheno_colors[model], **pheno_styles[model])
ax.set_xlabel(r"$E_\nu$ (MeV)")
ax.set_ylabel("Ratio to Estienne 2019")
ax.set_title("Phenomenological models — ratio to Estienne 2019")
ax.set_ylim(0.5, 1.5)
ax.legend(fontsize=9)

plt.tight_layout()
plt.show()

# %% [markdown]
# ## 4. Effect of fission fraction preset on the ratio
#
# Changing the fission fractions shifts $\bar{e}$ and changes the relative
# weights of each isotope's spectrum.  Here we check whether the
# pheno/interpolated ratio is stable across reactor types.

# %%
presets = ["typical", "ksnps", "conus", "daya_bay"]
preset_colors = {"typical": "C0", "ksnps": "C1", "conus": "C2", "daya_bay": "C3"}
preset_ls     = {"typical": "-",  "ksnps": "--", "conus": "-.", "daya_bay": ":"}

# Use huber_mueller as the representative pheno model
fig, ax = plt.subplots(figsize=(8, 5))
for preset in presets:
    pheno  = PhenomenologicalFlux(spectrum_model="huber_mueller",
                                   fission_fractions=preset)
    interp = make_flux("estienne2019", fission_fractions=preset)
    phi_p  = pheno(E_cmp, P, L)
    phi_i  = interp(E_cmp, P, L)
    ratio  = np.where(phi_i > 0, phi_p / phi_i, np.nan)
    ax.plot(E_cmp, ratio, label=preset,
            color=preset_colors[preset], ls=preset_ls[preset], lw=1.5)

ax.axhline(1.0, color="grey", lw=0.8, ls=":")
ax.set_xlabel(r"$E_\nu$ (MeV)")
ax.set_ylabel("Huber–Mueller / Estienne 2019")
ax.set_title("Huber–Mueller / Estienne 2019 ratio — sensitivity to reactor preset")
ax.set_ylim(0.7, 1.3)
ax.legend(fontsize=9)
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 5. Per-isotope spectral shape comparison
#
# The phenomenological models parameterise each actinide's spectrum as
# $S_i(E) = \exp\!\bigl(\sum_k a_k^{(i)} E^{k-1}\bigr)$.
# Here we overlay the pheno coefficients-based curves with the tabulated
# `isotope_spectrum` data from each isospec source for direct visual comparison.

# %%
iso_colors = {"U235": "C0", "U238": "C1", "Pu239": "C2", "Pu241": "C3"}
ISOTOPES   = ["U235", "U238", "Pu239", "Pu241"]

fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=False)
axes = axes.flatten()

for ax, iso in zip(axes, ISOTOPES):
    # isospec sources: tabulated S_i(E)
    for src, fx in isospec_fluxes.items():
        S = fx.isotope_spectrum(iso, E_cmp)
        mask = S > 0
        if not np.any(mask):
            continue
        ax.semilogy(E_cmp[mask], S[mask],
                    label=f"interp: {src}", color=src_colors[src], **src_styles[src])

    # pheno models: evaluate S_i(E) = exp(poly) directly
    for model, coeffs in SPECTRUM_MODELS.items():
        if iso not in coeffs:
            continue
        from nurecoil.flux.phenomenological import _eval_spectrum
        S_pheno = _eval_spectrum(coeffs[iso], E_cmp)
        ax.semilogy(E_cmp, S_pheno,
                    label=f"pheno: {model}", color=pheno_colors[model], **pheno_styles[model])

    iso_label = {"U235": r"$^{235}$U", "U238": r"$^{238}$U",
                 "Pu239": r"$^{239}$Pu", "Pu241": r"$^{241}$Pu"}[iso]
    ax.set_xlabel(r"$E_\nu$ (MeV)")
    ax.set_ylabel(r"$S_i(E_\nu)$  (# / MeV / fission)")
    ax.set_title(iso_label)
    ax.set_xlim(2.0, 8.0)
    ax.legend(fontsize=7)

fig.suptitle("Per-isotope spectrum: phenomenological coefficients vs. tabulated data", y=1.01)
plt.tight_layout()
plt.show()
