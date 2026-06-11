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
# # Interpolated Flux — Tabulated Data
#
# This notebook demonstrates how to use the interpolated flux classes with the
# bundled tabulated spectra.
#
# Two internal modes are available:
# - **Isotope mode** (`IsotopeInterpolatedFlux`): per-actinide spectra combined
#   on the fly with fission fractions (e.g. Estienne 2019, Mueller 2011)
# - **Composite mode** (`InterpolatedFlux`): pre-weighted total spectrum
#   (e.g. Kopeikin 2012)
#
# The unified constructor `make_flux` automatically detects which mode to use.
#
# Topics covered:
# 1. Discovering available data sources
# 2. Loading raw spectrum tables
# 3. `make_flux` — unified constructor
# 4. Per-isotope breakdown (`isotope_spectrum` / `isotope_flux`)
# 5. Composite spectrum (`InterpolatedFlux.from_composite`)
# 6. Fission fraction presets comparison

# %%
import numpy as np
import matplotlib.pyplot as plt

from nurecoil.flux import (
    InterpolatedFlux,
    IsotopeInterpolatedFlux,
    make_flux,
    load_spectrum,
    load_all_isotopes,
    list_sources,
    list_isotopes,
)

# %% [markdown]
# ## 1. Discover available data sources

# %%
print("Available sources:", list_sources())
print()
for src in list_sources():
    isos = list_isotopes(src)
    kind = "isospec  " if isos else "composite"
    print(f"  {src:<22}  {kind}  {isos if isos else ''}")

# %% [markdown]
# ## 2. Raw per-isotope spectra — all isospec sources
#
# All four isospec sources (CEA 2023, Estienne 2019, Mueller 2011, Vogel 1989)
# are plotted together.  Each panel shows one actinide; different sources are
# distinguished by line style.  Note the different energy ranges and point
# densities across sources.

# %%
ISOSPEC_SOURCES = [s for s in list_sources() if list_isotopes(s)]
COMPOSITE_SOURCES = [s for s in list_sources() if not list_isotopes(s)]
ISOTOPES = ["U235", "U238", "Pu239", "Pu241"]

# Load all isospec tables: raw_data[source][isotope] = (E, S)
raw_data = {src: {} for src in ISOSPEC_SOURCES}
for src in ISOSPEC_SOURCES:
    for iso in ISOTOPES:
        raw_data[src][iso] = load_spectrum(src, iso)

print("Isospec sources:", ISOSPEC_SOURCES)
print("Composite sources:", COMPOSITE_SOURCES)
print()
for src in ISOSPEC_SOURCES:
    for iso in ISOTOPES:
        E, S = raw_data[src][iso]
        print(f"  {src:<14} {iso}:  {len(E):>3} pts,  {E[0]:.3f} – {E[-1]:.3f} MeV")

# %%
src_styles = {
    "CEA2023":      dict(lw=1.5, ls="-"),
    "estienne2019": dict(lw=1.5, ls="--"),
    "mueller2011":  dict(lw=1.5, ls="-."),
    "vogel1989":    dict(lw=1.5, ls=":"),
}
src_colors = {
    "CEA2023":      "C0",
    "estienne2019": "C1",
    "mueller2011":  "C2",
    "vogel1989":    "C3",
}
iso_titles = {
    "U235":  r"$^{235}$U",
    "U238":  r"$^{238}$U",
    "Pu239": r"$^{239}$Pu",
    "Pu241": r"$^{241}$Pu",
}

fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=False)
axes = axes.flatten()

for ax, iso in zip(axes, ISOTOPES):
    for src in ISOSPEC_SOURCES:
        E, S = raw_data[src][iso]
        ax.semilogy(E, S, label=src,
                    color=src_colors[src], **src_styles[src])
    ax.set_xlabel(r"$E_\nu$ (MeV)")
    ax.set_ylabel(r"$S_i(E_\nu)$  (# / MeV / fission)")
    ax.set_title(iso_titles[iso])
    ax.set_xlim(0, 13)
    ax.legend(fontsize=8)

fig.suptitle("Per-isotope antineutrino spectra — all isospec sources", y=1.01)
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 3. Fission-fraction-weighted total spectra — all isospec sources
#
# `make_flux` combines the per-actinide spectra using the fission fractions
# $f_i$ (default: "typical" PWR).  The total flux is:
#
# $$\frac{d\Phi}{dE_\nu} = \frac{P_{\rm th}}{4\pi L^2 \bar{e}}
# \sum_i f_i \, S_i(E_\nu)$$
#
# Here we compare the weighted total spectrum from each isospec source,
# plus the two composite sources for reference.

# %%
P = 3.0   # GW thermal power
L = 10.0  # m baseline

# Build flux objects for every source (auto-detect mode)
fluxes = {src: make_flux(src) for src in list_sources()}

print("Flux objects created:")
for src, fx in fluxes.items():
    print(f"  {src:<22}  {type(fx).__name__:<30}  "
          f"E = {fx.E_min:.3f} – {fx.E_max:.2f} MeV")

# %%
# Two energy grids: high-energy (isospec range) and full range
E_hi  = np.linspace(1.8, 8.0,  400)   # all isospec sources valid here
E_full = np.linspace(0.01, 10.0, 600)  # wider view

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# --- left: isospec sources on a common high-energy grid ---
ax = axes[0]
for src in ISOSPEC_SOURCES:
    fx = fluxes[src]
    phi = fx(E_hi, P, L)
    ax.semilogy(E_hi, phi, label=src,
                color=src_colors[src], **src_styles[src])
ax.set_xlabel(r"$E_\nu$ (MeV)")
ax.set_ylabel(r"$d\Phi/dE_\nu$  (cm$^{-2}$ MeV$^{-1}$)")
ax.set_title(f"Isospec sources — weighted total flux\n(P = {P} GW, L = {L} m, preset = 'typical')")
ax.legend(fontsize=9)

# --- right: all sources on full range (log-log to show composite sub-MeV) ---
ax = axes[1]
composite_styles = {
    "kopeikin2012": dict(lw=2, ls="-",  color="C4"),
    "kopeikin1999": dict(lw=2, ls="--", color="C5"),
}
for src in ISOSPEC_SOURCES:
    fx = fluxes[src]
    E_src = np.linspace(fx.E_min + 1e-3, fx.E_max - 1e-3, 500)
    phi = fx(E_src, P, L)
    ax.loglog(E_src, phi, label=src,
              color=src_colors[src], **src_styles[src])
for src in COMPOSITE_SOURCES:
    fx = fluxes[src]
    E_src = np.linspace(fx.E_min + 1e-4, fx.E_max - 1e-4, 500)
    phi = fx(E_src, P, L)
    ax.loglog(E_src, phi, label=f"{src} (composite)", **composite_styles[src])
ax.set_xlabel(r"$E_\nu$ (MeV)")
ax.set_ylabel(r"$d\Phi/dE_\nu$  (cm$^{-2}$ MeV$^{-1}$)")
ax.set_title("All sources — full energy range (log-log)")
ax.legend(fontsize=8)

plt.tight_layout()
plt.show()

# %% [markdown]
# ## 4. Per-isotope breakdown
#
# `IsotopeInterpolatedFlux` exposes two per-isotope methods:
# - `isotope_spectrum(iso, E)` — raw S_i(E) [# / MeV / fission], no
#   fission-fraction weighting
# - `isotope_flux(iso, E, P, L)` — partial flux contribution of isotope i,
#   with fission-fraction weight applied
#
# The per-isotope fluxes sum exactly to the total flux.
#
# > **Note on the right panel below:** the fission fractions $f_i$ are
# > energy-independent constants (they depend only on reactor type, not on
# > $E_\nu$).  The ratio $f_i \Phi_i / \Phi_\mathrm{total}$ still varies with
# > energy because each actinide has a **different spectral shape** $S_i(E)$.
# > Harder-spectrum isotopes (U238, Pu241) contribute relatively more at high
# > energy even when their $f_i$ is small.

# %%
iso_colors = {"U235": "C0", "U238": "C1", "Pu239": "C2", "Pu241": "C3"}

# Use Estienne 2019 as the reference isospec source for this breakdown
flux_e = fluxes["estienne2019"]
E_range = np.linspace(1.8, 8.0, 300)

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# --- left: raw S_i(E) from each isospec source, one panel per source ---
# Show Estienne 2019 as reference; all 4 isotopes
ax = axes[0]
for iso in flux_e.isotopes:
    S = flux_e.isotope_spectrum(iso, E_range)
    ax.semilogy(E_range, S, label=iso, color=iso_colors[iso])
ax.set_xlabel(r"$E_\nu$ (MeV)")
ax.set_ylabel(r"$S_i(E_\nu)$  (# / MeV / fission)")
ax.set_title("Per-isotope raw spectra\n(Estienne 2019, no fission-fraction weight)")
ax.legend()
ax.set_xlim(1.8, 8.0)

# --- right: f_i * Phi_i and their sum for each isospec source ---
ax = axes[1]
phi_total = flux_e(E_range, P, L)
for iso in flux_e.isotopes:
    phi_iso = flux_e.isotope_flux(iso, E_range, P, L)
    ax.semilogy(E_range, phi_iso,
                label=rf"$f_i \cdot \Phi_i$  {iso}", color=iso_colors[iso])
ax.semilogy(E_range, phi_total, "k--", lw=2, label="Total")
ax.set_xlabel(r"$E_\nu$ (MeV)")
ax.set_ylabel(r"$d\Phi_i/dE_\nu$  (cm$^{-2}$ MeV$^{-1}$)")
ax.set_title(f"Fission-fraction-weighted contributions\n(Estienne 2019, P = {P} GW, L = {L} m)")
ax.legend(fontsize=8)

plt.tight_layout()
plt.show()

# %%
# Fractional contribution of each isotope to total flux
# f_i is constant; the variation with E comes entirely from different spectral shapes S_i(E)
fig, axes = plt.subplots(1, len(ISOSPEC_SOURCES), figsize=(14, 4), sharey=True)

for ax, src in zip(axes, ISOSPEC_SOURCES):
    fx = fluxes[src]
    E_src = np.linspace(max(fx.E_min + 1e-3, 1.8),
                        min(fx.E_max - 1e-3, 8.0), 300)
    phi_tot = fx(E_src, P, L)
    for iso in fx.isotopes:
        phi_iso = fx.isotope_flux(iso, E_src, P, L)
        frac = np.where(phi_tot > 0, phi_iso / phi_tot, np.nan)
        ax.plot(E_src, frac, label=iso, color=iso_colors[iso])
    ax.set_xlabel(r"$E_\nu$ (MeV)")
    ax.set_title(src)
    ax.set_ylim(0, 1)
    ax.legend(fontsize=7)

axes[0].set_ylabel(r"$f_i S_i / \sum_j f_j S_j$")
fig.suptitle(
    r"Per-isotope spectral fraction  ($f_i$ are constants; variation reflects different $S_i(E)$ shapes)",
    y=1.02
)
plt.tight_layout()
plt.show()

# %%
# Verify: sum of isotope_flux == total flux
phi_total_v = flux_e(E_range, P, L)
phi_parts_v = {iso: flux_e.isotope_flux(iso, E_range, P, L) for iso in flux_e.isotopes}
phi_sum = sum(phi_parts_v.values())
max_rel_err = np.max(np.abs(phi_sum - phi_total_v) / np.where(phi_total_v > 0, phi_total_v, 1))
print(f"Max relative error (sum of isotope_flux vs total): {max_rel_err:.2e}")

# %% [markdown]
# ## 5. Composite spectrum — `InterpolatedFlux.from_composite`
#
# `kopeikin2012` provides a pre-weighted total spectrum (Table 3, PAN 2012).
# Useful for sub-MeV and low-energy studies.

# %%
flux_k2012 = InterpolatedFlux.from_composite("kopeikin2012")
flux_k1999 = InterpolatedFlux.from_composite("kopeikin1999")
print(f"Kopeikin 2012 valid range: {flux_k2012.E_min:.3f} – {flux_k2012.E_max:.2f} MeV")
print(f"Kopeikin 1999 valid range: {flux_k1999.E_min:.3f} – {flux_k1999.E_max:.2f} MeV")

# %%
E_low = np.linspace(0.02, 1.4, 300)

phi_k2012 = flux_k2012(E_low, P, L)
phi_k1999 = flux_k1999(E_low, P, L)

fig, ax = plt.subplots(figsize=(7, 4))
ax.semilogy(E_low, phi_k2012, label="Kopeikin 2012")
ax.semilogy(E_low, phi_k1999, label="Kopeikin 1999", linestyle="--")
ax.set_xlabel(r"$E_\nu$ (MeV)")
ax.set_ylabel(r"$d\Phi/dE_\nu$  (cm$^{-2}$ MeV$^{-1}$)")
ax.set_title(f"Sub-MeV reactor flux  (P = {P} GW,  L = {L} m)")
ax.legend()
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 6. Reactor-specific total spectra
#
# Different reactor types (typical PWR, KSNPS, CONUS, Daya Bay) use different
# fission fraction mixes, shifting the spectral shape.  Here we show the effect
# for each isospec source.
#
# Each row is one data source; the left panel shows absolute flux for all
# presets, and the right panel shows the ratio to the "typical" PWR preset to
# highlight spectral shape differences.

# %%
presets = ["typical", "ksnps", "conus", "daya_bay"]
preset_colors = {"typical": "C0", "ksnps": "C1", "conus": "C2", "daya_bay": "C3"}
preset_ls     = {"typical": "-",  "ksnps": "--", "conus": "-.", "daya_bay": ":"}

E_test = np.linspace(2.0, 7.0, 300)

fig, axes = plt.subplots(len(ISOSPEC_SOURCES), 2,
                         figsize=(13, 3.5 * len(ISOSPEC_SOURCES)))

for row, src in enumerate(ISOSPEC_SOURCES):
    ax_left  = axes[row, 0]
    ax_right = axes[row, 1]

    phi_ref = None
    for preset in presets:
        fx = make_flux(src, fission_fractions=preset)
        E_src = np.linspace(max(fx.E_min + 1e-3, 2.0),
                            min(fx.E_max - 1e-3, 7.0), 300)
        phi = fx(E_src, P, L)
        ax_left.semilogy(E_src, phi, label=preset,
                         color=preset_colors[preset], ls=preset_ls[preset])
        if preset == "typical":
            phi_ref = (E_src, phi)

    ax_left.set_xlabel(r"$E_\nu$ (MeV)")
    ax_left.set_ylabel(r"$d\Phi/dE_\nu$  (cm$^{-2}$ MeV$^{-1}$)")
    ax_left.set_title(f"{src} — absolute flux  (P = {P} GW, L = {L} m)")
    ax_left.legend(fontsize=8)

    E_ref, phi_typical = phi_ref
    for preset in presets:
        if preset == "typical":
            continue
        fx = make_flux(src, fission_fractions=preset)
        phi = fx(E_ref, P, L)
        ratio = np.where(phi_typical > 0, phi / phi_typical, np.nan)
        ax_right.plot(E_ref, ratio, label=f"{preset} / typical",
                      color=preset_colors[preset], ls=preset_ls[preset])

    ax_right.axhline(1.0, color="grey", lw=0.8, ls=":")
    ax_right.set_xlabel(r"$E_\nu$ (MeV)")
    ax_right.set_ylabel("Ratio to 'typical'")
    ax_right.set_title(f"{src} — ratio to 'typical' preset")
    ax_right.legend(fontsize=8)

plt.tight_layout()
plt.show()
