"""
Load per-actinide or composite antineutrino spectra from the bundled data files.

File naming convention
----------------------
    {source}_isospec_{isotope}.csv      — per-actinide spectra
    {source}_composite_total.csv        — pre-weighted composite spectra

where *source* is the lower-case citation key (e.g. ``"estienne2019"``) and
*isotope* is one of ``"u235"``, ``"u238"``, ``"pu239"``, ``"pu241"``.

CSV format
----------
Lines starting with ``#`` are comments.  The first non-comment line is the
column-name header.  Data lines use comma+space as delimiter; columns are
``E_nu`` [MeV] and ``dN_dE`` [# / MeV / fission].  Optionally a third column
``dN_dE_unc`` may be present but is ignored by this loader.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from numpy.typing import NDArray

# --------------------------------------------------------------------------- #
# Internal helpers                                                             #
# --------------------------------------------------------------------------- #

_DATA_DIR = Path(__file__).parent / "data"

# Canonical isotope-name normalisation: accept any reasonable spelling and
# map to the lower-case file suffix used in the CSV filenames.
_ISO_ALIASES: dict[str, str] = {
    # U-235
    "u235": "u235", "U235": "u235", "u-235": "u235", "U-235": "u235",
    # U-238
    "u238": "u238", "U238": "u238", "u-238": "u238", "U-238": "u238",
    # Pu-239
    "pu239": "pu239", "Pu239": "pu239", "PU239": "pu239",
    "pu-239": "pu239", "Pu-239": "pu239",
    # Pu-241
    "pu241": "pu241", "Pu241": "pu241", "PU241": "pu241",
    "pu-241": "pu241", "Pu-241": "pu241",
}

# Upper-case API keys (used by ReactorMix / IsotopeInterpolatedFlux)
_ISO_TO_API_KEY: dict[str, str] = {
    "u235": "U235",
    "u238": "U238",
    "pu239": "Pu239",
    "pu241": "Pu241",
}


def _normalise_isotope(isotope: str) -> str:
    """Return the lower-case file-suffix form, e.g. ``'U235'`` → ``'u235'``."""
    key = _ISO_ALIASES.get(isotope)
    if key is None:
        raise ValueError(
            f"Unknown isotope '{isotope}'.  "
            f"Accepted values: {sorted(set(_ISO_ALIASES.values()))}"
        )
    return key


def _read_csv(path: Path) -> tuple[NDArray, NDArray]:
    """Read a two-column (E_nu, dN_dE) CSV, skipping comment and header rows.

    Duplicate energy points (e.g. from digitised figures) are resolved by
    averaging the corresponding dN_dE values before returning.
    """
    E_list, S_list = [], []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            # Skip the column-name header row (first token not a float)
            tokens = s.split(",")
            try:
                E_list.append(float(tokens[0]))
                S_list.append(float(tokens[1]))
            except ValueError:
                continue  # header row
    if not E_list:
        raise ValueError(f"No numeric data found in {path}")

    E_arr = np.asarray(E_list, dtype=float)
    S_arr = np.asarray(S_list, dtype=float)

    # --- resolve duplicate energy values by averaging dN_dE ---
    if not np.all(np.diff(E_arr) > 0):
        unique_E, inverse = np.unique(E_arr, return_inverse=True)
        unique_S = np.zeros_like(unique_E)
        counts = np.zeros(len(unique_E), dtype=int)
        for i, idx in enumerate(inverse):
            unique_S[idx] += S_arr[i]
            counts[idx] += 1
        E_arr = unique_E
        S_arr = unique_S / counts

    return E_arr, S_arr


# --------------------------------------------------------------------------- #
# Public catalogue helpers                                                     #
# --------------------------------------------------------------------------- #

def list_sources() -> list[str]:
    """Return the unique source keys available in the bundled data directory."""
    sources = set()
    for p in _DATA_DIR.glob("*.csv"):
        sources.add(p.stem.split("_")[0])
    return sorted(sources)


def list_isotopes(source: str) -> list[str]:
    """
    Return the API-key isotope names available for *source*.

    Returns an empty list if *source* only provides a composite spectrum.
    """
    result = []
    for p in _DATA_DIR.glob(f"{source}_isospec_*.csv"):
        suffix = p.stem.split("_isospec_")[-1]   # e.g. "u235"
        result.append(_ISO_TO_API_KEY.get(suffix, suffix))
    return sorted(result)


# --------------------------------------------------------------------------- #
# Core loader                                                                  #
# --------------------------------------------------------------------------- #

def load_spectrum(
    source: str,
    isotope: str | None = None,
) -> tuple[NDArray, NDArray]:
    """
    Load a bundled antineutrino spectrum table.

    Parameters
    ----------
    source : str
        Citation key matching the CSV filename prefix, e.g. ``"estienne2019"``,
        ``"mueller2011"``, ``"kopeikin2012"``.  Case-sensitive.
    isotope : str, optional
        Actinide name, e.g. ``"U235"``, ``"Pu239"``.  Required for
        per-actinide (``isospec``) datasets; omit for composite datasets.
        Case-insensitive (``"u235"``, ``"U-235"`` etc. all work).

    Returns
    -------
    E_nu : ndarray, shape (n,)
        Neutrino energy grid [MeV], strictly increasing.
    dN_dE : ndarray, shape (n,)
        Spectrum per fission [# / MeV / fission].

    Raises
    ------
    ValueError
        If the requested source / isotope combination is not found.

    Examples
    --------
    >>> E, S = load_spectrum("estienne2019", "U235")
    >>> E, S = load_spectrum("kopeikin2012")
    """
    if isotope is not None:
        iso_key = _normalise_isotope(isotope)
        path = _DATA_DIR / f"{source}_isospec_{iso_key}.csv"
        if not path.exists():
            available = [p.name for p in _DATA_DIR.glob(f"{source}_isospec_*.csv")]
            raise ValueError(
                f"No isospec file for source='{source}', isotope='{isotope}'.\n"
                f"Available files: {available if available else '(none)'}"
            )
    else:
        # Try composite first, then fall back to isospec wildcard
        path = _DATA_DIR / f"{source}_composite_total.csv"
        if not path.exists():
            candidates = list(_DATA_DIR.glob(f"{source}_*.csv"))
            raise ValueError(
                f"No composite file for source='{source}', and no isotope was "
                f"specified.\nAvailable files for this source: "
                f"{[p.name for p in candidates] if candidates else '(none)'}\n"
                f"Hint: pass isotope='U235' etc. to load a per-actinide spectrum."
            )

    return _read_csv(path)


def load_all_isotopes(source: str) -> dict[str, tuple[NDArray, NDArray]]:
    """
    Load all available per-actinide spectra for *source*.

    Returns
    -------
    dict mapping API-key isotope name (e.g. ``"U235"``) to ``(E_nu, dN_dE)``.

    Raises
    ------
    ValueError
        If no isospec files are found for *source*.

    Examples
    --------
    >>> spectra = load_all_isotopes("estienne2019")
    >>> flux = IsotopeInterpolatedFlux(spectra)
    """
    paths = sorted(_DATA_DIR.glob(f"{source}_isospec_*.csv"))
    if not paths:
        raise ValueError(
            f"No isospec files found for source='{source}'.\n"
            f"Available sources: {list_sources()}"
        )
    result: dict[str, tuple[NDArray, NDArray]] = {}
    for p in paths:
        suffix = p.stem.split("_isospec_")[-1]
        api_key = _ISO_TO_API_KEY.get(suffix, suffix)
        result[api_key] = _read_csv(p)
    return result
