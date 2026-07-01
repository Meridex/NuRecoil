r"""
TabulatedFlux — load a pre-computed absolute antineutrino flux from a CSV file
and scale it to arbitrary (P, L) at call time.

File format
-----------
Line 0 (header / comment):
    ``// ..., for {P_ref} GW at {L_ref} m``
    The parser extracts the two numbers from this line using a regex.

Remaining lines (whitespace-separated, first column is an integer index):
    ``<index>  <E_MeV>  <flux_value>``
    where flux_value is in  # / MeV / cm² / s  at (P_ref, L_ref).

Scaling
-------
At call time the stored spectrum is rescaled as::

    Φ(E; P, L) = Φ_table(E) × (P / P_ref) × (L_ref / L)²

so the object can be reused for any reactor power and baseline.

Units
-----
    E_nu   : MeV
    P      : GW
    L      : m
    output : # / MeV / cm² / s
"""

from __future__ import annotations

import re
from importlib.resources import files
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import PchipInterpolator

from nurecoil.flux.base import FluxBase

# Regex to pull "P GW at L m" out of the header comment.
# Matches e.g.: "for 4.6 GW at 30 m" or "For 3.0GW at 500.0m"
_HEADER_RE = re.compile(
    r"for\s+([\d.eE+\-]+)\s*GW\s+at\s+([\d.eE+\-]+)\s*m",
    re.IGNORECASE,
)

# Path to the built-in data directory (distributed with the package).
_DATA_DIR = files("nurecoil.flux") / "data"
_DEFAULT_SOURCE = "reactor_neutrino_flux_example.csv"


def _parse_csv(path: Path) -> tuple[NDArray, NDArray, float, float]:
    """
    Parse a TabulatedFlux CSV file.

    Returns
    -------
    E_table : NDArray
        Neutrino energy grid [MeV], strictly increasing.
    flux_table : NDArray
        Flux values [# / MeV / cm² / s] at (P_ref, L_ref).
    P_ref : float
        Reference reactor power [GW].
    L_ref : float
        Reference baseline [m].
    """
    with open(path, encoding="utf-8") as fh:
        lines = fh.readlines()

    if not lines:
        raise ValueError(f"TabulatedFlux: file is empty: {path}")

    # --- parse header ---
    header = lines[0].strip()
    m = _HEADER_RE.search(header)
    if m is None:
        raise ValueError(
            f"TabulatedFlux: cannot parse reference (P, L) from header line:\n"
            f"  {header!r}\n"
            f"Expected format: '// ..., for <P> GW at <L> m'"
        )
    P_ref = float(m.group(1))
    L_ref = float(m.group(2))

    # --- parse data rows ---
    E_list: list[float] = []
    flux_list: list[float] = []
    for lineno, line in enumerate(lines[1:], start=2):
        line = line.strip()
        if not line or line.startswith("//") or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 3:
            raise ValueError(
                f"TabulatedFlux: expected at least 3 columns on line {lineno}, "
                f"got: {line!r}"
            )
        # columns: index  E_MeV  flux_value
        E_list.append(float(parts[1]))
        flux_list.append(float(parts[2]))

    if len(E_list) < 2:
        raise ValueError(
            f"TabulatedFlux: need at least 2 data rows, got {len(E_list)}"
        )

    E_table = np.asarray(E_list, dtype=float)
    flux_table = np.asarray(flux_list, dtype=float)

    if not np.all(np.diff(E_table) > 0):
        raise ValueError("TabulatedFlux: energy column must be strictly increasing.")

    return E_table, flux_table, P_ref, L_ref


class TabulatedFlux(FluxBase):
    r"""
    Antineutrino flux loaded from a pre-computed CSV table.

    The CSV stores $d\Phi/dE_\nu$ [# / MeV / cm² / s] at a reference power
    $P_{\rm ref}$ [GW] and baseline $L_{\rm ref}$ [m] encoded in the header
    comment.  At evaluation time the spectrum is rescaled:

    .. math::

        \Phi(E_\nu;\, P, L)
        = \Phi_{\rm table}(E_\nu)
          \times \frac{P}{P_{\rm ref}}
          \times \left(\frac{L_{\rm ref}}{L}\right)^2

    Interpolation between table points uses PCHIP (shape-preserving cubic
    Hermite), which avoids the oscillations of cubic splines for reactor
    spectra that span several orders of magnitude.

    Parameters
    ----------
    source : str or Path, optional
        Path to a CSV file, **or** the bare filename of a file bundled in
        ``nurecoil/flux/data/`` (without directory).
        Defaults to the built-in ``reactor_neutrino_flux_example.csv``.

    Examples
    --------
    >>> flux = TabulatedFlux()                          # built-in example
    >>> flux = TabulatedFlux("my_reactor_flux.csv")     # absolute/relative path
    >>> phi  = flux(np.array([2.0, 4.0, 6.0]), P=3.0, L=500.0)
    """

    # E_min / E_max are set from the table in __init__.
    E_min: float = 0.0
    E_max: float = 0.0

    def __init__(self, source: str | Path | None = None) -> None:
        path = self._resolve_path(source)
        E_table, flux_table, P_ref, L_ref = _parse_csv(path)

        self._csv_path = path
        self.E_min = float(E_table[0])
        self.E_max = float(E_table[-1])
        self._P_ref = P_ref
        self._L_ref = L_ref
        self._interp = PchipInterpolator(E_table, flux_table, extrapolate=False)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _resolve_path(source: str | Path | None) -> Path:
        """Return an absolute Path for *source*, searching the data dir."""
        if source is None:
            # Use the built-in default file.
            ref = _DATA_DIR / _DEFAULT_SOURCE
            # importlib.resources traversable → convert to real path
            return Path(str(ref))

        p = Path(source)
        if p.is_absolute() or p.exists():
            return p.resolve()

        # Treat as a bare filename inside the bundled data directory.
        candidate = Path(str(_DATA_DIR / str(source)))
        if candidate.exists():
            return candidate

        raise FileNotFoundError(
            f"TabulatedFlux: cannot find '{source}'. "
            f"Pass an absolute path or a filename inside nurecoil/flux/data/."
        )

    # ------------------------------------------------------------------
    # FluxBase implementation
    # ------------------------------------------------------------------

    def _flux(self, E_nu: NDArray, P: float, L: float) -> NDArray:
        scale = (P / self._P_ref) * (self._L_ref / L) ** 2
        interpolated = self._interp(E_nu)
        # PCHIP with extrapolate=False returns NaN outside range; set to 0.
        result = np.where(np.isfinite(interpolated), interpolated, 0.0)
        return result * scale
