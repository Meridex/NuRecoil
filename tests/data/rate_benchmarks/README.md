# Rate Benchmark Data Format

Put external CEvNS/EvES benchmark files in this directory with the extension
`.json`.  The test suite will skip benchmark comparisons when no `.json` files
are present.

Each JSON file describes one benchmark case.  Energies are in MeV unless the
field name says otherwise.

```json
{
  "name": "ge76_cevns_example",
  "observable": "cevns",
  "target": {"Z": 32, "A": 76},
  "normalization": {
    "target_mass_g": 1000.0,
    "P_GW": 3.0,
    "L_m": 30.0,
    "exposure_s": null
  },
  "flux": {
    "model": "phenomenological",
    "fission_fractions": "typical",
    "spectrum_model": "huber_mueller",
    "energy_per_fission": "ma_2013"
  },
  "cross_section": {
    "model": "sm_cevns",
    "flavor": "tree",
    "form_factor": "helm"
  },
  "detector": {
    "resolution": {"model": "cdex"},
    "quenching": {"model": "lindhard", "k": null}
  },
  "integration": {
    "E_nu_range_MeV": [1.8, 8.0]
  },
  "reference": {
    "unit": "counts_per_s_per_bin",
    "rtol": 0.02,
    "atol": 0.0,
    "bins": [
      {"low_MeV": 0.0, "high_MeV": 0.00005, "value": 0.0}
    ]
  }
}
```

For EvES, use:

```json
{
  "observable": "eves",
  "cross_section": {
    "model": "sm_eves",
    "flavor": "nu_e",
    "antineutrino": false
  },
  "detector": {
    "resolution": {"model": "cdex"},
    "quenching": {"model": "constant", "value": 1.0}
  }
}
```

Rules:

- `observable` must be `"cevns"` or `"eves"`.
- `reference.unit` must be `"counts_per_s_per_bin"` or `"counts_per_bin"`.
- If `reference.unit` is `"counts_per_bin"`, set `normalization.exposure_s`.
- `bins` must be ordered and contiguous enough for the intended comparison.
- `value` should be the external reference value for that detected-energy bin.
- Set `rtol` to the reference precision. For literature plots digitized by eye,
  start with a loose value such as `0.05`; for table-level data, use a tighter
  value such as `0.01` or below.
