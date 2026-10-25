# PCM (Sb2S3) — Thermal / Optical / Python

This repository contains the PCM simulation project. The workflow is split into three branches:

| Branch | Contents |
|---|---|
| `python` | Python workflow: HEAT → `T_pcm` → `Xf` (JMA) → `neff` (FDTD mode expansion) |
| `thermal` | Lumerical DEVICE thermal simulation (`PCM_Heat.ldev`) |
| `optical` | Lumerical FDTD optical simulation (`PCM.fsp`) |

## Workflow

1. **HEAT** reads the transient temperature of the PCM patch.
2. **JMA** converts `T_pcm(t)` to crystallinity `Xf(t)`.
3. **FDTD** updates the `PCM_mix` refractive index for selected `Xf` values.
4. A `ModeExpansionMonitor` at `x = 0` (inside the PCM patch) gives `neff` at 550 nm.
5. `neff(t)` is interpolated over the full time vector.

## Results at 550 nm

| Xf | n_PCM | k_PCM | Re(neff) | Im(neff) |
|---:|---:|---:|---:|---:|
| 0.0 | 2.8300 | 0.00000 | 1.381398 | 0.000345 |
| 0.5 | 3.0732 | 0.02012 | 1.426107 | 0.005944 |
| 1.0 | 3.3700 | 0.05000 | 1.524540 | 0.020181 |

`Re(neff)` is the effective phase index; `Im(neff)` is the attenuation/loss part.

## Branches

```bash
git clone https://github.com/mahmoud-diaa-eldeen/PCM.git
cd PCM
git branch -a
git checkout python
```

- `python`: `pcm.py`, `pcm_results.csv`, `pcm_results.png`, README.
- `thermal`: `PCM_Heat.ldev`, `PCM_Heat_p0.log`, README.
- `optical`: `PCM.fsp`, `PCM_p0.log`, README.

## Notes

- `RUN_FDTD = False` in `pcm.py` computes `neff` without a full time-domain FDTD run. Set it to `True` to also recompute `T_optical` (requires more RAM).
- The current JMA model describes crystallization; it does not model melt-quench amorphization.
- The optical mode is a PCM-loaded hybrid mode (the core index is lower than the cladding index).
