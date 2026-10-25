# PCM Python workflow

`pcm.py` performs the full chain:

1. Reads the HEAT thermal result (`PCM_Heat.ldev`).
2. Computes crystallinity `Xf(t)` with the JMA model.
3. Updates `PCM_mix` n/k for the selected `Xf` targets.
4. Solves the FDTD mode-expansion monitor and extracts `neff` at 550 nm.
5. Saves `pcm_results.csv` and `pcm_results.png`.

## Files

- `pcm.py`
- `pcm_results.csv` — columns: `t_s, T_pcm_K, Xf, neff_real, neff_imag, T_optical`
- `pcm_results.png`

## Requirements

- Lumerical 2024 R1 (DEVICE + FDTD)
- Python 3, numpy, matplotlib
- `lumapi` from the Lumerical installation

## Run

```bash
python pcm.py
```

`RUN_FDTD = False` computes `neff` without a full time-domain FDTD run. Set it to `True` to also recompute `T_optical` (requires more RAM).

## Results at 550 nm

| Xf | n_PCM | k_PCM | Re(neff) | Im(neff) |
|---:|---:|---:|---:|---:|
| 0.0 | 2.8300 | 0.00000 | 1.381398 | 0.000345 |
| 0.5 | 3.0732 | 0.02012 | 1.426107 | 0.005944 |
| 1.0 | 3.3700 | 0.05000 | 1.524540 | 0.020181 |
