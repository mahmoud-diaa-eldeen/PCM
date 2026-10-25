# PCM optical simulation

`PCM.fsp` is the Lumerical FDTD optical simulation.

## Structure

- `Cladding` — SiO2
- `Core` — object-defined dielectric
- `PCM_Patch` — `PCM_mix` (n,k material)
- `Tungsten Heater` — TiN

## Monitors

- `T` — 2D X-normal DFT monitor, 400–700 nm.
- `mode_expansion` — ModeExpansionMonitor at `x = 0` inside `PCM_Patch`,
  2D X-normal, fundamental mode, **550 nm**, 1 frequency point.

## Files

- `PCM.fsp`
- `PCM_p0.log`

The mode expansion monitor is used to extract `neff` for each crystallinity value.
