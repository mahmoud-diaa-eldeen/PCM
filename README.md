# PCM thermal simulation

`PCM_Heat.ldev` is the Lumerical DEVICE thermal simulation.

## Settings

- Heat source: `power` boundary condition transient table.
- Power: **0.02 W**.
- Pulse: **200 ns on**, **10 ns fall**, **400 ns total**.
- `spatial results = "T"`.
- Transient minimum time step: **5 ns**.

## Files

- `PCM_Heat.ldev`
- `PCM_Heat_p0.log`

The result reaches `Tmax ≈ 883 K`, which drives `Xf` from 0 to 1 without melting (`Tm = 890 K`).
