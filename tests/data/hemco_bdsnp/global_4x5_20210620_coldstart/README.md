# HEMCO 3.12.1 global SoilNOx reference

This directory contains the independent gridded reference for CECE BDSNP
parity development.

## Reference identity

- HEMCO release: 3.12.1
- Upstream commit: `07da3c29fd85abc3824cb6288578b0b68c2395a3`
- Standalone deposited-N interface fix:
  `c62c250f9ab71f04759769d88100ae8877da0796`
- DRYCOEFF buffer fix:
  `0588005eae6d99bc7c7044f1706dc6ec8a52f33f`
- Ursa Slurm job: `19544461`
- Compiler: GNU Fortran 13.3.0
- NetCDF-C: 4.9.2
- NetCDF-Fortran: 4.6.1
- HEMCO precision: REAL8
- OpenMP: disabled

The two temporary patches repair standalone input plumbing. They do not
change the SoilNOx temperature, moisture, biome, canopy, pulsing,
fertilizer, deposited-N, or flux equations.

## Case

- Grid: official HEMCO global 4x5, 72 longitude × 46 latitude cells
- Time: 2021-06-20 12:00 to 13:00 UTC
- Emissions steps: one, 3600 seconds
- Meteorology: processed MERRA-2 4x5
- `UseSoilTemperature`: false
- `Use fertilizer NOx`: true
- Extension scaling: 1.0
- Cold start:
  - `GWET_PREV=0`
  - `PFACTOR=1`
  - `DRYPERIOD=0`
  - deposited-N reservoir initialized from `DepReservoirDefault.nc`
- New `DRY_TOTN` and `WET_TOTN` forcing: controlled zero

This is a deterministic cold-start engineering reference. It is suitable
for algorithm and restart parity, but is not a spun-up climatological
soil-NO product.

## Reference fields

From `hemco_soilnox_global_4x5_20210620_1200.nc`:

- `InvSOILNOX_NO`: total SoilNOx, kg NO m-2 s-1
- `EmisNO_Fert_a`: actual fertilizer plus deposited-N contribution
- Natural biome contribution:
  `InvSOILNOX_NO - EmisNO_Fert_a`

The explicit ext-104 `EmisNO_Fert` diagnostic duplicates the total
extension flux, so its name is misleading. Because that name was already
occupied, HEMCO registered its internal fertilizer plus deposited-N field
as `EmisNO_Fert_a`. CECE comparisons therefore use `InvSOILNOX_NO` for
total SoilNOx and `EmisNO_Fert_a` for the fertilizer plus deposited-N
contribution.

From `hemco_soilnox_restart_4x5_20210620_1300.nc`:

- `DEP_RESERVOIR`
- `GWET_PREV`
- `DRYPERIOD`
- `PFACTOR`

## Validated global totals

- Total SoilNOx: `6.29812993790932524e+02 kg NO s-1`
- Fertilizer plus deposited N:
  `1.69015832309265335e+02 kg NO s-1`
- Natural biome component:
  `4.60797161481667217e+02 kg NO s-1`

## Input hashes

- `DepReservoirDefault.nc`:
  `3289b982b9c5737da22693e9f28088f0df021b47ab6229b872c3329292b8397d`
- `soilNOx.climate.generic.05x05.nc`:
  `c48a2b06483609f67134796242ce29d66e967283a55c4f13a52a90a534c3cdfd`
- `soilNOx.fert_res.generic.05x05.nc`:
  `9de51b592f1da3dfbbb84ea11e62bd7c34013757cc9495117a60497dee1280c7`
- `soilNOx.landtype.generic.025x025.1L.nc`:
  `855088b7a77e602793ac410fa1cbef282921603b8ad6373b68c44cad22a58cac`
- MERRA-2 CN:
  `eafa5a37d19f4d8e26f5a29628f46e1ea3e94480e81c43410eb30b18e43c75c1`
- MERRA-2 A1:
  `5756866d478d91ece2040698cc3a4bd5b39a06e9e33b61c55e0402ef2237c0c5`

## Scope

This fixture is the HEMCO side of the gridded comparison. Passing the
existing scalar-vector tests alone does not establish gridded parity.
CECE must independently compute the 24-biome composition, environmental
terms, added-N contribution, state transition, and restart fields before
full parity can be claimed.
