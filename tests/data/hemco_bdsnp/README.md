# HEMCO SoilNOx scalar reference vectors

`hemco_3_12_1_soilnox_reference.csv` is an independently generated oracle for
the stateless scalar terms in HEMCO SoilNOx extension 104. It is not CECE
output and is not a full gridded HEMCO run.

## Provenance

- Repository: `https://github.com/geoschem/HEMCO`
- Release: `3.12.1`
- Commit: `07da3c29fd85abc3824cb6288578b0b68c2395a3`
- Science source: `src/Extensions/hcox_soilnox_mod.F90`
- Science-source Git blob: `a8712e3e89fc0034ccd9761b476ce07608ac4506`
- Pristine science-source SHA-256:
  `02196771880b6fdcb2cc0a81d4c18c730fec3f6a1f7a5adfe168e214ceb58270`
- Reference-generator patch SHA-256:
  `f1547e1f4a1db1f0546d5294c7f2d4676aa77ccc36857f39d5218087f0393462`
- Reference-driver SHA-256:
  `ce1c63aea184ccb283b093256a1776a0aac5d7d6be11c52d1349afe76c70ffa3`
- Reference executable SHA-256:
  `beda237d335d67c83436a78d735a27d6d0f25ef53b387a863d16eb8ac542038e`
- CSV SHA-256:
  `d44f065e19e446ceb9b0a92b3b63c736c25cd0b83d9246a552eb86734df7a783`
- Compiler: GNU Fortran 13.3.0
- Precision and threading: `USE_REAL8=ON`, `OMP=OFF`

The temporary generator exposed the original private `SoilTemp`, `SoilWet`,
`SoilCrf`, and `FertAdd` routines through test-only wrappers without modifying
their bodies. The generator reconstructed `unitconv` and `cell_flux` from
HEMCO's production expression, using the pinned `FERT_SCALE=0.0068`, pulse=1,
and land fraction=1. Two executions produced byte-identical CSV files.

## Coverage and limits

The fixture contains all 24 biome constants and 58 cases covering:

- air- and soil-temperature branches, boundaries, and caps;
- arid and non-arid moisture responses;
- canopy-reduction arithmetic;
- fertilizer/deposition arithmetic;
- unit conversion; and
- one-hot contributions for all 24 biomes.

The value 0.0068 is the fertilizer scale in this pinned HEMCO release; the
fixture does not establish it as a universal value for other configurations or
versions. Each cell contribution fixes the wetting pulse and land fraction to
one. This fixture does not validate wetting-pulse state/restart, deposited-N
reservoir evolution, land/water/snow masks, canopy-resistance generation,
meteorological I/O, regridding, temporal diagnostics, or a full HEMCO
standalone run.
