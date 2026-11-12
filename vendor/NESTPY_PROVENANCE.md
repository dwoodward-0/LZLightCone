# nestpy source provenance for the LZ-lightcone WIMP connector

This is the public source distribution already pinned by the project lockfile:

- Upstream package: `nestpy==2.1.3`.
- URL: https://files.pythonhosted.org/packages/4c/ed/e3b0465c42efee32946d08e632cb2c42693215d8189900c40c3daa4621aa/nestpy-2.1.3.tar.gz
- SHA256: `258155894377e24c02a763e8e58319e01fbd42f12aa7a8c960084848b1a84e84`
- Retrieval date: 2026-11-12.
- Local package version: `2.1.3+lz.wimp1`, identifying this source-only connector patch.

The sdist checksum was verified before extraction. Its upstream contents are
stored under `nestpy/` unchanged except for:

1. `nestpy/src/nestpy/spectra.cpp`: expose NEST's existing static
   `TestSpectra::WIMP_dRate`, read-only native preparation fields `integral`,
   `xMax`, and `divisor`, and `WIMP_with_metadata`, which performs one native
   preparation and returns its samples plus matching metadata.
2. `nestpy/pyproject.toml`: change only the local package version so a patched
   build cannot be mistaken for the upstream binary wheel or replaced by a
   newer unpatched release.

No file under `nestpy/lib/nest/` is modified. In particular the native WIMP
spectrum, Helm form factors, isotope sampling, halo parameters, reference
nucleon cross section of 1e-36 cm2, one-million-evaluation integral, Earth
velocity approximation, detector implementation and random-number generator
remain upstream code. Neither this patch nor the project uses NEST's reported
`0.0.0`/blank version strings as native-source revision information; the
verified archive checksum is the provenance record.

`NESTPY_WIMP_CONNECTOR.patch` records the patch relative to the verified
archive. The package remains subject to its upstream GPL-3.0 license,
retained at `nestpy/LICENSE`.
