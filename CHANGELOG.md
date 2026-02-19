CellListMap.jl Changelog
===========================
  
[badge-breaking]: https://img.shields.io/badge/BREAKING-red.svg
[badge-deprecation]: https://img.shields.io/badge/Deprecation-orange.svg
[badge-feature]: https://img.shields.io/badge/Feature-green.svg
[badge-experimental]: https://img.shields.io/badge/Experimental-yellow.svg
[badge-enhancement]: https://img.shields.io/badge/Enhancement-blue.svg
[badge-bugfix]: https://img.shields.io/badge/Bugfix-purple.svg
[badge-fix]: https://img.shields.io/badge/Fix-purple.svg
[badge-info]: https://img.shields.io/badge/Info-gray.svg

Version 0.10.0-DEV
--------------
- ![FEATURE][badge-feature] `pairwise!(f, x, sys)` maps `f` to the pairs of the particles defined in vector `x` to the cell list defined in `sys` (to reuse the same cell list for cross-computations with different set of particles).
- ![FEATURE][badge-feature] `pairwise!(; reset=[true(default)/false]`: the optional `reset` keyword of `pairwise!`, when set to `false`, avoids resetting the initial value of `output` to `zero(typeof(output))`.
- ![FEATURE][badge-feature] Automatic update tracking of coordinates and cell list updating using the ParticleSystemPositions internal coordinate representation. Views share the reference updated flag, such that mutations to views are tracked.
- ![ENHANCEMENT][badge-enhancement] Improve scaling of cell list construction and updating. 
- ![ENHANCEMENT][badge-enhancement] Improve performance of inner loops.
- ![ENHANCEMENT][badge-enhancement] Improve performance of neighborlist reduction.
- ![BREAKING][badge-breaking] `map_pairwise!` was renamed to `pairwise!`. 
- ![BREAKING][badge-breaking] `pairwise`, without the `!` was removed, to stress the fact that the function always mutates the `output` field of the `ParticleSystem` object. 
- ![BREAKING][badge-breaking] The initial value of `output` is not set to `zero(typeof(output))` by default, but retains the given value. Resetting occurs on the call to `pairwise!` by default, and can be skipped with `reset=false`. 
- ![BREAKING][badge-breaking] The previous "lower level interface" is now internal and not exported or public anymore.
- ![BREAKING][badge-breaking] The `autoswap` function was removed, in favor of the new `pairwise!(f, x, sys)` method. It not available anymore in for `neighborlist` functions. 
- ![BREAKING][badge-breaking] The python interface was currently discontinued.
- ![BREAKING][badge-breaking] `validate_coordinates` has to be a `Function` (`(x)  -> nothing` to skip all validation).
- ![INFO][badge-info] `AbstractParticleSystem`, `ParticleSystem1` and `ParticleSystem2` are documented as public, to control dispatch in advanced applications.
- ![INFO][badge-info] Renamed internal methods of `pairwise!` to `_pairwise!` to improve error messages. 
- ![INFO][badge-info] Remove custom linear algebra code, use Base methods.
- ![INFO][badge-info] Remove deprecated code, improve code coverage.
- ![INFO][badge-info] Remove internal coordinate swapping by size.
- ![INFO][badge-info] Internal implementation of neighborlist moved to ParticleSystem interface.
- ![INFO][badge-info] The internal representation of coordinantes is done with `ParticleSytemPositions` array type - which is not a subtype of abstract array. 
- ![BUGFIX][badge-bugfix] Fix cross-computation (`pairwise!(f, x, sys)`) with `NonPeriodicCell` failing to find pairs when particle coordinates span negative values or large coordinate ranges / do not modify input coordinates in `NonPeriodicCell`. This bug was never released, it was created and fixed in the development version.

Version 0.9.17
--------------
- ![ENHANCEMENT][badge-enhancement] Improve scaling when two sets of particles are used.
- ![ENHANCEMENT][badge-enhancement] Automatic updating of the number of batches if the number of particles change.
- ![INFO][badge-info] Some code simplification and add comments to code. Fix documentation typos.
- ![INFO][badge-info] Drop support for Julia 1.9 (requires 1.10)
- ![INFO][badge-info] The `autoswap` option was deprecated, and does not have any effect now. Kept for backwards compatibility. 

Version 0.9.16
--------------
- ![INFO][badge-info] This version was broken and was replaced by 0.9.17.

Version 0.9.15
--------------
- ![INFO][badge-info] Better error messages and documentation of unitcell requirements.
- ![INFO][badge-info] Use julia-actions/cache@v2 in CI

Version 0.9.14
--------------
- ![INFO][badge-info] Document the fact that the order of the pairs output  by `neighborlist` functions is not guaranteed.

Version 0.9.13
--------------
- ![INFO][badge-info] Document the fact that the unitcell matrix is column-major, thus columns are the lattice vectors.

Version 0.9.12
--------------
- ![ENHANCEMENT][badge-enhancement] Improve type propagation when no PBCs are used, fixing Float32 to Float64 conversion of sides when `unitcell=nothing`.
- ![BUGFIX][badge-bugfix] When no PBCs are defined, and the position limits are smaller than the cutoff, with some pathological coordinates,  double-counting could occur. Fixed.

Version 0.9.11
--------------
- ![INFO][badge-info] Fix typo in color definition of non-interface function `draw_computing_cell`. 
- ![INFO][badge-info] Reduce precompilation work by running smaller examples
- ![INFO][badge-info] Add tests to random cells with negative off-diagonal elements.

Previous versions
--------------
Previous version changes are listed in the [releases](https://github.com/m3g/CellListMap.jl/releases) page.

