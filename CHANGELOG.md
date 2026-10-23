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

Version 0.10.10-DEV
--------------
- ![INFO][badge-info] Requires Julia 1.11 or later. The `public` names are declared with the `public` keyword, and Compat is no longer a dependency.
- ![INFO][badge-info] The test files (NAMD and Gromacs simulations) are no longer part of the package installation: they are a lazy artifact (from [MolSimTestFiles](https://github.com/m3g/MolSimTestFiles)), downloaded when needed, and obtained with `CellListMap.test_file(name)`, which is available after `using LazyArtifacts` (a package extension). The `CellListMap.argon_pdb_file` constant used in the examples is replaced by `CellListMap.test_file("cubic.pdb")`.

Version 0.10.9
--------------
- ![ENHANCEMENT][badge-enhancement] `pairwise!` (and the rest of the `ParticleSystem` path) can be compiled with `juliac --trim`: `show_progress` defaults to `nothing` (`Union{Nothing,Bool}`), such that the progress meter code is not compiled unless requested (`show_progress = true`/`false` work as before), the cell list construction waits on its tasks with typed results instead of `@sync`, and the `getproperty(sys, ::Symbol)` and `getproperty(pair, ::Symbol)` methods are inlined.
- ![BUGFIX][badge-bugfix] With a custom `validate_coordinates` function, or with `validate_coordinates = nothing`, `NaN` coordinates no longer throw the error of the default validation while building the cell lists (`nothing` now disables the check, as documented).

Version 0.10.8
--------------
- ![ENHANCEMENT][badge-enhancement] Faster parallel cell list updates: the default validation of the coordinates (no `NaN` values) is done while the particles are located in the cells, which runs in parallel, instead of in a separate serial pass, and the copy of the coordinates in `update!` is run in parallel for large systems (32k particles or more). This reduces the time of `update!(sys; xpositions = x, rebuild = true)` in parallel by ~20-30% (periodic systems with 110k and 1M particles, 14 threads). As a consequence, in parallel, same-size `update!`s of large systems allocate the (few) task objects of the copy.

Version 0.10.7
--------------
- ![FEATURE][badge-feature] Unit cells can be of any size relative to the cutoff. The previous requirement that the distances between opposite faces of the unit cell be greater than twice the cutoff was removed. For smaller cells, pairs are reported once for each periodic image within the cutoff (thus the same pair `(i, j)` may appear more than once, and `i == j` pairs, of a particle with its own images, may appear). For cells satisfying the previous requirement, results are unchanged.
- ![FEATURE][badge-feature] The unit cell of a periodic system can be updated to a cell of any shape: an orthorhombic cell (given as a vector of sides) can be updated to a triclinic cell (given as a matrix), and vice versa. Previously, the type of unit cell could not be changed after the construction of the system.
- ![FEATURE][badge-feature] `NeighborPair` has a new `shift` field: the integer lattice translation of particle `j` relative to particle `i`, such that `y[pair.j] + unitcell * pair.shift - x[pair.i] ≈ pair.y - pair.x`. It identifies the image of each pair, and is always zero for non-periodic systems.
- ![ENHANCEMENT][badge-enhancement] New implementation of the cell lists: the cells are defined on a grid in fractional coordinates of the unit cell, each particle is stored once (no replication of periodic images), and periodic images are reached in the pairwise computation by wrapping the indices of neighboring cells. Orthorhombic, triclinic and non-periodic systems share the same implementation. Candidate pairs between cells with many particles are filtered by projection on the axis connecting the cell centers. Cell list construction is much faster (2-8x), while the speedup of complete applications depends on the cost of the mapped function and on the density of the system: relative to the previous version, the cell list part of Packmol.jl (100k water molecules) is ~45% faster, the computation of distribution functions in ComplexMixtures.jl (dense system, 20 Å cutoff) is ~5% faster, and a pairwise velocity-correlation computation (1M particles) is ~30% faster for a sparse system and has about the same speed for a dense system.
- ![ENHANCEMENT][badge-enhancement] `neighborlist` no longer copies the whole list of pairs when the number of pairs slightly exceeds the preallocated estimate (which happened in about half of the uniform systems), making it up to 3x faster in those cases.
- ![ENHANCEMENT][badge-enhancement] The threaded copies of the output are reset in parallel, each by the task that uses it in the `pairwise!` computation, instead of serially before the computation. This removes a cost that grew with the number of threads for large array outputs (e.g. 70 ms per call with 16 threads, and 133 ms with 32 threads, for the per-atom gradients of a 2.77M-atom system in Packmol.jl). Simple outputs are unaffected.
- ![INFO][badge-info] In two-set systems, only the cell lists of the sets whose positions changed are rebuilt (as in previous versions: this fixes a regression of the new cell list implementation, in which both lists were rebuilt when either set was updated). The list of the first set, which is only iterated over, is built at a cost independent of the number of cells when it has few particles relative to the number of cells (by sorting its particles by cell), making updates of a small first set (e.g. one molecule queried against a large fixed set) nearly free.
- ![INFO][badge-info] Degenerate (zero-volume) unit cells throw an `ArgumentError`.
- ![INFO][badge-info] The order of the pairs in the neighbor lists, and of the calls to the mapped function in `pairwise!`, changed. The `pair.x` coordinate of each pair is now always the coordinate of particle `i` wrapped into the unit cell (identical to the input coordinate if the particle is inside the cell), and `pair.y` the image of particle `j` within the cutoff of `pair.x`.
- ![INFO][badge-info] Internal changes: the `CellList` structure was rewritten, the `OrthorhombicCell` and `TriclinicCell` unit cell types were replaced by a single `PeriodicCell` type (non-periodic systems use `NonPeriodicCell`), the `Box` structure no longer contains the `aligned_unit_cell`, `rotation`, `cell_size` fields, and the `AuxThreaded` auxiliary structures were removed.

Version 0.10.6
--------------
- ![FEATURE][badge-feature] `output_name` can be given as a `Val` (e.g. `output_name = Val(:energy)`), in addition to a `Symbol`, making the `ParticleSystem` constructors fully inferable. Passing a `Symbol` works as before.
- ![ENHANCEMENT][badge-enhancement] `ParticleSystem` construction is type-stable when coordinates are vectors of `SVector`s, for any type of unit cell input (`Vector`, `Matrix`, static arrays, or `nothing`), allowing code using CellListMap to be compiled with `juliac --trim`.
- ![INFO][badge_info] Remove dependency on `DocStringExtensions`. 

Version 0.10.5
--------------
- ![FEATURE][badge-feature] Add `rebuild` keyword argument to `update!` to rebuild the cell lists immediately. 

Version 0.10.4
--------------
- ![ENHANCEMENT][badge-enhancement] Improve scaling of cell list construction phases.

Version 0.10.3
--------------
- ![ENHANCEMENT][badge-enhancement] Skip `inv_rotation` for `OrthorhombicCell`s. 
- ![INFO][badge-info] Clarifty the meaning of `NeighborPair` fields, specifically concerning the minimum-image consideration.
- ![INFO][badge-info] Document the possibility of modifying individual positions with `setindex!` and `SVector`s. 

Version 0.10.2
--------------
- ![FEATURE][badge-feature] Support empty arrays of coordinates for `unitcell=nothing`. 
- ![INFO][badge-info] Show methods of ParticleSystems now show the update status of the coordinates, expliciting the status of the cell lists (updated or not).
- ![BUGFIX][badge-bugfix] If the number of particles was reduced to less than the number of batches of parallel computations, stale data remained in auxiliary threading arrays. This manifested when testing removing all particles from the a system. 

Version 0.10.1
--------------
- ![ENHANCEMENT][badge-enhancement] Improve type stability of getter functions of ParticleSystems.
- ![INFO][badge-info] Update CI setup-julia to v3.

Version 0.10.0
--------------
- ![FEATURE][badge-feature] Automatic tracking of position mutations and updating of cell lists. 
- ![FEATURE][badge-feature] `pairwise!(; reset=[true(default)/false]`: the optional `reset` keyword of `pairwise!`, when set to `false`, avoids resetting the initial value of `output` to `zero(typeof(output))`.
- ![FEATURE][badge-feature] `update!(sys::AbstractParticleSystem;)` as cleaner and more convenient way to update system properties. 
- ![FEATURE][badge-feature] `reduce_output!` is exposed API for advanced custom reductions. 
- ![ENHANCEMENT][badge-enhancement] Improve scaling of cell list construction and updating. 
- ![ENHANCEMENT][badge-enhancement] Improve performance of inner loops.
- ![ENHANCEMENT][badge-enhancement] Improve performance of neighborlist reduction and cold start of lists.
- ![ENHANCEMENT][badge-enhancement] Update list of unchanched coordinates of `NonPeriodicCell` only if the limits expanded. 
- ![ENHANCEMENT][badge-enhancement] Split computation of `NonPeriodicCell`s into its own path, which then can be simplified. Should display better performance and scalling, particularly in the cell list construction step.
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
- ![INFO][badge-info] The internal representation of coordinantes is done with `ParticleSystemPositions` array type - which is not a subtype of abstract array. 
- ![INFO][badge-info] Add performance test.
- ![INFO][badge-info] Update examples.
- ![INFO][badge-info] Default name for output in ParticleSystem is `default_output_name`, but it can be accessed though `sys.output`, as before.
- ![INFO][badge-info] Internal fields of ParticleSystems are not part of `PrivateParticleSystemData` and are accessed through internal getter functions, or through `sys.private.field`. 
- ![BUGFIX][badge-bugfix] Fix cross-computation (`pairwise!(f, x, sys)`) with `NonPeriodicCell` failing to find pairs when particle coordinates span negative values or large coordinate ranges / do not modify input coordinates in `NonPeriodicCell`. This bug was never released, it was created and fixed in the development version.
- ![BUGFIX][badge-bugfix] If the number of particles is reduced by updating and becomes smaller than nbatches, there was a crash. Fixed. Also fixed automatic updating of nbatches when first set of positions change (bug never released).
- ![BUGFIX][badge-bugfix] Fix `nbatches` not being updated for the `:map` phase of `CellListPair` systems when the particle count changes (bug never released).

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

