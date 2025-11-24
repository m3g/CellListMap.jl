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

Version 0.9.15-DEV
--------------
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
- ![BUGFIX][badge-bugfix] When no PBCs are defined, and the position limits are smaller than the cutoff, with some pathological coordinates,  double-counting could occurr. Fixed.

Version 0.9.11
--------------
- ![INFO][badge-info] Fix typo in color definition of non-interface function `draw_computing_cell`. 
- ![INFO][badge-info] Reduce precompilation work by running smaller examples
- ![INFO][badge-info] Add tests to random cells with negative off-diagonal elements.

Previous versions
--------------
Previous version changes are listed in the [releases](https://github.com/m3g/CellListMap.jl/releases) page.

