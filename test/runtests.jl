using TestItemRunner
@run_package_tests

# Aqua
@testitem "Aqua.test_all" begin
    import Aqua
    Aqua.test_all(CellListMap)
end