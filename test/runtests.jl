using TestItemRunner
@run_package_tests filter=t -> !(:performance in t.tags)

# Aqua
@testitem "Aqua.test_all" begin
    import Aqua
    Aqua.test_all(CellListMap)
end