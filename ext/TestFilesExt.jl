module TestFilesExt

using LazyArtifacts: LazyArtifacts, @artifact_str
import CellListMap: test_file

test_file(name::AbstractString) = joinpath(artifact"celllistmap_test_files", name)

end
