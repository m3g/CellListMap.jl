import PrecompileTools
PrecompileTools.@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    import LinearAlgebra
    cutoff = 0.05
    x3d = rand(3, 3)
    x3d[:, 3] .= x3d[:, 2] .+ 1e-5
    y3d = copy(x3d)
    u3d = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    x2d = rand(2, 3)
    x2d[:, 3] .= x2d[:, 2] .+ 1e-5
    y2d = copy(x2d)
    u2d = [1.0 0.0; 0.0 1.0]
    x, box = CellListMap.xatomic(5000)
    f(pair, out) = out += pair.d2
    PrecompileTools.@compile_workload begin
        neighborlist(x3d, y3d, cutoff, unitcell=u3d)
        neighborlist(x3d, y3d, cutoff, unitcell=LinearAlgebra.diag(u3d))
        neighborlist(x3d, cutoff, unitcell=u3d)
        neighborlist(x3d, y3d, cutoff)
        neighborlist(x3d, cutoff)
        neighborlist(x2d, y2d, cutoff, unitcell=u2d)
        neighborlist(x2d, y2d, cutoff, unitcell=LinearAlgebra.diag(u2d))
        neighborlist(x2d, cutoff)
        neighborlist(x2d, y2d, cutoff)
        neighborlist(x2d, cutoff, unitcell=u2d)
        # 3D
        sys = ParticleSystem(positions=x3d, unitcell=u3d, cutoff=cutoff, output=0.0)
        map_pairwise!(f, sys)
        map_pairwise!(f, y3d, sys)
        sys = ParticleSystem(xpositions=x3d, ypositions=y3d, unitcell=u3d, cutoff=cutoff, output=0.0)
        map_pairwise!(f, sys)
        # 2D
        sys = ParticleSystem(positions=x2d, unitcell=u2d, cutoff=cutoff, output=0.0)
        map_pairwise!(f, sys)
        map_pairwise!(f, y2d, sys)
        sys = ParticleSystem(xpositions=x2d, ypositions=y2d, unitcell=u2d, cutoff=cutoff, output=0.0)
        map_pairwise!(f, sys)
    end
end