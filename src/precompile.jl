import PrecompileTools
if Base.VERSION >= v"1.9"
    PrecompileTools.@setup_workload begin
        # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
        # precompile file and potentially make loading faster.
        cutoff = 0.05
        x3d = rand(3, 10)
        y3d = rand(3, 10)
        u3d = [ 1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0 ]
        x2d = rand(2, 10)
        y2d = rand(2, 10)
        u2d = [ 1.0 0.0; 0.0 1.0 ]
        f(_,_,_,_,d2,out) = out += d2
        PrecompileTools.@compile_workload begin
            box = Box(u3d, cutoff)
            cl = CellList(x3d, box)
            box = Box(u2d, cutoff)
            cl = CellList(x2d, box)
            map_pairwise!(f, 0.0, box, cl)
            neighborlist(x3d, y3d, cutoff, unitcell=u3d)
            neighborlist(x3d, cutoff, unitcell=u3d)
            neighborlist(x3d, y3d, cutoff)
            neighborlist(x3d, cutoff)
            neighborlist(x2d, y2d, cutoff, unitcell=u2d)
            neighborlist(x2d, cutoff)
            neighborlist(x2d, y2d, cutoff)
            neighborlist(x2d, cutoff, unitcell=u2d)
        end
    end
end