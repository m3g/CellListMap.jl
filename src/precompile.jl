@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    if Base.VERSION >= v"1.9"
        cutoff = 0.05
        x3d = rand(3, 100)
        y3d = rand(3, 100)
        u3d = [ 1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0 ]
        x2d = rand(2, 100)
        y2d = rand(2, 100)
        u2d = [ 1.0 0.0; 0.0 1.0 ]
        x, box = CellListMap.xatomic(5000)
        f(x,y,i,j,d2,out) = out += d2
        @compile_workload begin
            cl = CellList(x,box)
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