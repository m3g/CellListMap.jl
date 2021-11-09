using CellListMap
using ForwardDiff
using Unitful
using Measurements
using LinearAlgebra: norm_sqr

#
# This function computes the sum of the squared distances of the particles
#
function sumsq(x,sides,cutoff;parallel=false)
    box = Box(sides,cutoff)
    cl = CellList(x,box,parallel=parallel)
    s = zero(eltype(x[1]^2))
    s = map_pairwise!(
        (x,y,i,j,d2,s) -> s += d2,
        s, box, cl, parallel=parallel
    )
    return s
end

#
# This is the analytical gradient of the above function
#
function sumsq_grad(x,sides,cutoff;parallel=false)
    g = similar(x)
    box = Box(sides,cutoff)
    cl = CellList(x,box)
    g .= zero(eltype(x))
    map_pairwise!(
        (x,y,i,j,d2,g) -> begin
            dx = x - y
            g[:,i] .+= 2*dx
            g[:,j] .-= 2*dx
            return g
        end,
        g, box, cl, parallel=parallel
    )
    return g
end

#
# Propagating measurements (or other more general types)
#
function sumsq_measurements(x_input,sides,cutoff;parallel=false)

    # The dual type of measurement does not propagate well on
    # all operations of CellListMap, and also comes with a signficant
    # cost, because they are mutable types. Therefore, it is better
    # to bypass the internal computations by closing over the 
    # input array and computing the result directly from the original
    # data. Here, `x_input` is the original array of coordinates
    # with uncertainties.  
    function sum_sqr_pair(i,j,s,x_input)
        d2 = norm_sqr(x_input[i] - x_input[j])
        s += d2
        return s
    end

    # Copy input matrix of measurements to matrix of floats
    x = getproperty.(x_input,:val)  
    
    # Build cell lists
    box = Box(sides,cutoff)
    cl = CellList(x,box)

    # The initial result is initialized with the appropriate type.
    s = measurement(0.,0.)

    # And instead of using the `x` and `y` coordinates provided by the
    # interface, we close over the `x_input` for the calculations.
    s = map_pairwise!(
        (x,y,i,j,d2,s) -> sum_sqr_pair(i,j,s,x_input),
        s, box, cl, parallel=parallel
    )
    return s
end

#
# Using generic types
#
function generic_types(iprint=true;parallel=false)

    # Data
    x = rand(3,1000)
    cutoff = 0.1 
    sides = [1.0, 1.0, 1.0]

    # Compare analtical and finite-difference gradient. Note that we convert the `sides`
    # and `cutoff` variables to the type of variable in `x`, to allow the propagation of the 
    # dual numbers required to ForwardDiff and Unitful
    function sumsq_generic(x)
        cutoff_generic = eltype(x).(cutoff) # cutoff is closed over
        sides_generic = eltype(x).(sides) # sides is closed over
        return sumsq(x, sides_generic, cutoff_generic)
    end
    check_grad = sumsq_grad(x,sides,cutoff) ≈ ForwardDiff.gradient(sumsq_generic, x)
    iprint && println("Gradients are correct: ", check_grad)

    # Now let us propagate units
    x = rand(3,1000)u"nm"
    cutoff = 0.1u"nm"
    sides = [1.0, 1.0, 1.0]u"nm"
    runit = sumsq(x,sides,cutoff)
    iprint && println("Result with units: ", runit)

    # Measurements, more general types, and optimal peformance
    x = [ measurement(rand(),0.01*rand()) for i in 1:3, j in 1:1000 ]
    cutoff = 0.1
    sides = [1.0, 1.0, 1.0]
    rmeasurement = sumsq_measurements(x,sides,cutoff,parallel=parallel)
    iprint && println("Result with units: ", rmeasurement)

    return check_grad, unit(runit), typeof(rmeasurement)
end
