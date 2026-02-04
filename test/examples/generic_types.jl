using CellListMap
using ForwardDiff
using Unitful
using Measurements
using StaticArrays
using LinearAlgebra: norm_sqr

#
# This function computes the sum of the squared distances of the particles
#
function sumsq(x, sides, cutoff; parallel = false)
    sys = ParticleSystem(
        positions = x,
        unitcell = sides,
        cutoff = cutoff,
        parallel = parallel,
        output = zero(typeof(cutoff))^2,
    )
    s = pairwise!((pair, s) -> s += pair.d2, sys)
    return s
end

#
# This is the analytical gradient of the above function
#
function sumsq_grad(x, sides, cutoff; parallel = false)
    sys = ParticleSystem(
        positions = x,
        unitcell = sides,
        cutoff = cutoff,
        parallel = parallel,
        output = similar(x),
    )
    g = pairwise!(
        (pair, g) -> begin
            (; i, j, x, y) = pair
            dx = x - y
            g[:, i] .+= 2 * dx
            g[:, j] .-= 2 * dx
            return g
        end,
        sys
    )
    return g
end

#
# Propagating measurements (or other more general types)
#
function sumsq_measurements(x_input, sides, cutoff; parallel = false)

    # The dual type of measurement does not propagate well on
    # all operations of CellListMap, and also comes with a significant
    # cost, because they are not isbits types. Therefore, it is better
    # to bypass the internal computations by closing over the
    # input array and computing the result directly from the original
    # data. Here, `x_input` is the original array of coordinates
    # with uncertainties.
    function sum_sqr_pair(i, j, s, x_input, box)
        xi = x_input[i]
        xj = CellListMap.wrap_relative_to(x_input[j], xi, box.input_unit_cell.matrix)
        s += norm_sqr(xi - xj)
        return s
    end

    # Copy input vector of coordinates with uncertainties to a vector
    # of coordinates with values only
    x = [getproperty.(v, :val) for v in x_input]

    # Build cell lists
    box = CellListMap.Box(sides, cutoff)
    cl = CellListMap.CellList(x, box)

    # The initial result is initialized with the appropriate type.
    s = measurement(0.0, 0.0)

    # And instead of using the `x` and `y` coordinates provided by the
    # interface, we close over the `x_input` and `box` for the calculations.
    s = CellListMap._pairwise!(
        (pair, s) -> sum_sqr_pair(pair.i, pair.j, s, x_input, box),
        s, box, cl, parallel = parallel
    )
    return s
end

#
# Using generic types
#
function generic_types(iprint = true; parallel = false)

    #
    # Compare analtical and finite-difference gradient. Note that we convert the `sides`
    # and `cutoff` variables to the type of variable in `x`, to allow the propagation of the
    # dual numbers required to ForwardDiff and Unitful
    #
    x = rand(3, 1000)
    cutoff = 0.1
    sides = [1.0, 1.0, 1.0]

    function sumsq_generic(x)
        cutoff_generic = eltype(x).(cutoff) # cutoff is closed over
        sides_generic = eltype(x).(sides) # sides is closed over
        return sumsq(x, sides_generic, cutoff_generic)
    end
    check_grad = sumsq_grad(x, sides, cutoff) ≈ ForwardDiff.gradient(sumsq_generic, x)
    iprint && println("Gradients are correct: ", check_grad)

    #
    # Now let us propagate units
    #
    x = rand(3, 1000)u"nm"
    cutoff = 0.1u"nm"
    sides = [1.0, 1.0, 1.0]u"nm"
    runit = sumsq(x, sides, cutoff)
    iprint && println("Result with units: ", runit)

    #
    # Propagating uncertainties with Measurements
    #
    x = [SVector{3}(measurement(rand(), 0.01 * rand()) for i in 1:3) for j in 1:1000]
    cutoff = 0.1
    sides = [1.0, 1.0, 1.0]
    rmeasurement = sumsq_measurements(x, sides, cutoff, parallel = parallel)
    iprint && println("Result with uncertainty: ", rmeasurement)

    return check_grad, unit(runit), typeof(rmeasurement)
end

#
# Using generic types
#
function generic_types_triclinic(iprint = true; parallel = false)

    #
    # Compare analtical and finite-difference gradient. Note that we convert the `sides`
    # and `cutoff` variables to the type of variable in `x`, to allow the propagation of the
    # dual numbers required to ForwardDiff and Unitful
    #
    x = rand(3, 1000)
    cutoff = 0.1
    sides = [
        1.0 0.0 0.0
        0.0 1.0 0.0
        0.0 0.0 1.0
    ]

    function sumsq_generic(x)
        cutoff_generic = eltype(x).(cutoff) # cutoff is closed over
        sides_generic = eltype(x).(sides) # sides is closed over
        return sumsq(x, sides_generic, cutoff_generic)
    end
    check_grad = sumsq_grad(x, sides, cutoff) ≈ ForwardDiff.gradient(sumsq_generic, x)
    iprint && println("Gradients are correct: ", check_grad)

    #
    # Now let us propagate units
    #
    x = rand(3, 1000)u"nm"
    cutoff = 0.1u"nm"
    sides = [
        1.0 0.0 0.0
        0.0 1.0 0.0
        0.0 0.0 1.0
    ]u"nm"
    runit = sumsq(x, sides, cutoff)
    iprint && println("Result with units: ", runit)

    #
    # Propagating uncertainties with Measurements
    #
    x = [SVector{3}(measurement(rand(), 0.01 * rand()) for i in 1:3) for j in 1:1000]
    cutoff = 0.1
    sides = [
        1.0 0.0 0.0
        0.0 1.0 0.0
        0.0 0.0 1.0
    ]
    rmeasurement = sumsq_measurements(x, sides, cutoff, parallel = parallel)
    iprint && println("Result with uncertainty: ", rmeasurement)

    return check_grad, unit(runit), typeof(rmeasurement)
end
