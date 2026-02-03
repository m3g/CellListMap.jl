#
# This interface is needed to generate random particle coordinates in ComplexMixtures.jl
#
"""
    get_computing_box(sys::AbstractParticleSystem)

Retrieves the computing box of the system. The computing box is large enough to
contain all coordinates of the particles, plus the cutoff.

"""
get_computing_box(sys::AbstractParticleSystem) = sys._box.computing_box