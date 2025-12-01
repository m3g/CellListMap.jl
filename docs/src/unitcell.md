# Unitcell requirements

The unitcell provided by the user must satisfy certain conditionsfor the cell list algorithm to work correctly.
If the conditions are not satisfied, the computation will throw an early error. These conditions guarantee
that the box is not too small for the given required cutoff, considering the perodic conditions. In the 
simple case of orthorhombic boxes, the conditions imply that each side must be greater than twice the cutoff.
These conditions are usually met except for very small systems, which are not within the current scope of this package.

## In 3D

Given an unitcell matrix of the form $\left[\vec{a}~\vec{b}~\vec{c}\right]$, where $\vec{a}$, $\vec{b}$, and 
$\vec{c}$ are the lattice vectors, the following conditions must be met:

```math
\vec{a} \cdot \frac{\vec{b} \times \vec{c}}{|\vec{b} \times \vec{c}|} \gt 2r_{cut}
```
```math
\vec{b} \cdot \frac{\vec{c} \times \vec{a}}{|\vec{c} \times \vec{a}|} \gt 2r_{cut}
```
```math
\vec{c} \cdot \frac{\vec{a} \times \vec{b}}{|\vec{a} \times \vec{b}|} \gt 2r_{cut}
```

where $r_{cut}$ is the cutoff. 

## In 2D

Given an unitcell matrix of the form $\left[\vec{a}~\vec{b}\right]$, where $\vec{a}$ and $\vec{b}$ 
are the lattice vectors, the following conditions must be met:

```math
\sqrt{|\vec{a}|^2 - \left(\vec{a}\cdot \frac{\vec{b}}{|\vec{b}|}\right)^2} \gt 2r_{cut} 
```
```math
\sqrt{|\vec{b}|^2 - \left(\vec{b}\cdot \frac{\vec{a}}{|\vec{a}|}\right)^2} \gt 2r_{cut} 
```
where $r_{cut}$ is the cutoff. 