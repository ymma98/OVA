# OVA
Orthogonal curvilinear coordinate Vector Analysis

# Motivation

Two general vector analysis packages GVA [1] and the SymFields [2] have been developed before using the Mathematica and Python respectively. 

Usually we want to normalize the basis vectors in orthogal curvilinear coordinates (for example the vector analysis formulas in cylindrical and spherical coordinates as what we usually see in many textbook appendices) and keep the basis vectors unnormalized in non-orthogonal coordinates (It seems meaningless to normalize the basis vectors in non-orthogonal coordinates because the co-variant and contra-variant vectors are still different with each other after the normalization, while the normalization process will make the formula derivation much more tedious). On the other hand, expanding the vector expressions including the 2-rank tensor in the component form is necessary in many senarios, for example the calculations of $\nabla \cdot \mathbf{T}$, $\vec{u} \cdot \mathbf{T}$ and $\vec{u} \cdot \nabla \vec{u}$ in a given coordinate system, and it remains to be enhanced in current packages.

In the OVA, we only perform the symbolic vector analysis for the orthogonal coordinate system. The basis vectors are normalized as we expected, and the package supports some 2-rank tensor calculations that are commonly used (for me). For the vector analysis in the non-normalized coordinate system, I recommand to use the general purpose vector analysis packages.

The examples can be found in the benchmark/ folder.

[1] Qin H *et al.* Symbolic vector analysis in plasma physics[J]. Computer physics communications, 1999, 116(1): 107-120.

[2] Chu N. SymFields: An Open Source Symbolic Fields Analysis Tool for General Curvilinear Coordinates in Python[C]//Proceedings of the Future Technologies Conference. Springer, Cham, 2021: 685-696.

# Notes

* When using the OVA package, the vector is represented as a `list` of `sympy` expressions. However, when performing the `+, -, *, \` operations on two vectors or a vector to a scalar, we have to either convert vector in type of `list` to `sympy.Array`, or manipulate each component by hand. For example, the basic operations on $\vec{u}$ is not as what we expected due to the nature of the `list` type.
```python
u = [expr1, expr2, expr3]  # u is a vector
u + u
u * 3

output:
[expr1, expr2, expr3, expr1, expr2, expr3]
[expr1, expr2, expr3, expr1, expr2, expr3, expr1, expr2, expr3]
```

Two possible approaches:
```python
u = [expr1, expr2, expr3]  # u is a vector
u = sympy.Array(u)
u+u

output:
[expr1 * 2, expr2 * 2, expr3 * 2] 
```
or
```python
u = [expr1, expr2, expr3]  # u is a vector
u_plus_u = [k+k for k in u]
u_plus_u

output：
[expr1 * 2, expr2 * 2, expr3 * 2] 
```

* Why the vectors are not used in the type of `sympy.Array` when doing the curl or other vector anaylses.


`Array` in `sympy` is an abbreviation for `ImmutableDenseNDimArray`! The `Array` in `sympy` is immutable, which prevents any modifications in the component level.