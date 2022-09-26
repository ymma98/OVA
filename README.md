# OVA
Orthogonal curvilinear coordinate Vector Analysis

# Note

* When using the OVA package, the vector is represented as a `list` of `sympy` expressions. However, when performing the `+, -, *, \` operations on two vectors or a vector to a scalar, we have to either convert vector in type of `list` to `sympy.Array`, or manipulate each component by hand. For example, the basic operations on $\vec{u}$ is not as what we expected due to the nature of the `list` type.
```python
u = [expr1, expr2, expr3]  # u is a vector
u + u
u * 3

output:
[expr1, expr2, expr3, expr1, expr2, expr3]
[expr1, expr2, expr3, expr1, expr2, expr3, expr1, expr2, expr3]
```

Two possible approach:
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