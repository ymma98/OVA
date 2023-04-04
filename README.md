# OVA
Orthogonal curvilinear coordinate Vector Analysis

# Motivation

* In plasma physics research, vector analysis is ubiquitous but often tedious, and the cumbersomeness has even become a trait of the field.

* The automatic symbolic derivation is necessary (at least for me).
* Two symbolic approches towards the vector analysis in general coordinate systems have been developed.
  * GVA [1] (General Vector Analysis), by Prof. Qin in 1997, a Mathematica package. Non-normalized basis is used, which has great generality and theoretical conciseness in any well-defined coordinates. However, we prefer to normalize the basis vectors in orthogonal curvilinear coordinates as is commonly shown in textbooks and cheatsheets.
  * SymFields [2], by Dr. Chu in 2020, a python package. It is also developed for vector analysis in the general coordinates. The package normalize basis vectors in both orthogonal and non-orthogonal coordinates, which complicates the mathematics and brings inconvenience for further developments (e.g. add more supports about 2-rank tensor related calculation).

# OVA

* The OVA can only be used for symbolic vector analysis in the orthogonal coordinate system.  For the vector analysis in the non-orthogonal coordinate system, I recommend to use the two general purpose vector analysis packages [1,2].
* The OVA (Orthogonal curvilinear coordinate Vector Analysis) package has been implemented in both Python and Mathematica versions.
  * The Python version OVA is implemented using the vector analysis formulas assuming the basis normalized by the Lamé coefficients (thus can only be used in orthogonal coordinates).
  * The GVA package, originally developed by Prof. Qin in 1997, has been modified to ensure compatibility with newer versions of Mathematica, beyond version 8.0.
  * The Mathematica version of the OVA package is built upon the modified GVA package, which includes updated vector definitions and the normalization of each component using Lamé coefficients.

# References

[1] Qin H *et al.* Symbolic vector analysis in plasma physics[J]. Computer physics communications, 1999, 116(1): 107-120.

[2] Chu N. SymFields: An Open Source Symbolic Fields Analysis Tool for General Curvilinear Coordinates in Python[C]//Proceedings of the Future Technologies Conference. Springer, Cham, 2021: 685-696.

