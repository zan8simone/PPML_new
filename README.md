# PPML_new
Periodically Patterned MultiLayer - Electromagnetic properties of patterned multilayers based on RCWA (Rigorous Coupled Wave Analysis, a.k.a. FMM - Fourier Modal Method)

This repository is the prosecution of the older [PPML---Periodically-Patterned-Multi-Layer](https://github.com/zan8simone/PPML---Periodically-Patterned-Multi-Layer). Due to major code style and numerical method revision, I created a new repository, but the versioning label starts from v4.0, to follow-up after the latest version in PPML---Periodically-Patterned-Multi-Layer.

## Requirements

Tested with MATLAB version R2023b.

## Electromagnetic properties of patterned multilayers based on RCWA

Rigorous coupled wave analysis (RCWA) based on the scattering matrix (SM) algorithm is one of the most powerful tools for the electromagnetic simulation of patterned multilayer structures.

PPML is a project which implements the SM-RCWA, based on the formalisms of [a,b]. Three groups of functions are available: one is for 1-d patterns under TM polarization, one is for 1-d patterns under TE polarization, another is for 2-d non-orthogonal unit cells with anisotropic dielectric. In all cases, conducting interfaces (e.g. 2d materials) between the layers in the stack are handled natively as a zero-thickness boundary condition.

### 1d_TE and 1d_TM PPML

For 1-d patterns under TE or TM polarization, currently available are functions for the calculation of
- reflectance and transmittance for all diffraction orders
- layer-by-layer absorptance
- scattering matrix of the zero-order channel (i.e., complex reflection and transmission coefficients)

Additional features for 1-d patterns under TM polarization:
- calculation of E and S fields inside the structure
- implementation of Fourier factorization rules [c,d], that make the code fully suitable for the simulation of metal components (plasmonic gratings)
- uniaxial materials (permittivity tensor diagonal in x-y-z reference system)

### 2d_basic_general PPML

For 2-d patterns, currently available are functions for the calculation of
- intensity reflectance and transmittance for all diffraction orders (i.e., diffraction efficiencies)
- layer-by-layer absorptance
- Jones matrix for all pairs of diffraction channels
- calculation of E and H fields inside the structure

The 2-d solver allows to work with arbitrarily shaped unit cells with arbitrary permittivity distribution (tensor with non-zero xx, xy, yx, and zz components). The drawback of such generality is that the solver is based on ordinary (fast) Fourier transform, without the use of factorization rules or local coordinate transformations. Hence, it may be inadequate when materials with strong imaginary permittivity are investigated. 


### Relation with previous versions of PPML

The 1d_TE and 1d_TM codes of PPML_new are essentially the same as what found in older versions ([github](https://github.com/zan8simone/PPML---Periodically-Patterned-Multi-Layer) and [Matlab File Exchange](https://it.mathworks.com/matlabcentral/fileexchange/55401-ppml-periodically-patterned-multi-layer)). The 2d_basic_general version has major differences with respect to the 2d_ functions available in older PPML versions. In detail:

| | Older PPML | PPML v4.0|
|-|------------|---------|
| Functions: |2d_rect, 2d_Lshape | 2d_basic_general |
| Unit cell: |rectangular |  arbitrary |
| Pattern: |rectangle or L-shape inclusion | arbitrary|
| Permittivity: | scalar-valued | tensor-valued |
| RCWA core: |Fourier factorization rules (FFR) and analytic Fourier transforms | No FFR, FTs based on FFT |

If you wish to study a rectangular unit cell with rectangular or L-shaped inclusions, it is advisable to use PPML v3.0.  

### Reference

The present code is distributed for free, but we kindly ask you to cite its source and, if applies, the publications below.
Several publications are based on PPML (see below). For some of them, you can find the corresponding tutorial in the software package.

## Publications based on various versions of PPML (not exhaustive list)

1. S. Zanotto, ... A. Pitanti, "Photonic bands, superchirality, and inverse design of a chiral minimal metasurface", Nanophotonics (2019), DOI: 10.1515/nanoph-2019-0321
2. S. Zanotto, ... A. Pitanti, "Optomechanics of Chiral Dielectric Metasurfaces", Advanced Optical Materials (2020), DOI: 10.1002/adom.201901507
3. S. Zanotto, ... D. S. Wiersma, "Multichannel remote polarization control enabled by nanostructured Liquid Crystalline Networks", Applied Physics Letters (2019), DOI 10.1063/1.5096648
4. S. Zanotto, G. C. La Rocca, and A. Tredicucci, “Understanding and overcoming fundamental limits of asymmetric light-light switches”, Optics Express 26, 3, 3618 (2018).
5. S. Zanotto, ..., A. Melloni, "Metasurface reconfiguration through lithium ion intercalation in a transition metal oxide", Advanced Optical Materials 2017, 5, 1600732 (2017).
6. S. Zanotto and A. Tredicucci, "Universal lineshapes at the crossover between weak and strong critical coupling in Fano-resonant coupled oscillators", Scientific Reports 6, 24592 (2016).
7. L. Baldacci, ..., A. Tredicucci, “Interferometric control of absorption in thin plasmonic metamaterials: general two port theory and broadband operation”, Optics Express 23, 9202 (2015).
8. S. Zanotto, ... A. Tredicucci, “Perfect energy feeding into strongly coupled systems and interferometric control of polariton absorption”, Nature Physics 10, 830 (2014).
9. J.-M. Manceau, ..., R. Colombelli, “Mid-infrared intersubband polaritons in dispersive metal-insulator-metal resonators”, Appl. Phys. Lett. 105, 081105 (2014).
10. J.-M. Manceau, ..., R. Colombelli, “Optical critical coupling into highly confining metal-insulator-metal resonators”, Appl. Phys. Lett. 103, 091110 (2013).
11. S. Zanotto, ..., A. Tredicucci, “Ultrafast optical bleaching of intersubband cavity polaritons”, Phys. Rev. B 86, 201302(R) (2012).
12. S. Zanotto, ... A. Tredicucci, “Analysis of line shapes and strong coupling with intersubband transitions in one-dimensional metallodielectric photonic crystal slabs”, Phys. Rev. B 85, 035307 (2012).
13. R. Degl'Innocenti, ... A. Tredicucci, “One-dimensional surface-plasmon gratings for the excitation of intersubband polaritons in suspended membranes”, Solid State Comm. 151, 1725-1727 (2011).
14. S. Zanotto, ... A. Tredicucci, “Intersubband polaritons in a one-dimensional surface plasmon photonic crystal”, Appl. Phys. Lett. 97, 231123 (2010).

## References for the method

a.     D. M. Whittaker & I. S. Culshaw, "Scattering-matrix treatment of patterned multilayer photonic structures",
Phys. Rev. B 60, 2610 (1999).

b.    M. Liscidini, D. Gerace, L. C. Andreani & J. E. Sipe, "Scattering-matrix analysis of periodically patterned multilayers with asymmetric unit cells and birefringent media", Phys. Rev. B 77, 035324 (2008).

c.     L. Li. "Use of Fourier series in the analysis of discontinuous periodic structures". J. Opt. Soc. Am. A 13, 1870 (1996).

d.    Lalanne, Philippe, and G. Michael Morris, "Highly improved convergence of the coupled-wave method for TM polarization", JOSA A 13, 779 (1996).
