# Foundation Drawing
## Deriving the Standard Model from Pure Geometry

**Author:** Sally Peck  
**Version:** 1.0, **License:** MIT  

[DOI:10.5281/zenodo.17918142](https://doi.org/10.5281/zenodo.17918142)

This code accompanies the paper: 
[DOI:10.5281/zenodo.17921745](https://doi.org/10.5281/zenodo.17921745)

## Brief Overview

**Foundation Drawing** is a scientific software repository that accompanies the research paper *"The Fine Structure Constant from Planck-Scale Geometry."*

It implements a novel framework based on geometry that derives the fundamental parameters of the Standard Model (couplings, masses, mixing angles) and uses absolutely **zero free parameters**.

It is called **FoundationDrawing**, a play on words referencing the foundational geometry of the universe and simple structures that look as if they were merely "hand-drawn."

## Key Achievements
This code verifies 32 physical quantities (covering the complete set of **26 Standard Model fundamental parameters, including 6 derived physical quantities such as - specific masses of the Proton and W boson)**. **27 of the 32 predictions** match experimental consensus to within 0.01% error.

Prior calculations for $\alpha$ relied on running differential equations from high energy to low energy, AdS/CFT methods, or unexplained numbers. Almost no one has a closed-form geometric derivation for $\sin^2\theta_W$ or the Cabibbo angle that matches experimental values.

Unlike standard numerical fits or ad hoc justifications, this code computes physical constants exclusively from the combinatorial geometry of a single primitive shape: **an equilateral triangle sharing vertices with a semicircular arc.**

| Parameter | Symbol | Geometric Prediction | Experimental (CODATA/PDG) | Precision |
| :--- | :--- | :--- | :--- | :--- |
| **Fine Structure Constant** | $\alpha^{-1}$ | **137.0359991765...** | 137.035999177 | **11 sig figs** |
| **Cabibbo Angle** | $\sin \theta_C$ | **0.2253** | 0.2250 | **0.01%** |
| **Weinberg Angle** | $\sin^2 \theta_W$ | **0.23124** | 0.23122 | **0.007%** |
| **Koide's Ratio** | $Q$ | **0.666666...** | 0.666661 | **0.0008%** |
| **Higgs-Z Ratio** | $m_H/m_Z$ | **1.3737** | 1.3737 | **0.0001%** |
| **Muon Mass** | $m_\mu$ | **105.66 MeV** | 105.66 MeV | **0.0002%** |

## The Geometric Primitive

The entire theory rests on counting the features of a single geometric object $\mathcal{P}$:

```text
            C
           / \           • 3 Vertices (A, B, C)
          /   \          • 2 Shared Vertices (A, B)
         /     \         • 1 Arc (Phase π)
        A───────B        • Tiling: 6 (Flat) vs 5 (Curved)
         \     /

          \   /
           ) (
            v
```

### Repository Overview

| File |Purpose |
|:--- | :--- | 
| ```constants.py``` | Geometric integers (2, 3, 5, 11) and mathematical constants (π, √3​). The core geometric configuration.|
|```foundationdrawing.py``` | Contains the FoundationDrawing class with the exact closed-form algebraic formulas for all 26 parameters and 6 derived quantities. |
| ```standardmodel.py``` | The main execution script. Loads the EXPERIMENTAL data and runs the verification suite. |
| ```checks.py``` | Utility functions for computing error percentage, significant figures, and σ deviation. |


### Prerequisites
* Python 3.8+
* **No external dependencies** (Standard Library only - No NumPy, no SciPy, just pure Python!)

### How to run
```
python3 standardmodel.py
```

### Sample Output:
```
--------------------------------------------------------------------------------
Foundation Drawing - Deriving the Standard Model from Pure Geometry.
(c) 2025 Sally Peck
--------------------------------------------------------------------------------
No free parameters! No magic numbers!
All exponents as named constants (DIMENSION_2D, DIMENSION_3D)
All integers as primitive shape constants (VERTICES_SHARED, VERTICES_TRIANGLE)
--------------------------------------------------------------------------------
...
Predicted (50 digits): 137.03599917650125533338797360700425817294654643249
CODATA 2022:           137.035999177 +/- 2.1E-8
Difference:            4.9874466661202639299574182705345356751E-10
Statistical:           0.02sigma (Perfect Match)
Significant figures:   11.4,
-------------------------------------------------------------------------------
...
Excellent (<0.01% error): 27/32
```
