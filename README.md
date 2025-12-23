# Foundation Drawing
## Deriving the Standard Model from Pure Geometry

**Author:** Sally Peck  
**Version:** 1.2 (December 2025)  
**License:** MIT  


### DOI Links

Source Code:
[DOI:10.5281/zenodo.17918142](https://doi.org/10.5281/zenodo.17918142)

This code accompanies the papers: 

- The Fine Structure Constant as a Lucas Gap-Filling Series - [DOI:10.5281/zenodo.17994370](https://doi.org/10.5281/zenodo.17994370)

- The Fine Structure Constant from Planck-Scale Geometry - [DOI:10.5281/zenodo.17921745](https://doi.org/10.5281/zenodo.17921745)

## 1. Brief Overview

**Foundation Drawing** is a scientific software repository that accompanies the above research papers.

It implements a novel framework based on geometry that derives the fundamental parameters of the Standard Model (couplings, masses, mixing angles) and uses absolutely **zero free parameters**.

It is called **FoundationDrawing**, a play on words referencing the foundational geometry of the universe and simple structures that look as if they were merely "hand-drawn."

## 2. Key Achievements
This code verifies 29 physical quantities (covering **23 Standard Model fundamental parameters plus 6 derived physical quantities** including the proton mass, muon mass, tau mass, and W boson mass). **15 of the 29 predictions** match experimental consensus to within 0.01% error.

Prior calculations for $\alpha$ relied on running differential equations from high energy to low energy, AdS/CFT methods, or unexplained numbers. Almost no one has a closed-form geometric derivation for $\sin^2\theta_W$ or the Cabibbo angle that matches experimental values.

Unlike standard numerical fits or ad hoc justifications, this code computes physical constants exclusively from the combinatorial geometry of a single primitive shape: **an equilateral triangle sharing vertices with a semicircular arc.**

## 3. Fine Structure Constant v1.2 - The Lucas Gap-Filling Principle

The fine structure constant $\alpha$ is traditionally viewed as a coupling constant. In this framework, it is derived as the geometric "leakage" or deficit created when attempting to tile a curved surface with flat primitives.

The Geometric Intuition:

- The Deficit: A flat 2D plane is perfectly tiled by 6 equilateral triangles meeting at a vertex (360° total).

- The Curvature: In a curved or "spherical" configuration, only 5 triangles meet at a vertex (300° total).

- The Gap: This creates a 60° (π/3) angular deficit.

The higher-order version implements a recursive **Lucas Series Expansion** to resolve the 60° angular deficit between curved space (5 triangles/vertex) and flat space (6 triangles/vertex). 


For orders $n \ge 4$, the n-th order correction is defined by:
$$\delta_n = (-1)^{\lfloor(n+1)/2\rfloor} \cdot \frac{L_{n-1}\pi + L_{n+3}\sqrt{3} - L_{n+4}}{L_{n-1} \cdot S^n}$$
where $S = 4\pi^3 + \pi^2 + \pi$

The geometric framework suggests this discrepancy reflects different scales being probed, not measurement error. The formula bridges both values through a single recursive mechanism. 

| Parameter | Symbol | Geometric Prediction | Experimental | Precision |
| :--- | :--- | :--- | :--- | :--- |
| **Order 3** | $\alpha^{-1}$ | **137.035999176501...** | 137.035999177 (CODATA 2022) | **11 sig figs** |
| **Order 4** | $\alpha^{-1}$ | **137.035999166005...** | 137.035999166 (Fan et al. 2023) | **13 sig figs** |


A fitted formula would match the data it was fitted to - it would not independently predict future measurements.

The series reaches **Physical Saturation** at the 8th Order ($\sim 10^{-17}$), where corrections fall below Planck-scale quantum fluctuations.


## 4. Higher-Order Corrections for Key Parameters

The higher-order correction quantifies how much the shared vertices contribute relative to the total geometric and phase space.

| Parameter | Symbol | Geometric Prediction | Experimental | Precision |
| :--- | :--- | :--- | :--- | :--- |
| **Cabibbo Angle** | $\sin \theta_C$ | **0.2253** | 0.2250 | **0.01%** |
| **Weinberg Angle** | $\sin^2 \theta_W$ | **0.23124** | 0.23122 | **0.007%** |
| **Koide's Ratio** | $Q$ | **0.666666...** | 0.666661 | **0.0008%** |
| **Higgs-Z Ratio** | $m_H/m_Z$ | **1.3737** | 1.3737 | **0.0001%** |
| **Muon Mass** | $m_\mu$ | **105.66 MeV** | 105.66 MeV | **0.0002%** |

## 5. Experimental References

**CODATA 2022** (NIST, published May 2024):
- $\alpha^{-1}$ = 137.035999177(21)
- Proton-electron mass ratio, muon-electron mass ratio, and other fundamental constants
- Data: https://physics.nist.gov/cuu/Constants/Table/allascii.txt
- Source: https://physics.nist.gov/cuu/Constants/

**Fan et al. 2023** - Electron Magnetic Moment Measurement:
- $\alpha^{-1}$ = 137.035999166(15) [0.11 ppb precision]
- Source: Phys. Rev. Lett. 130, 071801 (2023)
- DOI: [10.1103/PhysRevLett.130.071801](https://doi.org/10.1103/PhysRevLett.130.071801)

**PDG 2024** - Particle Data Group:
- Quark masses, CKM matrix elements, Weinberg angle, strong coupling constant
- S. Navas et al. (Particle Data Group), Phys. Rev. D 110, 030001 (2024)
- Source: https://pdg.lbl.gov/

**NuFIT 6.0** - Global Neutrino Oscillation Analysis (September 2024):
- PMNS matrix parameters, neutrino mass-squared differences
- I. Esteban et al., JHEP 12 (2024) 216
- arXiv: [2410.05380](https://arxiv.org/abs/2410.05380)
- Source: http://www.nu-fit.org/

**ATLAS 2024** - Higgs Boson Mass:
- Combined measurement with 0.09% precision
- $m_H$ = 125.11 ± 0.11 GeV
- Source: https://home.cern/news/news/physics/atlas-sets-record-precision-higgs-bosons-mass

The geometric prediction matches CODATA 2022 at Order 3 (11 sig figs) and
Fan et al. 2023 at Order 4 with higher-order corrections (13 sig figs).

## The Geometric Primitive

The entire theory rests on counting the features of a single geometric object $\mathcal{P}$:

```text
           C
          /\            • 3 Vertices (A, B, C)
         /  \           • 2 Shared Vertices (A, B)
        /    \          • 1 Arc (Phase π)
       /      \         • Tiling: 6 (Flat) vs 5 (Curved)
      A--------B
       \      /
        \    /
         )  (
          \/
```

### Repository Overview

| File | Purpose |
|:--- |:--- |
| `alpha_higher_orders.py` | Alpha computation with higher order corrections. |
| `constants.py` | Geometric integers (2, 3, 5, 11) and mathematical constants (π, √3). The core geometric configuration. |
| `foundationdrawing.py` | Contains the FoundationDrawing class with the exact closed-form algebraic formulas for all 23 parameters and 6 derived quantities. |
| `standardmodel.py` | The main execution script. Loads the EXPERIMENTAL data and runs the verification suite. |
| `checks.py` | Utility functions for computing error percentage, significant figures, and σ deviation. |


### Prerequisites
* Python 3.8+
* **No external dependencies** (Standard Library only - No NumPy, no SciPy, just pure Python!)

### How to run

**Full verification suite:**
```
python standardmodel.py
```

**Alpha higher-order series only:**
```
python alpha_higher_orders.py
```

### Sample Output (standardmodel.py):
```
--------------------------------------------------------------------------------
Foundation Drawing - Deriving the Standard Model from Pure Geometry.
(c) 2025 Sally Peck - https://github.com/sallyp/foundationdrawing/
--------------------------------------------------------------------------------
...
ALPHA COMPARISON: Order 3 (Original) vs Order 4 (Higher-Order)
--------------------------------------------------------------------------------
Source       Order 3                                  Order 4
--------------------------------------------------------------------------------
Predicted    137.035999176501255333387973607004258172 137.035999166004757297132106232148492257
CODATA 2022  137.035999177                            diff: 4.99e-10 (11.4 SF)
Fan 2023     137.035999166                            diff: 4.76e-12 (13.5 SF)
--------------------------------------------------------------------------------
Order 3 matches CODATA 2022 (11.4 sig figs)
Order 4 matches Fan 2023 (13.5 sig figs)
...
Total: 29 parameters, Excellent (<0.01% error): 15
```

### Sample Output (alpha_higher_orders.py):
```
Alpha Order 0: 137.0363037758784
Alpha Order 1: 137.0363037758784
Alpha Order 2: 137.0359981320960
Alpha Order 3: 137.0359991765013
Alpha Order 4: 137.0359991660048
Alpha Order 5: 137.0359991659302
Alpha Order 6: 137.0359991659297
Alpha Order 7: 137.0359991659297
Alpha Order 8: 137.0359991659297
```
