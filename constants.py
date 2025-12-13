#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Foundation Drawing - Deriving the Standard Model from Pure Geometry.
# Author: (c) 2025 Sally Peck.
# -----------------------------------------------------------------------------

from decimal import Decimal, getcontext
import math

getcontext().prec = 50  # 50 decimal places so 10+ dp digits compute correctly

# -----------------------------------------------------------------------------
# Basic constants
# -----------------------------------------------------------------------------
ZERO = 0                      # Additive identity
DECIMAL_ZERO = 0.0            # Float zero for returns
UNITY = 1                     # Multiplicative identity / one complete shape
DECIMAL_UNITY = 1.0           # Float unity for division

PI = math.pi
PI_DECIMAL = Decimal('3.14159265358979323846264338327950288419716939937510')

SQRT3 = math.sqrt(3)
SQRT3_DECIMAL = Decimal(3).sqrt()

# -----------------------------------------------------------------------------
# Level 1: Basic Counting (directly from primitives)
# -----------------------------------------------------------------------------
VERTICES_SHARED = 2           # Half-circle endpoints A, B shared with triangle
VERTICES_TRIANGLE = 3         # Triangle has 3 vertices: A, B, C

# Decimal versions for high-precision calculations
VERTICES_SHARED_D = Decimal(VERTICES_SHARED)
VERTICES_TRIANGLE_D = Decimal(VERTICES_TRIANGLE)

# -----------------------------------------------------------------------------
# Level 2: Dimensional Factors (for exponents)
# -----------------------------------------------------------------------------
DIMENSION_1D = 1              # 1D phase space (line/arc)
DIMENSION_2D = 2              # 2D phase space (plane) - USE FOR **2
DIMENSION_3D = 3              # 3D phase space (volume) - USE FOR **3
DIMENSION_FACTOR = 4          # 2² = phase space dimension factor for EM

# Decimal versions
DIMENSION_FACTOR_D = Decimal(DIMENSION_FACTOR)

# -----------------------------------------------------------------------------
# Level 3: Triangles per Vertex (Tiling)
# -----------------------------------------------------------------------------
TRIANGLE_CURVED = 5           # In curved space: 5 triangles/vertex (300° total)
TRIANGLE_FLAT = 6             # In flat space: 6 triangles/vertex (360° total)
TRIANGLES_PER_EDGE = 2        # Two triangles share each edge

# Decimal versions
TRIANGLE_CURVED_D = Decimal(TRIANGLE_CURVED)
FLAT_TRIANGLES_D = Decimal(TRIANGLE_FLAT)
TRIANGLES_PER_EDGE_D = Decimal(TRIANGLES_PER_EDGE)

# -----------------------------------------------------------------------------
# Level 4: First-Order Combinations
# -----------------------------------------------------------------------------
VERTICES_PLUS_CURVED = VERTICES_TRIANGLE + TRIANGLE_CURVED    # 3 + 5 = 8
FLAT_PLUS_CURVED = TRIANGLE_FLAT + TRIANGLE_CURVED            # 6 + 5 = 11
CURVED_PLUS_SHARED = TRIANGLE_CURVED + VERTICES_SHARED        # 5 + 2 = 7

# Decimal versions
VERTICES_PLUS_CURVED_D = Decimal(VERTICES_PLUS_CURVED)
FLAT_PLUS_CURVED_D = Decimal(FLAT_PLUS_CURVED)

# -----------------------------------------------------------------------------
# Level 5: Second-Order Combinations
# -----------------------------------------------------------------------------
FLAT_TIMES_VERTICES = TRIANGLE_FLAT * VERTICES_TRIANGLE       # 6 × 3 = 18
EDGE_TIMES_FLAT = TRIANGLES_PER_EDGE * TRIANGLE_FLAT          # 2 × 6 = 12
FLAT_SQUARED = TRIANGLE_FLAT * TRIANGLE_FLAT                  # 6² = 36
VERTICES_CUBED = VERTICES_TRIANGLE ** VERTICES_TRIANGLE       # 3³ = 27
VERTICES_FOURTH = VERTICES_CUBED * VERTICES_TRIANGLE          # 3⁴ = 81

# Combinations for specific formulas
FLAT_PLUS_CURVED_PLUS_FLAT = FLAT_PLUS_CURVED + TRIANGLE_FLAT            # 11 + 6 = 17
FLAT_PLUS_VERTICES = TRIANGLE_FLAT + VERTICES_TRIANGLE                   # 6 + 3 = 9
EDGE_TIMES_CURVED_PLUS_SHARED = TRIANGLES_PER_EDGE * CURVED_PLUS_SHARED  # 2 × 7 = 14
ALL_CONFIGS_PLUS_VERT_CURVED = FLAT_PLUS_CURVED + VERTICES_PLUS_CURVED   # 11 + 8 = 19
FLAT_SQUARED_PLUS_CURVED = FLAT_SQUARED + TRIANGLE_CURVED                # 36 + 5 = 41

# For mass ratios
MASS_OFFSET_7 = CURVED_PLUS_SHARED                # 7 = 5 + 2
MASS_OFFSET_9 = FLAT_PLUS_VERTICES                # 9 = 6 + 3
MASS_OFFSET_12 = EDGE_TIMES_FLAT                  # 12 = 2 × 6
MASS_OFFSET_14 = EDGE_TIMES_CURVED_PLUS_SHARED    # 14 = 2 × 7

# For top/bottom ratio
TOP_BOTTOM_BASE = FLAT_SQUARED_PLUS_CURVED        # 41 = 36 + 5

# Special combinations
ICOSAHEDRON_FACES = 20                            # Platonic solid: 20 faces
VERT_CURVED_TIMES_ICOSA = VERTICES_PLUS_CURVED * ICOSAHEDRON_FACES    # 8 × 20 = 160
VERTICES_FOURTH_PLUS_VERT_CURVED = VERTICES_FOURTH + VERTICES_PLUS_CURVED  # 81 + 8 = 89

# For CKM matrix
CKM_OFFSET_4 = DIMENSION_FACTOR                   # 4 = 2²
CKM_OFFSET_160 = VERT_CURVED_TIMES_ICOSA          # 160 = 8 × 20

# For PMNS neutrino mixing
PMNS_OFFSET_27 = VERTICES_CUBED                   # 27 = 3³

# For neutrino mass hierarchy
NEUTRINO_89 = VERTICES_FOURTH_PLUS_VERT_CURVED    # 89 = 81 + 8

# For Higgs sector
HIGGS_DENOMINATOR = ALL_CONFIGS_PLUS_VERT_CURVED  # 19 = 11 + 8

# For α⁻¹ formula
ALPHA_SHARING_COEFF = VERTICES_SHARED             # 2 in numerator
ALPHA_SHARING_COEFF_D = VERTICES_SHARED_D
ALPHA_VERTEX_COEFF = VERTICES_TRIANGLE            # 3 in numerator
ALPHA_VERTEX_COEFF_D = VERTICES_TRIANGLE_D
ALPHA_EDGE_NORM = TRIANGLES_PER_EDGE              # 2 in denominator
ALPHA_EDGE_NORM_D = TRIANGLES_PER_EDGE_D
ALPHA_CURVED_COEFF = TRIANGLE_CURVED              # 5 in numerator
ALPHA_CURVED_COEFF_D = TRIANGLE_CURVED_D
ALPHA_CURVED_VERTEX = VERTICES_PLUS_CURVED        # 8 in numerator
ALPHA_CURVED_VERTEX_D = VERTICES_PLUS_CURVED_D
ALPHA_CONFIG_NORM = FLAT_PLUS_CURVED              # 11 in denominator
ALPHA_CONFIG_NORM_D = FLAT_PLUS_CURVED_D

# For αs formula
STRONG_OFFSET = FLAT_PLUS_CURVED_PLUS_FLAT        # 17 = 11 + 6

# For sin²θ_W formula
WEAK_ENDPOINT_COEFF = VERTICES_SHARED             # 2π coefficient
WEAK_DENOMINATOR_BASE = FLAT_TIMES_VERTICES       # 18 = 6 × 3

# For percentage calculations
PERCENT_FACTOR = 100
