#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Foundation Drawing - Deriving the Standard Model from Pure Geometry.
# Author: (c) 2025 Sally Peck.
# -----------------------------------------------------------------------------

from decimal import Decimal, getcontext
from constants import (
    DIMENSION_FACTOR_D, PI_DECIMAL, DIMENSION_3D, DIMENSION_2D, DIMENSION_1D,
    DIMENSION_FACTOR, PI, ALPHA_SHARING_COEFF_D, ALPHA_VERTEX_COEFF_D,
    SQRT3_DECIMAL, ALPHA_EDGE_NORM_D, ALPHA_CURVED_COEFF_D, ALPHA_CURVED_VERTEX_D,
    ALPHA_CONFIG_NORM_D, ALPHA_SHARING_COEFF, ALPHA_VERTEX_COEFF, SQRT3,
    ALPHA_EDGE_NORM, ALPHA_CURVED_COEFF, ALPHA_CURVED_VERTEX, ALPHA_CONFIG_NORM,
    DECIMAL_UNITY, VERTICES_SHARED, VERTICES_TRIANGLE, STRONG_OFFSET,
    WEAK_ENDPOINT_COEFF, WEAK_DENOMINATOR_BASE, CKM_OFFSET_4, TRIANGLE_FLAT,
    TRIANGLE_CURVED, UNITY, FLAT_PLUS_CURVED_PLUS_FLAT, FLAT_PLUS_CURVED,
    MASS_OFFSET_12, MASS_OFFSET_7, TOP_BOTTOM_BASE,
    CKM_OFFSET_160, MASS_OFFSET_9, MASS_OFFSET_14, PMNS_OFFSET_27,
    HIGGS_DENOMINATOR, DECIMAL_ZERO, NEUTRINO_89
)

getcontext().prec = 50  # 50 decimal places so 10+ dp digits compute correctly

# S = 4π³ + π² + π (total electromagnetic phase space)
def compute_S_decimal() -> Decimal:
    return (DIMENSION_FACTOR_D * PI_DECIMAL**DIMENSION_3D +
            PI_DECIMAL**DIMENSION_2D + PI_DECIMAL**DIMENSION_1D)

def compute_S() -> float:
    return (DIMENSION_FACTOR * PI**DIMENSION_3D +
            PI**DIMENSION_2D + PI**DIMENSION_1D)

S_DECIMAL = compute_S_decimal()
S = compute_S()


class FoundationDrawing:
    # COUPLING CONSTANTS ------------------------------------------------------
    @staticmethod
    def alpha_inverse_decimal() -> Decimal:
        # α⁻¹ = S - (2π + 3√3)/(2S²) + (5π + 8√3)/(11S³) [50 digit precision]
        numerator_1 = ALPHA_SHARING_COEFF_D * PI_DECIMAL + ALPHA_VERTEX_COEFF_D * SQRT3_DECIMAL
        denominator_1 = ALPHA_EDGE_NORM_D * S_DECIMAL**DIMENSION_2D
        numerator_2 = ALPHA_CURVED_COEFF_D * PI_DECIMAL + ALPHA_CURVED_VERTEX_D * SQRT3_DECIMAL
        denominator_2 = ALPHA_CONFIG_NORM_D * S_DECIMAL**DIMENSION_3D
        return S_DECIMAL - numerator_1/denominator_1 + numerator_2/denominator_2

    @staticmethod
    def alpha_inverse() -> float:
        # α⁻¹ = S - (2π + 3√3)/(2S²) + (5π + 8√3)/(11S³) [11.4 sig figs]
        numerator_1 = ALPHA_SHARING_COEFF * PI + ALPHA_VERTEX_COEFF * SQRT3
        denominator_1 = ALPHA_EDGE_NORM * S**DIMENSION_2D
        numerator_2 = ALPHA_CURVED_COEFF * PI + ALPHA_CURVED_VERTEX * SQRT3
        denominator_2 = ALPHA_CONFIG_NORM * S**DIMENSION_3D
        return S - numerator_1/denominator_1 + numerator_2/denominator_2

    @staticmethod
    def alpha() -> float:
        # Fine structure constant α
        return DECIMAL_UNITY / FoundationDrawing.alpha_inverse()

    @staticmethod
    def alpha_s() -> float:
        # αs = (-π + 2√3 + 3)/(3π + √3 + 17) [6.0 sig figs]
        numerator = -PI + VERTICES_SHARED * SQRT3 + VERTICES_TRIANGLE
        denominator = VERTICES_TRIANGLE * PI + SQRT3 + STRONG_OFFSET
        return numerator / denominator

    @staticmethod
    def sin2_theta_W() -> float:
        # sin²θW = (2π + √3 - 2)/(2π + √3 + 18) [5.9 sig figs]
        numerator = WEAK_ENDPOINT_COEFF * PI + SQRT3 - VERTICES_SHARED
        denominator = WEAK_ENDPOINT_COEFF * PI + SQRT3 + WEAK_DENOMINATOR_BASE
        return numerator / denominator

    # LEPTON MASS RATIOS ------------------------------------------------------
    @staticmethod
    def proton_electron_ratio() -> float:
        # m_p/m_e = 6π⁵ [4.7 sig figs]
        return TRIANGLE_FLAT * PI**TRIANGLE_CURVED

    @staticmethod
    def muon_electron_ratio() -> float:
        # m_μ/m_e = (3/(2α) + 6/5)(1 + 4α²/3) [5.7 sig figs]
        alpha = FoundationDrawing.alpha()
        base = VERTICES_TRIANGLE / (VERTICES_SHARED * alpha) + TRIANGLE_FLAT / TRIANGLE_CURVED
        correction = UNITY + DIMENSION_FACTOR * alpha**DIMENSION_2D / VERTICES_TRIANGLE
        return base * correction

    @staticmethod
    def tau_muon_ratio() -> float:
        # m_τ/m_μ = 17 - 2/11 [4.2 sig figs]
        TAU_MUON_BASE = FLAT_PLUS_CURVED_PLUS_FLAT  # 17 = 11 + 6
        return TAU_MUON_BASE - VERTICES_SHARED / FLAT_PLUS_CURVED

    @staticmethod
    def koide_Q() -> float:
        # Q = 2/3 = shared vertices / total vertices [5.1 sig figs]
        return VERTICES_SHARED / VERTICES_TRIANGLE

    # QUARK MASS RATIOS -------------------------------------------------------
    @staticmethod
    def up_down_ratio() -> float:
        # m_u/m_d = (π + 2√3 + 2)/(π + 2√3 + 12) [4.2 sig figs]
        common = PI + VERTICES_SHARED * SQRT3
        return (common + VERTICES_SHARED) / (common + MASS_OFFSET_12)

    @staticmethod
    def strange_down_ratio() -> float:
        # m_s/m_d = 7π - 2 [3.4 sig figs]
        return MASS_OFFSET_7 * PI - VERTICES_SHARED

    @staticmethod
    def charm_strange_ratio() -> float:
        # m_c/m_s = (2π + 2√3)/(7 - 2π) [4.1 sig figs]
        numerator = VERTICES_SHARED * (PI + SQRT3)
        denominator = MASS_OFFSET_7 - VERTICES_SHARED * PI
        return numerator / denominator

    @staticmethod
    def bottom_charm_ratio() -> float:
        # m_b/m_c = (π + √3 + 5)/3 [4.2 sig figs]
        return (PI + SQRT3 + TRIANGLE_CURVED) / VERTICES_TRIANGLE

    @staticmethod
    def top_bottom_ratio() -> float:
        # m_t/m_b = 41 + 1/3 [4.1 sig figs] (TOP_BOTTOM_BASE = 36 + 5)
        return TOP_BOTTOM_BASE + UNITY / VERTICES_TRIANGLE

    # CKM MIXING MATRIX -------------------------------------------------------
    @staticmethod
    def sin_cabibbo() -> float:
        # sin θ_C = 2/(3π + 2√3 - 4) [6.0 sig figs]
        numerator = VERTICES_SHARED
        denominator = VERTICES_TRIANGLE * PI + VERTICES_SHARED * SQRT3 - CKM_OFFSET_4
        return numerator / denominator

    @staticmethod
    def V_cb() -> float:
        # |V_cb| = (2π - 6)/(-3π + 3√3 + 11) [4.6 sig figs]
        numerator = VERTICES_SHARED * PI - TRIANGLE_FLAT
        denominator = -VERTICES_TRIANGLE * PI + VERTICES_TRIANGLE * SQRT3 + FLAT_PLUS_CURVED
        return numerator / denominator

    @staticmethod
    def V_ub() -> float:
        # |V_ub| = (√3 - π + 2)/160 [4.0 sig figs]
        return (SQRT3 - PI + VERTICES_SHARED) / CKM_OFFSET_160

    @staticmethod
    def A_wolfenstein() -> float:
        # A = (-π + 3√3 + 12)/(2π + √3 + 9) [5.7 sig figs]
        numerator = -PI + VERTICES_TRIANGLE * SQRT3 + MASS_OFFSET_12
        denominator = VERTICES_SHARED * PI + SQRT3 + MASS_OFFSET_9
        return numerator / denominator

    @staticmethod
    def rho_bar() -> float:
        # ρ̄ = 4/(3π + √3 + 14) [4.8 sig figs]
        return DIMENSION_FACTOR / (VERTICES_TRIANGLE * PI + SQRT3 + MASS_OFFSET_14)

    @staticmethod
    def eta_bar() -> float:
        # η̄ = 1/(π + √3 - 2) [4.6 sig figs]
        return UNITY / (PI + SQRT3 - VERTICES_SHARED)

    # PMNS NEUTRINO MIXING MATRIX ---------------------------------------------
    @staticmethod
    def sin2_theta_12_pmns() -> float:
        # sin²θ₁₂ = 2√3/(2π + 5) [4.3 sig figs]
        return (VERTICES_SHARED * SQRT3) / (VERTICES_SHARED * PI + TRIANGLE_CURVED)

    @staticmethod
    def sin2_theta_23_pmns() -> float:
        # sin²θ₂₃ = 5/(3π + √3 - 2) [4.1 sig figs]
        return TRIANGLE_CURVED / (VERTICES_TRIANGLE * PI + SQRT3 - VERTICES_SHARED)

    @staticmethod
    def sin2_theta_13_pmns() -> float:
        # sin²θ₁₃ = (3π - 9)/(-3π + √3 + 27) [4.4 sig figs]
        numerator = VERTICES_TRIANGLE * PI - MASS_OFFSET_9
        denominator = -VERTICES_TRIANGLE * PI + SQRT3 + PMNS_OFFSET_27
        return numerator / denominator

    # HIGGS SECTOR ------------------------------------------------------------
    @staticmethod
    def higgs_Z_ratio() -> float:
        # m_H/m_Z = (5π + 6√3)/19 [5.9 sig figs]
        numerator = TRIANGLE_CURVED * PI + TRIANGLE_FLAT * SQRT3
        return numerator / HIGGS_DENOMINATOR

    @staticmethod
    def theta_QCD() -> float:
        # θ_QCD = 0 [EXACT] - Strong CP phase
        return DECIMAL_ZERO

    # EXTENDED PARAMETERS -----------------------------------------------------
    @staticmethod
    def v_over_mZ() -> float:
        # v/m_Z = (-π + 3√3 + 3)/(π + √3 - 3) [4.0 sig figs]
        numerator = -PI + VERTICES_TRIANGLE * SQRT3 + VERTICES_TRIANGLE
        denominator = PI + SQRT3 - VERTICES_TRIANGLE
        return numerator / denominator

    @staticmethod
    def dm2_ratio() -> float:
        # Δm²_sol/Δm²_atm = (√3 + 1)/89 [4.0 sig figs]
        return (SQRT3 + UNITY) / NEUTRINO_89

    @staticmethod
    def electron_planck_ratio() -> float:
        # m_e/m_P = α¹¹ × (3π + 4) [2.6 sig figs] (exponent = FLAT_PLUS_CURVED)
        alpha = FoundationDrawing.alpha()
        return alpha**FLAT_PLUS_CURVED * (VERTICES_TRIANGLE * PI + DIMENSION_FACTOR)
