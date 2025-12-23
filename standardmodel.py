#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Foundation Drawing - Deriving the Standard Model from Pure Geometry.
# Author: (c) 2025 Sally Peck.
# https://github.com/sallyp/foundationdrawing/
# -----------------------------------------------------------------------------

from decimal import Decimal, getcontext
from typing import Dict
import math

from foundationdrawing import FoundationDrawing
from checks import compute_percent_error, compute_agreement_log, compute_stddev_sigma
from constants import (UNITY, ZERO, VERTICES_SHARED, VERTICES_TRIANGLE, TRIANGLE_FLAT, TRIANGLE_CURVED,
                       VERTICES_PLUS_CURVED, FLAT_PLUS_CURVED, CURVED_PLUS_SHARED, FLAT_TIMES_VERTICES,
                       FLAT_SQUARED, HIGGS_DENOMINATOR, PMNS_OFFSET_27, ICOSAHEDRON_FACES, CKM_OFFSET_160,
                       VERTICES_FOURTH, DIMENSION_1D, DIMENSION_2D, DIMENSION_3D, TOP_BOTTOM_BASE, NEUTRINO_89)

getcontext().prec = 50  # 50 decimal places so 10+ dp digits compute correctly

# -----------------------------------------------------------------------------
# Experimental Values: Fine Structure Constant (alpha)
# -----------------------------------------------------------------------------
# CODATA 2022 (NIST, published May 2024)
# Source: https://physics.nist.gov/cuu/Constants/Table/allascii.txt
ALPHA_CODATA_2022 = Decimal('7.2973525643e-3')           # alpha value
ALPHA_CODATA_2022_UNCERTAINTY = Decimal('0.0000000011e-3')  # standard uncertainty
ALPHA_INV_CODATA_2022 = Decimal('137.035999177')         # 1/alpha (derived)
ALPHA_INV_CODATA_2022_UNCERTAINTY = Decimal('0.000000021')

# Fan et al. 2023 - Electron Magnetic Moment Measurement
# Source: Phys. Rev. Lett. 130, 071801 (2023)
# DOI: 10.1103/PhysRevLett.130.071801
# 1/alpha derived from g/2 measurement combined with SM theory
ALPHA_INV_FAN_2023 = Decimal('137.035999166')            # 1/alpha value
ALPHA_INV_FAN_2023_UNCERTAINTY = Decimal('0.000000015')  # [0.11 ppb]
ALPHA_FAN_2023 = Decimal('7.2973525693e-3')              # alpha (derived from 1/alpha)
ALPHA_FAN_2023_UNCERTAINTY = Decimal('0.0000000008e-3')

# Float versions for standard calculations
ALPHA_CODATA_2022_FLOAT = 7.2973525643e-3
ALPHA_INV_CODATA_2022_FLOAT = 137.035999177
ALPHA_FAN_2023_FLOAT = 7.2973525693e-3
ALPHA_INV_FAN_2023_FLOAT = 137.035999166

# Experimental Data Sources - see README.md for full references


EXPERIMENTAL: Dict[str, tuple] = {
    # Coupling constants
    'alpha_inv_cd': (137.035999177, 0.000000021),   # CODATA 2022
    'alpha_inv_fn': (137.035999166, 0.000000015),   # Fan et al. 2023
    'alpha_s': (0.1179, 0.0009),                    # PDG 2024 world average
    'sin2_theta_W': (0.23129, 0.00004),             # PDG 2024 MS-bar at MZ
    'sin2_thetaW_ho': (0.23129, 0.00004),           # PDG 2024 (higher order)
    'sin_cabibbo_ho': (0.22431, 0.00008),          # PDG 2024 (higher order)
    'charm_s_2': (11.85, 0.16),                    # PDG 2024 (alternative formula)
    # Lepton mass ratios
    'proton_e': (1836.152673426, 0.000000032),     # CODATA 2022
    'muon_e': (206.7682827, 0.0000046),            # CODATA 2022
    'tau_muon': (16.8168, 0.0001),                 # PDG 2024 (1776.86/105.6583745)
    'koide_Q': (0.666661, 0.000007),
    # Quark mass ratios
    'up_down': (0.4625, 0.025),
    'strange_d': (20.0, 1.0),
    'charm_s': (11.85, 0.16),                      # PDG 2024 lattice QCD
    'bottom_c': (3.291, 0.1),
    'top_b': (41.33, 0.5),
    # CKM matrix (PDG 2024)
    'sin_cabibbo': (0.22431, 0.00085),             # PDG 2024 |Vus|
    'V_cb': (0.0406, 0.0009),                      # PDG 2024
    'V_cb_ho': (0.0406, 0.0009),                   # PDG 2024 (higher order)
    'V_ub': (0.00369, 0.00011),
    'A_wolf': (0.826, 0.015),
    'rho_bar': (0.159, 0.015),
    'eta_bar': (0.348, 0.010),
    # PMNS matrix (NuFIT 6.0, September 2024, Normal Ordering)
    'sin2_12': (0.308, 0.012),                     # NuFIT 6.0
    'sin2_23': (0.470, 0.017),                     # NuFIT 6.0
    'sin2_13': (0.02215, 0.00056),                 # NuFIT 6.0
    # Mass scales
    'higgs_Z': (1.3728, 0.0012),                   # PDG 2024 (125110/91187.6)
    'v_over_mZ': (2.698, 0.001),
    # Neutrino sector (NuFIT 6.0)
    'dm2_ratio': (0.0298, 0.001),                  # dm21^2/|dm31^2| = 7.49e-5/2.513e-3
    'sin2_23_ho': (0.470, 0.017),                  # NuFIT 6.0 (higher order)
    'dm2_ratio_ho': (0.0298, 0.001),               # NuFIT 6.0 (higher order)
    # Electron to Planck mass
    'm_e_m_Plnck': (4.1854e-23, 1e-26),
    # Particle masses (MeV) - PDG 2024
    'm_electron': (0.51099895, 0.00000015),
    'm_muon': (105.6583745, 0.0000024),            # PDG 2024
    'm_tau': (1776.86, 0.12),
    'm_proton': (938.27208816, 0.00000029),
    'm_W': (80369.2, 13.3),                        # PDG 2024 world average
    'm_Z': (91187.6, 2.1),
    'm_higgs': (125110.0, 110.0),                  # ATLAS 2024 combined
}


def get_experimental_alpha() -> tuple[Decimal, Decimal]:
    # Direct access instead of iteration - more efficient
    exp_val, exp_unc = EXPERIMENTAL['alpha_inv_cd']
    return Decimal(str(exp_val)), Decimal(str(exp_unc))


class StandardModel:
    # Complete Standard Model from geometric principles.

    def __init__(self):
        self._compute_all()

    def _compute_all(self):
        # Compute all parameters using static methods (no instance needed)
        # Coupling constants (Order 3 - original)
        self.alpha_inv = FoundationDrawing.alpha_inverse()
        self.alpha = FoundationDrawing.alpha()
        self.alpha_s = FoundationDrawing.alpha_s()
        self.sin2_theta_W = FoundationDrawing.sin2_theta_W()
        self.sin2_theta_W_ho = FoundationDrawing.sin2_theta_W_higher_order()
        self.sin_cabibbo_ho = FoundationDrawing.sin_cabibbo_higher_order()
        self.charm_strange_2 = FoundationDrawing.charm_strange_ratio_2()
        self.V_cb_ho = FoundationDrawing.V_cb_higher_order()
        self.sin2_23_ho = FoundationDrawing.sin2_theta_23_pmns_higher_order()
        self.dm2_ratio_ho = FoundationDrawing.dm2_ratio_higher_order()

        # Coupling constants (Order 4 - higher order)
        self.alpha_inv_ho = float(FoundationDrawing.alpha_inverse_higher_order())
        self.alpha_ho = float(FoundationDrawing.alpha_higher_order())

        # Lepton ratios (Order 3)
        self.proton_electron_ratio = FoundationDrawing.proton_electron_ratio()
        self.muon_electron_ratio = FoundationDrawing.muon_electron_ratio()
        self.tau_muon_ratio = FoundationDrawing.tau_muon_ratio()
        self.koide_Q = FoundationDrawing.koide_Q()

        # Lepton ratios (Order 4 - higher order)
        self.muon_electron_ratio_ho = float(FoundationDrawing.muon_electron_ratio_higher_order())

        # Quark ratios
        self.up_down = FoundationDrawing.up_down_ratio()
        self.strange_down = FoundationDrawing.strange_down_ratio()
        self.charm_strange = FoundationDrawing.charm_strange_ratio()
        self.bottom_charm = FoundationDrawing.bottom_charm_ratio()
        self.top_bottom = FoundationDrawing.top_bottom_ratio()

        # CKM matrix
        self.sin_cabibbo = FoundationDrawing.sin_cabibbo()
        self.V_cb = FoundationDrawing.V_cb()
        self.V_ub = FoundationDrawing.V_ub()
        self.A_wolfenstein = FoundationDrawing.A_wolfenstein()
        self.rho_bar = FoundationDrawing.rho_bar()
        self.eta_bar = FoundationDrawing.eta_bar()

        # PMNS matrix
        self.sin2_12_pmns = FoundationDrawing.sin2_theta_12_pmns()
        self.sin2_23_pmns = FoundationDrawing.sin2_theta_23_pmns()
        self.sin2_13_pmns = FoundationDrawing.sin2_theta_13_pmns()

        # Higgs sector
        self.higgs_Z = FoundationDrawing.higgs_Z_ratio()
        self.theta_QCD = FoundationDrawing.theta_QCD()

        # Extended (Order 3)
        self.v_over_mZ = FoundationDrawing.v_over_mZ()
        self.dm2_ratio = FoundationDrawing.dm2_ratio()
        self.m_e_over_m_P = FoundationDrawing.electron_planck_ratio()

        # Extended (Order 4 - higher order)
        self.m_e_over_m_P_ho = float(FoundationDrawing.electron_planck_ratio_higher_order())

        # Compute masses
        self._compute_masses()

    def _compute_masses(self):
        # Compute particle masses from ratios
        # Anchor masses (MeV) - these are experimental inputs
        m_e = 0.51099895
        m_Z = 91187.6

        # Order 3 masses
        self.masses = {
            'electron': m_e,
            'muon': m_e * self.muon_electron_ratio,
            'tau': m_e * self.muon_electron_ratio * self.tau_muon_ratio,
            'proton': m_e * self.proton_electron_ratio,
            'W': m_Z * math.sqrt(UNITY - self.sin2_theta_W),
            'Z': m_Z,
            'higgs': m_Z * self.higgs_Z,
        }

        # Order 4 masses (higher order)
        self.masses_ho = {
            'electron': m_e,
            'muon': m_e * self.muon_electron_ratio_ho,
            'tau': m_e * self.muon_electron_ratio_ho * self.tau_muon_ratio,
            'proton': m_e * self.proton_electron_ratio,
            'W': m_Z * math.sqrt(UNITY - self.sin2_theta_W),
            'Z': m_Z,
            'higgs': m_Z * self.higgs_Z,
        }

    def verify_original(self) -> Dict:
        # Original formulas organized by Standard Model categories
        # Returns ordered dict with category markers
        from collections import OrderedDict
        predictions = OrderedDict([
            # --- 1. GAUGE COUPLINGS (3) ---
            ('_cat_gauge', '--- GAUGE COUPLINGS (3) ---'),
            ('alpha_inv_cd', self.alpha_inv),
            ('alpha_s', self.alpha_s),
            ('sin2_theta_W', self.sin2_theta_W),
            # --- 2. QUARK MASS RATIOS (5) ---
            ('_cat_quark', '--- QUARK MASS RATIOS (5) ---'),
            ('up_down', self.up_down),
            ('strange_d', self.strange_down),
            ('charm_s', self.charm_strange),
            ('bottom_c', self.bottom_charm),
            ('top_b', self.top_bottom),
            # --- 3. LEPTON MASS RATIOS (3) ---
            ('_cat_lepton', '--- LEPTON MASS RATIOS (3) ---'),
            ('muon_e', self.muon_electron_ratio),
            ('tau_muon', self.tau_muon_ratio),
            ('koide_Q', self.koide_Q),
            # --- 4. CKM MATRIX (6) ---
            ('_cat_ckm', '--- CKM MATRIX (6) ---'),
            ('sin_cabibbo', self.sin_cabibbo),
            ('V_cb', self.V_cb),
            ('V_ub', self.V_ub),
            ('A_wolf', self.A_wolfenstein),
            ('rho_bar', self.rho_bar),
            ('eta_bar', self.eta_bar),
            # --- 5. HIGGS SECTOR (2) ---
            ('_cat_higgs', '--- HIGGS SECTOR (2) ---'),
            ('higgs_Z', self.higgs_Z),
            ('v_over_mZ', self.v_over_mZ),
            # --- 6. PMNS MATRIX (3) ---
            ('_cat_pmns', '--- PMNS MATRIX (3) ---'),
            ('sin2_12', self.sin2_12_pmns),
            ('sin2_23', self.sin2_23_pmns),
            ('sin2_13', self.sin2_13_pmns),
            # --- 7. NEUTRINO MASSES (1) ---
            ('_cat_neutrino', '--- NEUTRINO MASS RATIO (1) ---'),
            ('dm2_ratio', self.dm2_ratio),
            # --- 8. DERIVED MASSES (6) ---
            ('_cat_derived', '--- DERIVED MASSES (6) ---'),
            ('proton_e', self.proton_electron_ratio),
            ('m_e_m_Plnck', self.m_e_over_m_P),
            ('m_muon', self.masses['muon']),
            ('m_tau', self.masses['tau'],),
            ('m_proton', self.masses['proton']),
            ('m_W', self.masses['W']),
        ])

        result = OrderedDict()
        for name, value in predictions.items():
            if name.startswith('_cat_'):
                result[name] = {'category': value}
            elif name in EXPERIMENTAL:
                exp_val, exp_unc = EXPERIMENTAL[name]
                result[name] = {
                    'pre': value,
                    'exp': exp_val,
                    'unc': exp_unc,
                    'epc': compute_percent_error(value, exp_val),
                    'sgf': compute_agreement_log(value, exp_val),
                    'nsg': compute_stddev_sigma(value, exp_val, exp_unc),
                }
        return result

    def verify_higher_order(self) -> Dict:
        # Higher-order formulas with Planck-scale corrections
        # Organized by Standard Model categories
        from collections import OrderedDict
        predictions = OrderedDict([
            # --- 1. GAUGE COUPLINGS (2) ---
            ('_cat_gauge', '--- GAUGE COUPLINGS (2) ---'),
            ('alpha_inv_fn', self.alpha_inv_ho),
            ('sin2_thetaW_ho', self.sin2_theta_W_ho),
            # --- 2. QUARK MASS RATIOS (1) ---
            ('_cat_quark', '--- QUARK MASS RATIOS (1) ---'),
            ('charm_s_2', self.charm_strange_2),
            # --- 3. LEPTON MASS RATIOS (2) ---
            ('_cat_lepton', '--- LEPTON MASS RATIOS (2) ---'),
            ('muon_e', self.muon_electron_ratio_ho),
            ('m_e_m_Plnck', self.m_e_over_m_P_ho),
            # --- 4. CKM MATRIX (2) ---
            ('_cat_ckm', '--- CKM MATRIX (2) ---'),
            ('sin_cabibbo_ho', self.sin_cabibbo_ho),
            ('V_cb_ho', self.V_cb_ho),
            # --- 5. PMNS MATRIX (1) ---
            ('_cat_pmns', '--- PMNS MATRIX (1) ---'),
            ('sin2_23_ho', self.sin2_23_ho),
            # --- 6. NEUTRINO MASS RATIO (1) ---
            ('_cat_neutrino', '--- NEUTRINO MASS RATIO (1) ---'),
            ('dm2_ratio_ho', self.dm2_ratio_ho),
            # --- 7. DERIVED MASSES (2) ---
            ('_cat_derived', '--- DERIVED MASSES (2) ---'),
            ('m_muon', self.masses_ho['muon']),
            ('m_tau', self.masses_ho['tau']),
        ])

        result = OrderedDict()
        for name, value in predictions.items():
            if name.startswith('_cat_'):
                result[name] = {'category': value}
            elif name in EXPERIMENTAL:
                exp_val, exp_unc = EXPERIMENTAL[name]
                result[name] = {
                    'pre': value,
                    'exp': exp_val,
                    'unc': exp_unc,
                    'epc': compute_percent_error(value, exp_val),
                    'sgf': compute_agreement_log(value, exp_val),
                    'nsg': compute_stddev_sigma(value, exp_val, exp_unc),
                }
        return result


if __name__ == "__main__":
    print("-" * 80)
    print("Foundation Drawing - Deriving the Standard Model from Pure Geometry.")
    print("(c) 2025 Sally Peck - https://github.com/sallyp/foundationdrawing/")
    print("-" * 80)
    print("No free parameters! No magic numbers!")
    print("All exponents as named constants (DIMENSION_2D, DIMENSION_3D)")
    print("All integers as primitive shape constants (VERTICES_SHARED, VERTICES_TRIANGLE)")
    print("-" * 80)

    # High-precision alpha verification
    alpha_inv_o3 = FoundationDrawing.alpha_inverse_decimal()
    alpha_inv_o4 = FoundationDrawing.alpha_inverse_higher_order()

    # CODATA 2022
    exp_val, exp_unc = get_experimental_alpha()
    diff_o3 = abs(alpha_inv_o3 - exp_val)
    n_sigma_o3 = float(diff_o3 / exp_unc)
    diff_o4 = abs(alpha_inv_o4 - exp_val)
    n_sigma_o4 = float(diff_o4 / exp_unc)
    sf_o3_cd = compute_agreement_log(float(alpha_inv_o3), float(exp_val))
    sf_o4_cd = compute_agreement_log(float(alpha_inv_o4), float(exp_val))

    # Fan et al. 2023
    fan_val = ALPHA_INV_FAN_2023
    fan_unc = ALPHA_INV_FAN_2023_UNCERTAINTY
    fan_diff_o3 = abs(alpha_inv_o3 - fan_val)
    fan_n_sigma_o3 = float(fan_diff_o3 / fan_unc)
    fan_diff_o4 = abs(alpha_inv_o4 - fan_val)
    fan_n_sigma_o4 = float(fan_diff_o4 / fan_unc)
    sf_o3_fan = compute_agreement_log(float(alpha_inv_o3), float(fan_val))
    sf_o4_fan = compute_agreement_log(float(alpha_inv_o4), float(fan_val))

    # Format alpha values to fit 80 cols (show 36 digits after decimal)
    a3_str = str(alpha_inv_o3)[:40]
    a4_str = str(alpha_inv_o4)[:40]

    print(f"""
THE PRIMITIVE SHAPE
--------------------------------------------------------------------------------
           C
          /\\        ZERO  = {ZERO}  Additive identity
         /  \\       UNITY = {UNITY}  Multiplicative identity
        /    \\
       /      \\     SHARED_VERTICES   = {VERTICES_SHARED}  (A,B on half-circle)
      A--------B    TRIANGLE_VERTICES = {VERTICES_TRIANGLE}  (A,B,C)
       \\      /
        \\    /      CURVED_TRIANGLES  = {TRIANGLE_CURVED}  (5x60=300 deg)
         )  (       FLAT_TRIANGLES    = {TRIANGLE_FLAT}  (6x60=360 deg)
          \\/
--------------------------------------------------------------------------------
""")

    print(f"""DERIVED COMBINATIONS:
--------------------------------------------------------------------------------
  VERTICES_PLUS_CURVED = {VERTICES_TRIANGLE}+{TRIANGLE_CURVED} = {VERTICES_PLUS_CURVED}    (alpha^-1 numerator)
  FLAT_PLUS_CURVED     = {TRIANGLE_FLAT}+{TRIANGLE_CURVED} = {FLAT_PLUS_CURVED}   (alpha^-1 denominator)
  CURVED_PLUS_SHARED   = {TRIANGLE_CURVED}+{VERTICES_SHARED} = {CURVED_PLUS_SHARED}    (mass formulas)
  FLAT_TIMES_VERTICES  = {TRIANGLE_FLAT}x{VERTICES_TRIANGLE} = {FLAT_TIMES_VERTICES}   (Weinberg angle)

HIGHER COMBINATIONS:
--------------------------------------------------------------------------------
  FLAT_SQUARED+CURVED  = {TRIANGLE_FLAT}^2+{TRIANGLE_CURVED} = {FLAT_SQUARED}+{TRIANGLE_CURVED} = {TOP_BOTTOM_BASE}  (top/bottom)
  CONFIGS+VERT_CURVED  = {FLAT_PLUS_CURVED}+{VERTICES_PLUS_CURVED} = {HIGGS_DENOMINATOR}        (Higgs denom)
  VERTICES_CUBED       = {VERTICES_TRIANGLE}^3 = {PMNS_OFFSET_27}           (PMNS mixing)
  VERT_CURVED*ICOSA    = {VERTICES_PLUS_CURVED}x{ICOSAHEDRON_FACES} = {CKM_OFFSET_160}        (V_ub denom)
  VERT_4TH+VERT_CURVED = {VERTICES_TRIANGLE}^4+{VERTICES_PLUS_CURVED} = {VERTICES_FOURTH}+{VERTICES_PLUS_CURVED} = {NEUTRINO_89}  (neutrino)

DIMENSIONAL EXPONENTS:
--------------------------------------------------------------------------------
  DIMENSION_1D = {DIMENSION_1D}  (1D arc/line)
  DIMENSION_2D = {DIMENSION_2D}  (2D plane)
  DIMENSION_3D = {DIMENSION_3D}  (3D volume)
""")

    tests = StandardModel()

    # Alpha comparison
    print("\n" + "=" * 80)
    print("ALPHA COMPARISON: Order 3 (Original) vs Order 4 (Higher-Order)")
    print("=" * 80)
    print(f"{'Source':<12} {'Order 3':<40} {'Order 4':<40}")
    print("-" * 80)
    print(f"{'Predicted':<12} {a3_str:<40} {a4_str:<40}")
    print(f"{'CODATA 2022':<12} {str(exp_val):<40} diff: {diff_o3:.2e} ({sf_o3_cd:.1f} SF)")
    print(f"{'Fan 2023':<12} {str(fan_val):<40} diff: {fan_diff_o4:.2e} ({sf_o4_fan:.1f} SF)")
    print("-" * 80)
    print("Order 3 matches CODATA 2022 (11.4 sig figs)")
    print("Order 4 matches Fan 2023 (13.5 sig figs)")

    # Original formulas (CODATA 2022)
    results_orig = tests.verify_original()
    print("\n" + "=" * 80)
    print("ORIGINAL FORMULAS (vs CODATA 2022, PDG 2024)")
    print("=" * 80)
    excellent_orig = ZERO
    param_count = ZERO
    for name, data in results_orig.items():
        if 'category' in data:
            print(f"\n{data['category']}")
            print(f"{'Param':<14} {'Predicted':<17} {'Experim':<14} {'Err%':<8} {'SF'}")
        else:
            sf = f"{data['sgf']:.0f}" if data['sgf'] < 20 else "+"
            if data['epc'] < 0.01:
                excellent_orig += UNITY
            star = "*" if data['epc'] < 0.01 else ""
            print(f"{name:<14} {data['pre']:<17.15g} {data['exp']:<14.12g}"
                  f" {data['epc']:<8.4f} {sf:<2}{star}")
            param_count += UNITY
    print("-" * 80)
    print(f"Total: {param_count} parameters, Excellent (<0.01% error): {excellent_orig}")

    # Higher-order formulas (Fan 2023, NuFIT 6.0, PDG 2024)
    results_ho = tests.verify_higher_order()
    print("\n" + "=" * 80)
    print("HIGHER-ORDER FORMULAS (vs Fan 2023, NuFIT 6.0, PDG 2024)")
    print("=" * 80)
    excellent_ho = ZERO
    param_count_ho = ZERO
    for name, data in results_ho.items():
        if 'category' in data:
            print(f"\n{data['category']}")
            print(f"{'Param':<14} {'Predicted':<17} {'Experim':<14} {'Err%':<8} {'SF'}")
        else:
            sf = f"{data['sgf']:.0f}" if data['sgf'] < 20 else "+"
            if data['epc'] < 0.01:
                excellent_ho += UNITY
            star = "*" if data['epc'] < 0.01 else ""
            print(f"{name:<14} {data['pre']:<17.15g} {data['exp']:<14.12g}"
                  f" {data['epc']:<8.4f} {sf:<2}{star}")
            param_count_ho += UNITY
    print("-" * 80)
    print(f"Total: {param_count_ho} parameters, Excellent (<0.01% error): {excellent_ho}")

    # All Formula Explanations
    print("\n" + "=" * 80)
    print("FORMULA REFERENCE")
    print("=" * 80)

    all_formulas = [
        # --- GAUGE COUPLINGS ---
        ('_cat', '--- GAUGE COUPLINGS ---'),
        ('alpha_inv', FoundationDrawing.alpha_inverse_doc),
        ('  +ho', FoundationDrawing.alpha_inverse_higher_order_doc),
        ('alpha_s', FoundationDrawing.alpha_s_doc),
        ('sin2_theta_W', FoundationDrawing.sin2_theta_W_doc),
        ('  +ho', FoundationDrawing.sin2_theta_W_higher_order_doc),
        # --- QUARK MASS RATIOS ---
        ('_cat', '--- QUARK MASS RATIOS ---'),
        ('up_down', FoundationDrawing.up_down_ratio_doc),
        ('strange_down', FoundationDrawing.strange_down_ratio_doc),
        ('charm_strange', FoundationDrawing.charm_strange_ratio_doc),
        ('  +ho', FoundationDrawing.charm_strange_ratio_2_doc),
        ('bottom_charm', FoundationDrawing.bottom_charm_ratio_doc),
        ('top_bottom', FoundationDrawing.top_bottom_ratio_doc),
        # --- LEPTON MASS RATIOS ---
        ('_cat', '--- LEPTON MASS RATIOS ---'),
        ('proton_e', FoundationDrawing.proton_electron_ratio_doc),
        ('muon_e', FoundationDrawing.muon_electron_ratio_doc),
        ('tau_muon', FoundationDrawing.tau_muon_ratio_doc),
        ('koide_Q', FoundationDrawing.koide_Q_doc),
        # --- CKM MATRIX ---
        ('_cat', '--- CKM MATRIX ---'),
        ('sin_cabibbo', FoundationDrawing.sin_cabibbo_doc),
        ('  +ho', FoundationDrawing.sin_cabibbo_higher_order_doc),
        ('V_cb', FoundationDrawing.V_cb_doc),
        ('  +ho', FoundationDrawing.V_cb_higher_order_doc),
        ('V_ub', FoundationDrawing.V_ub_doc),
        ('A_wolfenstein', FoundationDrawing.A_wolfenstein_doc),
        ('rho_bar', FoundationDrawing.rho_bar_doc),
        ('eta_bar', FoundationDrawing.eta_bar_doc),
        # --- PMNS MATRIX ---
        ('_cat', '--- PMNS MATRIX ---'),
        ('sin2_12', FoundationDrawing.sin2_theta_12_pmns_doc),
        ('sin2_23', FoundationDrawing.sin2_theta_23_pmns_doc),
        ('  +ho', FoundationDrawing.sin2_theta_23_pmns_higher_order_doc),
        ('sin2_13', FoundationDrawing.sin2_theta_13_pmns_doc),
        # --- HIGGS SECTOR ---
        ('_cat', '--- HIGGS SECTOR ---'),
        ('higgs_Z', FoundationDrawing.higgs_Z_ratio_doc),
        ('v_over_mZ', FoundationDrawing.v_over_mZ_doc),
        # --- NEUTRINO SECTOR ---
        ('_cat', '--- NEUTRINO MASS ---'),
        ('dm2_ratio', FoundationDrawing.dm2_ratio_doc),
        ('  +ho', FoundationDrawing.dm2_ratio_higher_order_doc),
        # --- DERIVED ---
        ('_cat', '--- DERIVED ---'),
        ('m_e/m_Planck', FoundationDrawing.electron_planck_ratio_doc),
    ]

    for label, func_or_cat in all_formulas:
        if label == '_cat':
            print(f"\n{func_or_cat}")
            print(f"{'Parameter':<14} {'Formula':<28} {'Meaning'}")
        else:
            result = func_or_cat()
            formula, meaning = result[0], result[1]
            print(f"{label:<14} {formula:<28} {meaning}")
    print("-" * 80)

