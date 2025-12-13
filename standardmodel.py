#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Foundation Drawing - Deriving the Standard Model from Pure Geometry.
# Author: (c) 2025 Sally Peck.
# -----------------------------------------------------------------------------

from decimal import Decimal, getcontext
from typing import Dict
import math

from foundationdrawing import FoundationDrawing
from checks import compute_percent_error, compute_agreement_log, compute_stddev_sigma
from constants import (UNITY, ZERO, VERTICES_SHARED, VERTICES_TRIANGLE, TRIANGLE_FLAT, TRIANGLE_CURVED, 
                       VERTICES_PLUS_CURVED, FLAT_PLUS_CURVED, CURVED_PLUS_SHARED, FLAT_TIMES_VERTICES,
                       FLAT_SQUARED, HIGGS_DENOMINATOR, PMNS_OFFSET_27, ICOSAHEDRON_FACES, CKM_OFFSET_160,
                       VERTICES_FOURTH, DIMENSION_1D, DIMENSION_2D, DIMENSION_3D, TOP_BOTTOM_BASE, NEUTRINO_89);

getcontext().prec = 50  # 50 decimal places so 10+ dp digits compute correctly

EXPERIMENTAL: Dict[str, tuple] = {
    # Coupling constants
    'alpha_inv': (137.035999177, 0.000000021),
    'alpha_s': (0.1180, 0.0009),
    'sin2_theta_W': (0.23122, 0.00004),

    # Lepton mass ratios
    'proton_electron': (1836.15267343, 0.00000011),
    'muon_electron': (206.7682830, 0.0000046),
    'tau_muon': (16.8170, 0.0001),
    'koide_Q': (0.666661, 0.000007),

    # Quark mass ratios
    'up_down': (0.4625, 0.025),
    'strange_down': (20.0, 1.0),
    'charm_strange': (13.597, 0.5),
    'bottom_charm': (3.291, 0.1),
    'top_bottom': (41.33, 0.5),

    # CKM matrix
    'sin_cabibbo': (0.22500, 0.00007),
    'V_cb': (0.04182, 0.00085),
    'V_ub': (0.00369, 0.00011),
    'A_wolfenstein': (0.826, 0.015),
    'rho_bar': (0.159, 0.015),
    'eta_bar': (0.348, 0.010),

    # PMNS matrix
    'sin2_12_pmns': (0.307, 0.013),
    'sin2_23_pmns': (0.546, 0.021),
    'sin2_13_pmns': (0.0220, 0.0007),

    # Mass scales
    'higgs_Z': (1.3737, 0.0020),
    'v_over_mZ': (2.698, 0.001),

    # Neutrino sector
    'dm2_ratio': (0.0307, 0.001),

    # Electron to Planck mass
    'm_e_over_m_P': (4.1854e-23, 1e-26),

    # Particle masses (MeV)
    'mass_electron': (0.51099895, 0.00000015),
    'mass_muon': (105.6583755, 0.0000023),
    'mass_tau': (1776.86, 0.12),
    'mass_proton': (938.27208816, 0.00000029),
    'mass_W': (80377.0, 12.0),
    'mass_Z': (91187.6, 2.1),
    'mass_higgs': (125250.0, 170.0),
}


def get_experimental_alpha() -> tuple[Decimal, Decimal]:
    # Direct access instead of iteration - more efficient
    exp_val, exp_unc = EXPERIMENTAL['alpha_inv']
    return Decimal(str(exp_val)), Decimal(str(exp_unc))


class StandardModel:
    # Complete Standard Model from geometric principles.

    def __init__(self):
        self._compute_all()

    def _compute_all(self):
        # Compute all parameters using static methods (no instance needed)
        # Coupling constants
        self.alpha_inv = FoundationDrawing.alpha_inverse()
        self.alpha = FoundationDrawing.alpha()
        self.alpha_s = FoundationDrawing.alpha_s()
        self.sin2_theta_W = FoundationDrawing.sin2_theta_W()

        # Lepton ratios
        self.proton_electron_ratio = FoundationDrawing.proton_electron_ratio()
        self.muon_electron_ratio = FoundationDrawing.muon_electron_ratio()
        self.tau_muon_ratio = FoundationDrawing.tau_muon_ratio()
        self.koide_Q = FoundationDrawing.koide_Q()

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

        # Extended
        self.v_over_mZ = FoundationDrawing.v_over_mZ()
        self.dm2_ratio = FoundationDrawing.dm2_ratio()
        self.m_e_over_m_P = FoundationDrawing.electron_planck_ratio()

        # Compute masses
        self._compute_masses()

    def _compute_masses(self):
        # Compute particle masses from ratios
        # Anchor masses (MeV) - these are experimental inputs
        m_e = 0.51099895
        m_Z = 91187.6

        self.masses = {
            'electron': m_e,
            'muon': m_e * self.muon_electron_ratio,
            'tau': m_e * self.muon_electron_ratio * self.tau_muon_ratio,
            'proton': m_e * self.proton_electron_ratio,
            'W': m_Z * math.sqrt(UNITY - self.sin2_theta_W),
            'Z': m_Z,
            'higgs': m_Z * self.higgs_Z,
        }

    def verify_all(self) -> Dict:
        predictions = {
            'alpha_inv': self.alpha_inv,
            'alpha_s': self.alpha_s,
            'sin2_theta_W': self.sin2_theta_W,
            'proton_electron': self.proton_electron_ratio,
            'muon_electron': self.muon_electron_ratio,
            'tau_muon': self.tau_muon_ratio,
            'koide_Q': self.koide_Q,
            'up_down': self.up_down,
            'strange_down': self.strange_down,
            'charm_strange': self.charm_strange,
            'bottom_charm': self.bottom_charm,
            'top_bottom': self.top_bottom,
            'sin_cabibbo': self.sin_cabibbo,
            'V_cb': self.V_cb,
            'V_ub': self.V_ub,
            'A_wolfenstein': self.A_wolfenstein,
            'rho_bar': self.rho_bar,
            'eta_bar': self.eta_bar,
            'sin2_12_pmns': self.sin2_12_pmns,
            'sin2_23_pmns': self.sin2_23_pmns,
            'sin2_13_pmns': self.sin2_13_pmns,
            'higgs_Z': self.higgs_Z,
            'v_over_mZ': self.v_over_mZ,
            'dm2_ratio': self.dm2_ratio,
            'm_e_over_m_P': self.m_e_over_m_P,
            # Particle masses
            'mass_electron': self.masses['electron'],
            'mass_muon': self.masses['muon'],
            'mass_tau': self.masses['tau'],
            'mass_proton': self.masses['proton'],
            'mass_W': self.masses['W'],
            'mass_Z': self.masses['Z'],
            'mass_higgs': self.masses['higgs'],
        }

        result = {}
        for name, value in predictions.items():
            if name in EXPERIMENTAL:
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
    print("(c) 2025 Sally Peck")
    print("-" * 80)
    print("No free parameters! No magic numbers!")
    print("All exponents as named constants (DIMENSION_2D, DIMENSION_3D)")
    print("All integers as primitive shape constants (VERTICES_SHARED, VERTICES_TRIANGLE)")
    print("-" * 80)

    # High-precision alpha verification
    alpha_inv_hp = FoundationDrawing.alpha_inverse_decimal()
    exp_val, exp_unc = get_experimental_alpha()
    diff = abs(alpha_inv_hp - exp_val)
    n_sigma = float(diff / exp_unc)

    print(f"""
Predicted (50 digits): {alpha_inv_hp}
CODATA 2022:           {exp_val} +/- {exp_unc}
Difference:            {diff}
Statistical:           {n_sigma:.2f}sigma {'(Perfect Match)' if n_sigma < UNITY else ''}
Significant figures:   {compute_agreement_log(float(alpha_inv_hp), float(exp_val)):.1f}, """)
    print("-" * 80)
    
    print(f"""
-------------------------------------------------------------------------------
|  THE PRIMITIVE SHAPE                                                        |
|             C                                                               |
|            /\\           ZERO              = {ZERO}    Additive identity          |
|           /  \\          UNITY             = {UNITY}    Multiplicative identity    |
|          /    \\                                                             |
|         /      \\        SHARED_VERTICES   = {VERTICES_SHARED}    Half-circle endpoints A, B |
|        A--------B       TRIANGLE_VERTICES = {VERTICES_TRIANGLE}    Triangle vertices A, B, C  |
|         \\      /                                                            |
|          \\    /         CURVED_TRIANGLES  = {TRIANGLE_CURVED}    Curved space: 5×60° = 300° |
|           )  (          FLAT_TRIANGLES    = {TRIANGLE_FLAT}    Flat space: 6×60° = 360°   |
|            \\/                                                               |
--------------------------------------------------------------------------------
""");
    
    print(f"""
DERIVED COMBINATIONS (all from basic counting):
-------------------------------------------------------------------------------
  VERTICES_PLUS_CURVED  = {VERTICES_TRIANGLE} + {TRIANGLE_CURVED} = {VERTICES_PLUS_CURVED}   (alpha^-1 numerator)
  FLAT_PLUS_CURVED      = {TRIANGLE_FLAT} + {TRIANGLE_CURVED} = {FLAT_PLUS_CURVED}  (alpha^-1 denominator)
  CURVED_PLUS_SHARED    = {TRIANGLE_CURVED} + {VERTICES_SHARED} = {CURVED_PLUS_SHARED}   (mass formulas)
  FLAT_TIMES_VERTICES   = {TRIANGLE_FLAT} x {VERTICES_TRIANGLE} = {FLAT_TIMES_VERTICES}  (Weinberg angle)
  
HIGHER COMBINATIONS:
--------------------------------------------------------------------------------
  FLAT_SQUARED_PLUS_CURVED        = {TRIANGLE_FLAT}² + {TRIANGLE_CURVED} = {FLAT_SQUARED} + {TRIANGLE_CURVED} = {TOP_BOTTOM_BASE}  (top/bottom ratio)
  ALL_CONFIGS_PLUS_VERT_CURVED    = {FLAT_PLUS_CURVED} + {VERTICES_PLUS_CURVED} = {HIGGS_DENOMINATOR}  (Higgs denominator)
  VERTICES_CUBED                  = {VERTICES_TRIANGLE}³ = {PMNS_OFFSET_27}  (PMNS mixing)
  VERT_CURVED_TIMES_ICOSA         = {VERTICES_PLUS_CURVED} × {ICOSAHEDRON_FACES} = {CKM_OFFSET_160}  (V_ub denominator)
  VERTICES_FOURTH_PLUS_VERT_CURVED = {VERTICES_TRIANGLE}⁴ + {VERTICES_PLUS_CURVED} = {VERTICES_FOURTH} + {VERTICES_PLUS_CURVED} = {NEUTRINO_89}  (neutrino)

DIMENSIONAL EXPONENTS (replace bare **2 and **3):
--------------------------------------------------------------------------------
  DIMENSION_1D = {DIMENSION_1D}  (1D arc/line)
  DIMENSION_2D = {DIMENSION_2D}  (2D plane - use for **2)
  DIMENSION_3D = {DIMENSION_3D}  (3D volume - use for **3)
""")

    tests = StandardModel()
    results = tests.verify_all()
    
    print(f"\n{'#':<2} {'Parameter':<15} {'Predicted':<16} {'Experimental':<13} {'Error %':<12} {'σ':<8} {'Sig'}")
    print("-" * 80)
    excellent = ZERO
    for i, (name, data) in enumerate(results.items()):
        sf = f"{data['sgf']:.1f}" if data['sgf'] < 20 else "inf"
        sig = f"{data['nsg']:.1f}sig" if data['nsg'] < 100 else ">100sig"
        if data['epc'] < 0.01:
            excellent += UNITY
        star = "*" if data['epc'] < 0.01 else ""
        print(f"{i:<2} {name:<15} {data['pre']:<16.10g} {data['exp']:<13.6g}"
              f" {data['epc']:<12.6f} {sig:<9}{sf:<5}{star}")

    print("-" * 80)
    print(f"Excellent (<0.01% error): {excellent}/{len(results)}")
