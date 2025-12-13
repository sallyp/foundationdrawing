#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Foundation Drawing - Deriving the Standard Model from Pure Geometry.
# Author: (c) 2025 Sally Peck.
# -----------------------------------------------------------------------------

from decimal import getcontext
import math
from constants import ZERO, PERCENT_FACTOR

getcontext().prec = 50  # 50 decimal places so 10+ dp digits compute correctly

def compute_agreement_log(predicted: float, actual: float) -> float:
    # check agreement using log - returns significant figures of agreement
    if actual == ZERO:
        return float('inf') if predicted == ZERO else float('-inf')
    if predicted == actual:
        return float('inf')
    relative_error = abs((predicted - actual) / actual)
    if relative_error == ZERO:
        return float('inf')
    return -math.log10(relative_error)

def compute_percent_error(predicted: float, actual: float) -> float:
    # compute error as percentage
    if actual == ZERO:
        return float('inf') if predicted != ZERO else ZERO
    return abs(predicted - actual) / actual * PERCENT_FACTOR

def compute_stddev_sigma(predicted: float, exp_value: float, exp_uncertainty: float) -> float:
    # compute standard deviations from experiment
    if exp_uncertainty == ZERO:
        return float('inf') if predicted != exp_value else ZERO
    return abs(predicted - exp_value) / exp_uncertainty
