#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Foundation Drawing - Deriving the Standard Model from Pure Geometry.
# Author: (c) 2025 Sally Peck.
# https://github.com/sallyp/foundationdrawing/
# -----------------------------------------------------------------------------

# Higher-order corrections to alpha^-1 using Lucas numbers
from foundationdrawing import FoundationDrawing

# -----------------------------------------------------------------------------
# Formatting utilities
# -----------------------------------------------------------------------------
def format_decimal_spaced(value, decimal_groups=None):
    """
    Format a number with spaces every 6 digits after the decimal point.
    Example: 137.035999166 -> '137.035999 166'

    Args:
        value: Number to format (Decimal, float, or str)
        decimal_groups: Number of 6-digit groups to show (None = all)
    """
    s = str(value)
    if '.' not in s:
        return s

    integer_part, decimal_part = s.split('.')

    # Group decimal digits in six's
    groups = [decimal_part[i:i+6] for i in range(0, len(decimal_part), 6)]

    if decimal_groups is not None:
        groups = groups[:decimal_groups]

    return integer_part + '.' + ' '.join(groups)


# Use the method from FoundationDrawing
alpha_inverse = FoundationDrawing.alpha_inverse_higher_order


if __name__ == "__main__":
    # Order 4 matches Fan et al. 2023 to 12+ significant figures
    for i in range(9):
        value = round(alpha_inverse(i), 13)
        print(f"Alpha Order {i}: {value}")
