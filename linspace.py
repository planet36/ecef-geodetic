# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: OSL-3.0

# pylint: disable=fixme
# pylint: disable=invalid-name
# pylint: disable=missing-function-docstring
# pylint: disable=missing-module-docstring

__author__ = 'Steven Ward'
__license__ = 'OSL-3.0'

from decimal import Decimal as D
import numpy as np

# https://numpy.org/doc/stable/reference/generated/numpy.linspace.html
# adapted from
# https://github.com/numpy/numpy/blob/main/numpy/_core/function_base.py#L26
# The start and stop values in numpy.linspace must be float type.
#def linspace(start, stop, num=50, endpoint=True):
def linspace(start: D, stop: D, num: int = 50, endpoint: bool = True):
    start = D(start)
    # XXX: stop is mandatory
    stop = D(stop)

    if num is None:
        num = 50 # default value in numpy.linspace
    else:
        num = int(num)

    if num < 0:
        raise ValueError(f'Number of samples ({num}) must be non-negative')

    #div = (num - 1) if endpoint else num
    div = num - endpoint

    y = np.arange(0, num, dtype=D)

    delta = stop - start

    if num > 1:
        step = delta / div
        y *= step
    #else:
    #    # 0 and 1 item long sequences have an undefined step
    #    step = np.nan
    #    # Multiply with delta to allow possible override of output class.
    #    y *= delta

    y += start

    if endpoint and num > 1:
        y[-1] = stop

    return y
