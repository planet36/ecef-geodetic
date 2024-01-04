# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: OSL-3.0

# pylint: disable=invalid-name
# pylint: disable=missing-function-docstring
# pylint: disable=missing-module-docstring

__author__ = 'Steven Ward'
__license__ = 'OSL-3.0'

from decimal import Decimal as D
import numpy as np

# https://numpy.org/doc/stable/reference/generated/numpy.arange.html
def arange(start: D, stop: D = None, step: D = None, endpoint: bool = False):
    start = D(start)

    if stop is None:
        stop = start
        start = D(0) # default value in numpy.arange
    else:
        stop = D(stop)

    if step is None:
        step = D(1) # default value in numpy.arange
    else:
        step = D(step)

    if step == 0:
        raise ValueError(f'Step ({step}) must be non-zero')

    a = np.arange(start, stop, step)

    if endpoint:
        if len(a) == 0 or a[-1] != stop:
            # Make it a closed interval by appending the stop value.
            a = np.append(a, stop)

    # pylint: disable=pointless-string-statement
    '''
    while start < stop:
        yield start
        start += step

    if endpoint:
        yield stop
    '''

    return a
