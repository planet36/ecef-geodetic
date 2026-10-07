# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: MPL-2.0

# pylint: disable=invalid-name
# pylint: disable=missing-function-docstring
# pylint: disable=missing-module-docstring

__author__ = 'Steven Ward'
__license__ = 'MPL-2.0'

from decimal import Decimal as D

import numpy as np


# https://numpy.org/doc/stable/reference/generated/numpy.arange.html
def arange(start: D | str, stop: D | str | None = None, step: D | str | None = None,
           endpoint: bool = False) -> np.ndarray:
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

    # NumPy's stubs omit Decimal, which np.arange accepts with an object dtype.
    a = np.arange(start, stop, step, dtype=object) # type: ignore[call-overload]

    if endpoint:
        if len(a) == 0 or a[-1] != stop:
            # Make it a closed interval by appending the stop value.
            a = np.append(a, [stop])

    # pylint: disable=pointless-string-statement
    '''
    while start < stop:
        yield start
        start += step

    if endpoint:
        yield stop
    '''

    return a
