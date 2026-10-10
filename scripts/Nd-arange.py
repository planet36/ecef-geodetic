# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: MPL-2.0

# pylint: disable=invalid-name

'''
Usage:
Nd-arange x_start x_stop x_step [y_start y_stop y_step [z_start z_stop z_step [...]]]
Nd-arange -90 90 1 -180 180
Nd-arange -90 90 1 -180 180 1 0
Nd-arange -90 90 0.25 -180 180 0.25
Nd-arange -90 90 0.25 -180 180 0.25 0
https://docs.scipy.org/doc/numpy/reference/generated/numpy.arange.html
'''

__author__ = 'Steven Ward'
__license__ = 'MPL-2.0'

import itertools
import sys

import arange
import remove_exponent

l = []

# Group the arguments in threes, and pad a short last group with None.
for g in itertools.zip_longest(*[iter(sys.argv[1:])] * 3):
    (start, stop, step) = g
    a = arange.arange(start, stop, step, endpoint=True)
    d = map(remove_exponent.remove_exponent, a)
    #l.append(map(str, d))
    # pylint: disable=consider-using-f-string
    l.append(map('{:f}'.format, d))

for x in itertools.product(*l):
    print(' '.join(x))
