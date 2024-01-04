# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: OSL-3.0

# pylint: disable=invalid-name

'''
Usage:
Nd-linspace x_start x_stop x_num [y_start y_stop y_num [z_start z_stop z_num [...]]]
Nd-linspace -90 90 181 -180 180 361
Nd-linspace -90 90 181 -180 180 361 0 0 1
Nd-linspace -90 90 721 -180 180 1441
Nd-linspace -90 90 721 -180 180 1441 0 0 1
https://docs.scipy.org/doc/numpy/reference/generated/numpy.linspace.html
'''

__author__ = 'Steven Ward'
__license__ = 'OSL-3.0'

import itertools
import sys

import linspace
import grouper
import remove_exponent

l = []

for g in grouper.grouper(sys.argv[1:], 3):
    (start, stop, num) = g
    a = linspace.linspace(start, stop, num, endpoint=True)
    a = map(remove_exponent.remove_exponent, a)
    #l.append(map(str, a))
    # pylint: disable=consider-using-f-string
    l.append(map('{:f}'.format, a))

for x in itertools.product(*l):
    print(' '.join(x))
