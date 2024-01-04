# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: OSL-3.0

# pylint: disable=invalid-name
# pylint: disable=no-else-return

##### TODO-COMPILER: use itertools.batched when python 3.12 is available
# https://docs.python.org/3/library/itertools.html?highlight=itertools#itertools.batched

'''
This module contains the grouper function from the Python itertools recipes.
Identical to <https://more-itertools.readthedocs.io/en/stable/api.html#more_itertools.grouper>.
Copied from <https://docs.python.org/3/library/itertools.html#itertools-recipes>.
'''

__author__ = 'Steven Ward'
__license__ = 'OSL-3.0'

from itertools import zip_longest

def grouper(iterable, n, *, incomplete='fill', fillvalue=None):
    "Collect data into non-overlapping fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, fillvalue='x') --> ABC DEF Gxx
    # grouper('ABCDEFG', 3, incomplete='strict') --> ABC DEF ValueError
    # grouper('ABCDEFG', 3, incomplete='ignore') --> ABC DEF
    args = [iter(iterable)] * n
    if incomplete == 'fill':
        return zip_longest(*args, fillvalue=fillvalue)
    if incomplete == 'strict':
        return zip(*args, strict=True)
    if incomplete == 'ignore':
        return zip(*args)
    else:
        raise ValueError('Expected fill, strict, or ignore')
