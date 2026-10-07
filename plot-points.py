# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: MPL-2.0

# pylint: disable=missing-function-docstring
# pylint: disable=missing-module-docstring
# pylint: disable=no-else-return
# pylint: disable=bad-indentation
# pylint: disable=fixme
# pylint: disable=invalid-name
# pylint: disable=pointless-string-statement
# pylint: disable=trailing-newlines

__author__ = 'Steven Ward'
__license__ = 'MPL-2.0'
__version__ = '2024-01-08'

import sys
from enum import Enum, auto, unique
from types import FrameType
from typing import NoReturn

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from matplotlib.patches import Ellipse

from ellipsoid_np import WGS84

'''

Examples:

# region 0
python3 Nd-arange.py 0 90 0.0001 0 0 1 |
python3 plot-points.py -v -g --ell --evo --lim --km

# region 1
python3 Nd-arange.py 0 90 0.001 -100 1000 100 |
python3 plot-points.py -v -g --ell --evo --lim --km

# region 2
python3 Nd-arange.py 0 90 0.01 -10_000 100_000 1_000 |
python3 plot-points.py -v -g --ell --evo --lim --km

# region 3
python3 Nd-arange.py 0 90 0.1 -1_000_000 10_000_000 10_000 |
python3 plot-points.py -v -g --ell --evo --lim --km

# region 4
python3 Nd-arange.py 0 90 1 -5_000_000 500_000_000 100_000 |
python3 plot-points.py -v -g --ell --evo --lim --km

# negative height, inside evolute
python3 Nd-arange.py 0 90 5 -6_383_000 0 1_000 |
python3 plot-points.py -v -g --ell --evo --km

# negative height, inside evolute
python3 Nd-arange.py 0 90 1 -6_383_000 0 10_000 |
python3 plot-points.py -v -g --ell --evo --km

'''

# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
# pylint: disable=too-many-positional-arguments
def plot_points_2d(points: npt.NDArray[np.float64], plot_ellipse: bool = False,
                   plot_evolute: bool = False, limit_extents: bool = False,
                   dpi: float = plt.rcParams["figure.dpi"], km: bool = False) -> None:

    # https://matplotlib.org/stable/gallery/color/named_colors.html

    #plt.style.use('dark_background')

    # https://matplotlib.org/stable/api/figure_api.html#module-matplotlib.figure
    fig = plt.figure(dpi=dpi)
    # https://matplotlib.org/stable/api/figure_api.html#matplotlib.figure.Figure.gca
    ax = fig.gca()

    # The WGS-84 ellipsoid is used by default.
    a = WGS84.a
    b = WGS84.b

    if km:
        a /= 1000
        b /= 1000

    a2 = a*a
    b2 = b*b

    if plot_ellipse:
        # https://numpy.org/doc/stable/reference/generated/numpy.linspace.html
        #num = 32768
        #t = np.linspace(-np.pi, np.pi, num+1)
        #x = a * np.cos(t)
        #y = b * np.sin(t)
        #ax.plot(x, y, marker=None, color='tab:blue', linestyle='solid')
        # https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Ellipse.html
        ellipse = Ellipse(xy=(0, 0), width=2*a, height=2*b, facecolor='None', edgecolor='tab:blue')
        ax.add_patch(ellipse)

    if plot_evolute:
        # https://numpy.org/doc/stable/reference/generated/numpy.linspace.html
        num = 64
        t = np.linspace(-np.pi, np.pi, num+1)
        # https://mathworld.wolfram.com/EllipseEvolute.html
        x = (a2 - b2) / a * np.cos(t)**3
        y = (b2 - a2) / b * np.sin(t)**3
        ax.plot(x, y, marker=None, color='tab:red', linestyle='solid')

    # https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set_xlabel.html
    # https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set_ylabel.html
    if km:
        ax.set_xlabel('W (km)')
        ax.set_ylabel('Z\n(km)', rotation='horizontal')
    else:
        ax.set_xlabel('W (m)')
        ax.set_ylabel('Z\n(m)', rotation='horizontal')

    # https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set_aspect.html
    ax.set_aspect('equal')

    # https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.margins.html
    ax.margins(x=0, y=0)

    # https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.tick_params.html
    ax.tick_params(rotation=30)

    # https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.axhline.html
    ax.axhline(color='tab:purple', linestyle='--', linewidth=1)

    # https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.axvline.html
    ax.axvline(color='tab:purple', linestyle='--', linewidth=1)

    # XXX: matplotlib: the ',' marker (pixel) doesn't work
    # https://github.com/matplotlib/matplotlib/issues/11460
    # https://stackoverflow.com/questions/39753282/scatter-plot-with-single-pixel-marker-in-matplotlib
    ax.scatter(points[:, 0], points[:, 1], marker='.', s=(72/fig.dpi)**2,
               color='tab:orange', alpha=0.8)

    if limit_extents:
        x_min = points[:, 0].min()
        x_max = points[:, 0].max()
        y_min = points[:, 1].min()
        y_max = points[:, 1].max()

        # https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set_xlim.html
        ax.set_xlim(left=x_min, right=x_max)

        # https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set_ylim.html
        ax.set_ylim(bottom=y_min, top=y_max)

    plt.show()

# pylint: disable=missing-class-docstring
@unique
class INPUT_DATA_FORMAT(Enum):
    ECEF = auto()
    GEODETIC = auto()

def main(argv: list[str] | None = None) -> int:

    # pylint: disable=import-outside-toplevel
    import argparse
    import signal
    from pathlib import Path

    if argv is None:
        argv = sys.argv

    program_name = Path(argv[0]).name

    program_authors = [__author__]

    # pylint: disable=unused-argument
    def signal_handler(signal_num: int, execution_frame: FrameType | None) -> NoReturn:
        print()
        sys.exit(128 + signal_num)

    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    parser = argparse.ArgumentParser(
        prog=program_name,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''Plot 2D points read from stdin.

The default input data format is ECEF and can be changed to Geodetic with the '-g' option.

ECEF (W, Z) data is in meters.  Geodetic (latitude, height) data is in degrees and meters.''')

    version = f"{program_name} {__version__}\nWritten by {', '.join(program_authors)}"

    parser.add_argument('-V', '--version', action='version', version=version)
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print diagnostics.')
    parser.add_argument('-g', dest='input_data_format', action='store_const',
                        const=INPUT_DATA_FORMAT.GEODETIC, default=INPUT_DATA_FORMAT.ECEF,
                        help='Set the input data format to Geodetic instead of ECEF.')
    parser.add_argument('--ell', dest='plot_ellipse', action='store_true',
                        help='Plot the 2D ellipse of the WGS-84 datum.')
    parser.add_argument('--evo', dest='plot_evolute', action='store_true',
                        help='Plot the evolute of the 2D ellipse of the WGS-84 datum.')
    parser.add_argument('--lim', dest='limit_extents', action='store_true',
                        help='Limit the plot area to the extents of the input data.')
    parser.add_argument('--dpi', type=float, default=plt.rcParams["figure.dpi"],
                        help='Specify the DPI value.  (default: %(default)s)')
    parser.add_argument('-k', '--km', action='store_true',
                        help='Convert the input data from meters to kilometers.')

    args = parser.parse_args(argv[1:])

    def print_verbose(s: str) -> None:
        """Print the message if verbose mode is on"""
        if args.verbose:
            print(f"# {s}", file=sys.stderr)

    def print_error(s: str) -> None:
        """Print the error message"""
        print(f"Error: {s}", file=sys.stderr)
        print(f"Try '{program_name} --help' for more information.", file=sys.stderr)

    print_verbose(f'{args=}')

    # https://numpy.org/doc/stable/reference/generated/numpy.loadtxt.html
    points = np.loadtxt(sys.stdin, ndmin=2)

    print_verbose(f'{points=}')
    print_verbose(f'{points.shape=}')

    if points.size == 0:
        print_error('No points were read')
        return 1

    if points.shape[1] != 2:
        print_error(f'Expected 2 columns, got {points.shape[1]}')
        return 1

    if args.input_data_format == INPUT_DATA_FORMAT.GEODETIC:
        points = WGS84.geodetic_2d_to_ecef(points[:,0], points[:,1])
        print_verbose(f'{points=}')
        print_verbose(f'{points.shape=}')

    if args.km:
        print_verbose('(convert m to km)')
        points /= 1000
        print_verbose(f'{points=}')
        print_verbose(f'{points.shape=}')

    plot_points_2d(points, args.plot_ellipse, args.plot_evolute, args.limit_extents, args.dpi,
                   args.km)

    return 0

if __name__ == '__main__':
    sys.exit(main())
