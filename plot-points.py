# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: OSL-3.0

# pylint: disable=missing-function-docstring
# pylint: disable=missing-module-docstring
# pylint: disable=no-else-return
# pylint: disable=bad-indentation
# pylint: disable=fixme
# pylint: disable=invalid-name
# pylint: disable=pointless-string-statement
# pylint: disable=trailing-newlines

__author__ = 'Steven Ward'
__license__ = 'OSL-3.0'
__version__ = '2024-01-08'

from enum import Enum, auto, unique
import sys

from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import numpy as np

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
def plot_points_2d(points, plot_ellipse=False, plot_evolute=False,
                   limit_extents=False, dpi=plt.rcParams["figure.dpi"], km=False):

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

# pylint: disable=too-many-branches
# pylint: disable=too-many-statements
def main(argv = None):

    # pylint: disable=import-outside-toplevel
    import getopt
    import os
    import signal

    if argv is None:
        argv = sys.argv

    program_name = os.path.basename(argv[0])

    program_authors = [__author__]

    # default values
    default_verbose = False
    default_input_data_format = INPUT_DATA_FORMAT.ECEF
    default_plot_ellipse = False
    default_plot_evolute = False
    default_limit_extents = False
    default_dpi = plt.rcParams["figure.dpi"]
    default_km = False

    # mutable values
    verbose = default_verbose
    input_data_format = default_input_data_format
    plot_ellipse = default_plot_ellipse
    plot_evolute = default_plot_evolute
    limit_extents = default_limit_extents
    dpi = default_dpi
    km = default_km

    # pylint: disable=unused-argument
    def signal_handler(signal_num, execution_frame):
        print()
        sys.exit(0)

    signal.signal(signal.SIGINT, signal_handler) # Interactive attention signal. (Ctrl-C)
    signal.signal(signal.SIGTERM, signal_handler) # Termination request. (kill default signal)

    def print_version():
        """Print the version information"""
        print(f"{program_name} {__version__}")
        print(f"Written by {', '.join(program_authors)}")

    def print_verbose(s):
        """Print the message if verbose mode is on"""
        if verbose:
            print(f"# {s}", file=sys.stderr)

    # pylint: disable=unused-variable
    def print_warning(s):
        """Print the warning message"""
        print(f"Warning: {s}", file=sys.stderr)

    def print_error(s):
        """Print the error message"""
        print(f"Error: {s}", file=sys.stderr)
        print(f"Try '{program_name} --help' for more information.", file=sys.stderr)

    def print_help():
        """Print the help message"""

        print(f'''Plot 2D points read from stdin.
Usage: {program_name} [OPTION]...

The default input data format is ECEF and can be changed to Geodetic with the '-g' option.

ECEF (W, Z) data is in meters.  Geodetic (latitude, height) data is in degrees and meters.

OPTIONS

-V, --version
    Print the version information, then exit.

-h, --help
    Print this message, then exit.

-v, --verbose
    Print diagnostics.
    (default: {default_verbose})

-g
    Set the input data format to Geodetic instead of ECEF.

--ell
    Plot the 2D ellipse of the WGS-84 datum.

--evo
    Plot the evolute of the 2D ellipse of the WGS-84 datum.

--lim
    Limit the plot area to the extents of the input data.

--dpi DPI
    Specify the DPI value.
    (default: {default_dpi})

-k, --km
    Convert the input data from meters to kilometers.
''')

    short_options = 'Vhvgk'
    long_options = ['version', 'help', 'verbose', 'ell', 'evo', 'lim', 'dpi=', 'km']

    try:
        (options, remaining_args) = getopt.getopt(argv[1:], short_options, long_options)
    except getopt.GetoptError as err:
        print_error(err)
        return 1

    for (option, value) in options:
        try:
            if option in ['-V', '--version']:
                print_version()
                return 0
            elif option in ['-h', '--help']:
                print_help()
                return 0
            elif option in ['-v', '--verbose']:
                verbose = True
            elif option in ['-g']:
                input_data_format = INPUT_DATA_FORMAT.GEODETIC
            elif option in ['--ell']:
                plot_ellipse = True
            elif option in ['--evo']:
                plot_evolute = True
            elif option in ['--lim']:
                limit_extents = True
            elif option in ['--dpi']:
                dpi = float(value)
            elif option in ['-k', '--km']:
                km = True
            else:
                print_error(f"Unhandled option: {option}")
                return 1

        except ValueError as err:
            print_error(err)
            return 1

    print_verbose(f'{remaining_args=}')

    print_verbose(f'{input_data_format=}')
    print_verbose(f'{plot_ellipse=}')
    print_verbose(f'{plot_evolute=}')
    print_verbose(f'{limit_extents=}')
    print_verbose(f'{dpi=}')

    # https://numpy.org/doc/stable/reference/generated/numpy.loadtxt.html
    points = np.loadtxt(sys.stdin)

    print_verbose(f'{points=}')
    print_verbose(f'{points.shape=}')

    if points.size == 0: # empty
        return 1

    if input_data_format == INPUT_DATA_FORMAT.GEODETIC:
        points = WGS84.geodetic_2d_to_ecef(points[:,0], points[:,1])
        print_verbose(f'{points=}')
        print_verbose(f'{points.shape=}')

    if km:
        print_verbose('(convert m to km)')
        points /= 1000
        print_verbose(f'{points=}')
        print_verbose(f'{points.shape=}')

    if points.shape[1] != 2:
        raise ValueError("Must give points of 2 dimensions")

    plot_points_2d(points, plot_ellipse, plot_evolute, limit_extents, dpi, km)

if __name__ == '__main__':
    sys.exit(main())
