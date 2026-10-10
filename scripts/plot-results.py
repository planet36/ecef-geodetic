# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: MPL-2.0

"""Plot speed against accuracy from an acc-speed.filtered.csv file."""

# The module name comes from the file name, which has a hyphen.
# pylint: disable=invalid-name

__author__ = 'Steven Ward'
__license__ = 'MPL-2.0'
__version__ = '2026-10-10'

import argparse
import csv
import math
import sys
from typing import Literal

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.transforms import Bbox

# Columns of the filtered CSV
NAME_COL = 'Name'
MEAN_ERR_COL = 'Mean dist. error (nm)'
TIME_COL = 'Time/call (ns)'

# process-results.bash keeps only the algorithms whose mean error is under this many nm.
MAX_MEAN_ERR_NM = 10

# Colors of the dark theme
SURFACE = '#1a1a19'
GRID = '#3a3a37'
TEXT_PRIMARY = '#ffffff'
TEXT_SECONDARY = '#c3c2b7'
LEADER = '#8a8983'
SERIES = '#3987e5'

# Each axis ends on the first whole step that leaves at least this fraction of the axis past
# its largest value, so that no marker is cut off at the edge.
EDGE_MARGIN = 0.03

# The marker is this many points across.  It has a 2-point ring in the surface color.
MARKER_SIZE_PT = 10

# A direction from a point, with the horizontal and vertical text alignment of a label there
HAlign = Literal['left', 'center', 'right']
VAlign = Literal['bottom', 'center', 'top']
Direction = tuple[tuple[int, int], HAlign, VAlign]

# The directions a label may sit in from its point, in order of preference, with the text
# alignment that keeps the label on that side
LABEL_DIRECTIONS: tuple[Direction, ...] = (
    ((0, 1), 'center', 'bottom'),
    ((1, 0), 'left', 'center'),
    ((-1, 0), 'right', 'center'),
    ((0, -1), 'center', 'top'),
    ((1, 1), 'left', 'bottom'),
    ((-1, 1), 'right', 'bottom'),
    ((1, -1), 'left', 'top'),
    ((-1, -1), 'right', 'top'),
)

# The distances a label may sit at from its point, in points, in order of preference.  A
# label beyond the first distance gets a leader line back to its point.
LABEL_DISTANCES_PT = (8, 20, 34, 50, 70, 95)

# A point is more crowded the more points lie within this many typographic points of it.  The
# labels of the most crowded points are placed first.
CROWDING_RADIUS_PT = 150

# Fixed placements for labels that the search leaves next to another point, with the direction,
# alignment, and distance in the forms above.  These labels are placed before the others.
LABEL_OVERRIDES: dict[str, tuple[tuple[int, int], HAlign, VAlign, int]] = {
    'Bowring 1976': ((1, 0), 'left', 'center', 8),
    'Halley': ((-1, 0), 'right', 'center', 8),
    'Heikkinen 1982': ((-1, 0), 'right', 'center', 8),
    'Householder': ((-1, -1), 'right', 'top', 8),
    'Lin-Wang 1995': ((-1, 1), 'right', 'bottom', 8),
    'Lin-Wang 1995 (c.h.)': ((0, 1), 'center', 'bottom', 20),
    'Schroder': ((1, 1), 'left', 'bottom', 8),
    'Shu 2010': ((0, -1), 'center', 'top', 8),
}


def load_results(path: str) -> tuple[list[str], list[float], list[float]]:
    """Return the names, mean errors (nm), and speeds (M conversions/s) in a CSV file."""
    with open(path, newline='', encoding='utf-8') as f:
        # The first two lines name the input files, and the third is the header.
        f.readline()
        f.readline()
        reader = csv.DictReader(f)
        rows = list(reader)

    fieldnames = reader.fieldnames or []
    for col in (NAME_COL, MEAN_ERR_COL, TIME_COL):
        if col not in fieldnames:
            raise ValueError(f"no {col!r} column (have: {', '.join(fieldnames)})")

    if not rows:
        raise ValueError('no data rows')

    mean_errs_nm = []
    times_ns = []

    for r in rows:
        # DictReader fills the fields missing from a short row with None.
        if None in r.values():
            raise ValueError(f"{r[NAME_COL]!r} has too few fields")
        try:
            mean_err_nm = float(r[MEAN_ERR_COL])
            time_ns = float(r[TIME_COL])
        except ValueError as e:
            raise ValueError(f"{r[NAME_COL]!r}: {e}") from e
        if not 0 <= mean_err_nm < MAX_MEAN_ERR_NM:
            raise ValueError(f"{r[NAME_COL]!r} has a mean error of {mean_err_nm} nm, which is "
                             f"not under {MAX_MEAN_ERR_NM} nm (is this the filtered file?)")
        if not 0 < time_ns < math.inf:
            raise ValueError(f"{r[NAME_COL]!r} has a time of {time_ns} ns, which is not "
                             "positive and finite")
        mean_errs_nm.append(mean_err_nm)
        times_ns.append(time_ns)

    # Every iterative algorithm does 2 iterations, as the README says.
    names = [r[NAME_COL].removesuffix(' (x2)') for r in rows]
    speeds = [1000 / time_ns for time_ns in times_ns]

    return names, mean_errs_nm, speeds


# pylint: disable-next=too-many-locals
def place_labels(ax: Axes, names: list[str], xs: list[float], ys: list[float]) -> None:
    """Label every point where the label overlaps no other label or point, if it can.

    Each label takes the first placement, in order of preference, that overlaps nothing placed
    so far.  When every placement overlaps something, it takes the one with the fewest
    overlaps.  A label in LABEL_OVERRIDES takes its fixed placement instead, unless that
    placement leaves the axes.
    """
    px_per_pt = ax.figure.dpi / 72
    marker_r_px = MARKER_SIZE_PT / 2 * px_per_pt
    axes_bbox = ax.get_window_extent()

    points = [tuple(ax.transData.transform((x, y))) for (x, y) in zip(xs, ys)]
    obstacles = [Bbox.from_extents(px - marker_r_px, py - marker_r_px,
                                   px + marker_r_px, py + marker_r_px) for (px, py) in points]

    # Label the fixed placements first, then the most crowded points, while the most room is
    # left around them.
    def label_order(i: int) -> tuple[bool, int]:
        return (names[i] not in LABEL_OVERRIDES,
                -sum(math.dist(points[i], q) < CROWDING_RADIUS_PT * px_per_pt for q in points))

    for i in sorted(range(len(names)), key=label_order):
        (px, py) = points[i]
        # https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.annotate.html
        text = ax.annotate(names[i], (xs[i], ys[i]), xytext=(0, 0), textcoords='offset points',
                           color=TEXT_SECONDARY, fontsize=11, zorder=4)
        extent = text.get_window_extent()
        (w, h) = (extent.width, extent.height)

        search = [(dist_pt, d) for dist_pt in LABEL_DISTANCES_PT for d in LABEL_DIRECTIONS]
        searches: tuple[list[tuple[int, Direction]], ...]
        if names[i] in LABEL_OVERRIDES:
            (direction, ha, va, dist_pt) = LABEL_OVERRIDES[names[i]]
            # An override that leaves the axes falls back to the search.
            searches = ([(dist_pt, (direction, ha, va))], search)
        else:
            searches = (search,)

        placements = []
        for candidates in searches:
            for (dist_pt, ((dx, dy), ha, va)) in candidates:
                scale = dist_pt / math.hypot(dx, dy)
                offset = (dx * scale, dy * scale)
                (anchor_x, anchor_y) = (px + offset[0] * px_per_pt, py + offset[1] * px_per_pt)
                x0 = {'left': anchor_x, 'center': anchor_x - w / 2, 'right': anchor_x - w}[ha]
                y0 = {'bottom': anchor_y, 'center': anchor_y - h / 2, 'top': anchor_y - h}[va]
                bbox = Bbox.from_bounds(x0, y0, w, h)
                if not (axes_bbox.contains(bbox.x0, bbox.y0) and
                        axes_bbox.contains(bbox.x1, bbox.y1)):
                    continue
                overlaps = sum(bbox.overlaps(o) for o in obstacles)
                placements.append((overlaps, dist_pt, offset, ha, va, bbox))
                if overlaps == 0:
                    break
            if placements:
                break

        if not placements:
            raise ValueError(f"no placement of the label {names[i]!r} fits inside the axes")

        # min returns the first of equal placements, which is the one most preferred.
        (_, dist_pt, offset, ha, va, bbox) = min(placements, key=lambda c: c[0])

        # https://matplotlib.org/stable/api/text_api.html#matplotlib.text.Annotation
        text.xyann = offset
        text.set_horizontalalignment(ha)
        text.set_verticalalignment(va)
        obstacles.append(bbox)

        if dist_pt > LABEL_DISTANCES_PT[0]:
            # An annotation without text draws its arrow from exactly the offset to the point.
            ax.annotate('', (xs[i], ys[i]), xytext=offset, textcoords='offset points',
                        arrowprops={'arrowstyle': '-', 'color': LEADER, 'linewidth': 0.8,
                                    'shrinkA': 0, 'shrinkB': 0}, zorder=1)


def axis_limit(largest: float, step: int) -> int:
    """Return the first multiple of step that leaves EDGE_MARGIN of the axis past largest."""
    return step * (math.floor(largest / (1 - EDGE_MARGIN) / step) + 1)


def plot_results(names: list[str], mean_errs_nm: list[float], speeds: list[float]) -> Figure:
    """Return a scatter plot of speed against mean error."""
    fig, ax = plt.subplots(figsize=(16, 11.5))

    fig.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    ax.set_title('ECEF-to-Geodetic\nSpeed vs. accuracy', color=TEXT_PRIMARY, fontsize=18,
                 fontweight='bold', pad=20)

    # https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set_xlabel.html
    ax.set_xlabel('Mean distance error (nm)', color=TEXT_SECONDARY, fontsize=13, labelpad=10)
    # https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set_ylabel.html
    ax.set_ylabel('Million conversions per second', color=TEXT_SECONDARY, fontsize=13,
                  labelpad=10)

    ax.set_xlim(0, axis_limit(max(mean_errs_nm), 1))
    ax.set_ylim(0, axis_limit(max(speeds), 10))

    ax.grid(color=GRID, linewidth=1)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_color(GRID)
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=11, length=0, pad=8)

    # https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.scatter.html
    ax.scatter(mean_errs_nm, speeds, s=MARKER_SIZE_PT**2, color=SERIES, edgecolors=SURFACE,
               linewidths=2, zorder=3)

    # The labels are placed in display coordinates, so the layout must be final first.
    fig.tight_layout()
    place_labels(ax, names, mean_errs_nm, speeds)

    return fig


def main() -> None:
    """Parse the options, then show or save the plot."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results_file',
                        help='an acc-speed.filtered.csv file from process-results.bash')
    parser.add_argument('-o', '--output',
                        help='save the plot to this file (for example, a .png) instead of '
                             'showing it')
    args = parser.parse_args()

    try:
        results = load_results(args.results_file)
    except OSError as e:
        parser.error(f"{args.results_file}: {e.strerror}")
    except ValueError as e:
        parser.error(f"{args.results_file}: {e}")

    for name in sorted(LABEL_OVERRIDES.keys() - set(results[0])):
        print(f"{parser.prog}: warning: no point is named {name!r}, so its label override is "
              "unused", file=sys.stderr)

    if args.output:
        # Saving a file needs no window.  An interactive backend would still load its GUI
        # toolkit, which can print graphics warnings.
        plt.switch_backend('agg')

    fig = plot_results(*results)
    if args.output:
        fig.savefig(args.output, dpi=200)
        print(f"Created file:\n{args.output}")
    else:
        plt.show()


if __name__ == '__main__':
    main()
