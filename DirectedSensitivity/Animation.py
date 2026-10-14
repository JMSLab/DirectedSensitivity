"""An animation of the path from one certificate to the other."""
import math

import numpy as np

INK, MUTED = '#0b0b0b', '#52514e'
LOWER, UPPER = '#2a78d6', '#eb6834'
FRAMES_PER_SECOND = 12
GIF_DPI = 100
GIF_COLORS = 64


def DeformationPath(problem, points, curve=None):
    """The object, the fit, and the primitive at each point of a path through the baseline.

    ``points`` is the path the caller has chosen, running from the lower certificate
    through the baseline estimate to the upper one. Every point must be feasible, and one
    of them must be the baseline. ``curve`` maps a point to the values the animation
    plots; without it there is nothing to animate.
    """
    points = [np.atleast_1d(np.asarray(point, dtype=float)).ravel() for point in points]
    baseline = _Baseline(problem, points)
    frames = [_Frame(problem, point, position, curve)
              for point, position in zip(points, _Positions(points, baseline))]
    infeasible = [index for index, frame in enumerate(frames) if not frame['feasible']]
    if infeasible:
        raise ValueError(f"the path is not feasible at {len(infeasible)} of "
                         f"{len(frames)} points")
    return {'budget': _Number(problem.criterion_at_estimate),
            'baseline': baseline,
            'frames': frames}


def WriteAnimation(path, deformation, grid, labels=None, markers=None, title=None,
                   figure=None, annotate=None):
    """Write the path as a loop that runs out to one certificate and back.

    The loop closes on itself, so it reads as a deformation rather than a reel that
    snaps back at the end. ``path`` names a GIF.

    With ``figure`` the baseline and the two certificates are also written as a static
    figure, in the format its suffix names. ``annotate`` is called with that figure's
    axes before it is saved, and does not reach the loop: it marks the certificates,
    which is not what the loop shows.
    """
    frames = deformation['frames']
    grid = [_Number(value) for value in np.ravel(grid)]
    _CheckCurves(frames, grid)
    payload = {
        'title': title,
        'labels': dict({'x': 'x', 'y': 'primitive', 'object': 'c', 'criterion': r'$\hat{q}$',
                        'lower': 'lower certificate', 'upper': 'upper certificate'},
                       **(labels or {})),
        'grid': grid,
        'markers': _Markers(markers),
        'baseline': deformation['baseline'],
        'frames': frames,
    }
    path = _WriteGif(path, payload)
    if figure is not None:
        _WriteFigure(figure, payload, annotate)
    return path


def _Baseline(problem, points, tolerance=1e-9):
    """Where the baseline estimate sits on the path the loop runs along."""
    estimate = np.atleast_1d(np.asarray(problem.estimate, dtype=float)).ravel()
    distances = [float(np.linalg.norm(point - estimate)) for point in points]
    closest = int(np.argmin(distances))
    if distances[closest] > tolerance:
        raise ValueError("the path must pass through the baseline estimate")
    return closest


def _Positions(points, baseline):
    """Each point's signed share of the way to its end of the path, from -1 to 1."""
    legs = [float(np.linalg.norm(later - earlier))
            for earlier, later in zip(points, points[1:])]
    travelled = np.concatenate([[0.0], np.cumsum(legs)]) - sum(legs[:baseline])
    below, above = -travelled[0] or 1.0, travelled[-1] or 1.0
    return [value / (below if value < 0 else above) for value in travelled]


def _Frame(problem, point, position, curve):
    value = float(problem.target(point))
    slack = float(problem.Slack(point))
    return {'weight': _Number(position),
            'c': _Number(value),
            'criterion': _Number(problem.criterion_at_estimate - slack),
            'deviation': _Number(value - problem.target_at_estimate),
            'slack': _Number(slack),
            'feasible': bool(problem.IsFeasible(point)),
            'point': [_Number(value) for value in point],
            'curve': None if curve is None else [_Number(value)
                                                 for value in np.ravel(curve(point))]}


def _CheckCurves(frames, grid):
    widths = {len(frame['curve']) for frame in frames if frame['curve'] is not None}
    if not widths:
        raise ValueError("the animation needs the frames to carry curves")
    if widths != {len(grid)}:
        raise ValueError(f"the curves do not match the grid: {sorted(widths)} against "
                         f"{len(grid)}")


def _Markers(markers):
    if markers is None:
        return None
    horizontal = [_Number(value) for value in np.ravel(markers['x'])]
    vertical = [_Number(value) for value in np.ravel(markers['y'])]
    if len(horizontal) != len(vertical):
        raise ValueError("markers['x'] and markers['y'] must have the same length")
    return {'x': horizontal, 'y': vertical, 'label': markers.get('label')}


def _WriteFigure(path, payload, annotate):
    """The baseline and the two certificates, as a static figure."""
    chart = _Figure(figsize=(7.2, 4.6), constrained_layout=True)
    axes = chart.subplots()
    _Backdrop(axes, payload)
    if annotate is not None:
        annotate(axes)
    chart.savefig(path, dpi=200)
    return path


def _WriteGif(path, payload, fps=FRAMES_PER_SECOND):
    """The frames, as a loop that runs out to one certificate and back."""
    labels, frames = payload['labels'], payload['frames']
    chart = _Figure(figsize=(7.2, 4.6), dpi=GIF_DPI, constrained_layout=True)
    axes = chart.subplots()
    _Backdrop(axes, payload, alpha=0.4)
    if payload['title']:
        axes.set_title(payload['title'], fontsize=11, loc='left')
    axes.set_xlim(*_Span([payload['grid']]))
    axes.set_ylim(*_Span([frame['curve'] for frame in frames]
                         + ([payload['markers']['y']] if payload['markers'] else [])))
    moving, = axes.plot([], [], linewidth=2.4, zorder=6)
    readout = axes.text(0.98, 0.96, '', transform=axes.transAxes, ha='right', va='top',
                        fontsize=12, color=INK)

    def Render(index):
        frame = frames[index]
        moving.set_data(payload['grid'], _Curve(frame['curve']))
        moving.set_color(_Side(index, payload))
        readout.set_text(f"{labels['object']} = {frame['c']:.4g}\n"
                         f"{labels['criterion']} = {frame['criterion']:.4g}")
        chart.canvas.draw()
        return np.asarray(chart.canvas.buffer_rgba())[..., :3].copy()

    images = [Render(index) for index in range(len(frames))]
    out_and_back = list(range(len(frames))) + list(range(len(frames) - 2, 0, -1))
    return _SaveGif(path, [images[index] for index in out_and_back], fps)


def _SaveGif(path, images, fps, colors=GIF_COLORS):
    """One palette for every frame, so pixels unchanged since the last frame are stored
    as transparent."""
    from PIL import Image

    sample = np.concatenate(images[::max(1, len(images) // 16)], axis=0)
    palette = Image.fromarray(sample).quantize(colors, method=Image.Quantize.MAXCOVERAGE,
                                               dither=Image.Dither.NONE)
    frames = [Image.fromarray(image).quantize(palette=palette, dither=Image.Dither.NONE)
              for image in images]
    frames[0].save(path, save_all=True, append_images=frames[1:], loop=0,
                   duration=round(1000 / fps), optimize=True)
    return path


def _Backdrop(axes, payload, alpha=1.0):
    """The baseline and the two certificates, under whatever else the axes carry."""
    labels, frames = payload['labels'], payload['frames']
    for frame, colour, name, style in (
            (frames[payload['baseline']], INK, 'baseline', (0, (5, 2))),
            (frames[0], LOWER, labels['lower'], 'solid'),
            (frames[-1], UPPER, labels['upper'], 'solid')):
        axes.plot(payload['grid'], _Curve(frame['curve']), color=colour, linewidth=1.8,
                  linestyle=style, alpha=alpha,
                  label=f"{name}: {labels['object']} = {frame['c']:.6g}")
    if payload['markers']:
        axes.scatter(payload['markers']['x'], payload['markers']['y'], s=28,
                     facecolor='none', edgecolor=MUTED, zorder=5,
                     label=payload['markers']['label'])
    axes.set_xlabel(labels['x'])
    axes.set_ylabel(labels['y'])
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.legend(frameon=False, fontsize=9, loc='best')
    return axes


def _Side(index, payload):
    """The colour of the stop: ink at the baseline, one hue for each direction."""
    if index == payload['baseline']:
        return INK
    return LOWER if payload['frames'][index]['weight'] < 0 else UPPER


def _Span(series, pad=0.06):
    """A padded range covering every finite value, so the axes hold still."""
    values = [value for values in series for value in values
              if value is not None and math.isfinite(value)]
    low, high = min(values), max(values)
    room = (high - low or abs(high) or 1.0) * pad
    return low - room, high + room


def _Curve(values):
    """Curve values as floats, with anything non-finite left as a gap."""
    return np.array([np.nan if value is None else value for value in values],
                    dtype=float)


def _Figure(**options):
    """A figure drawn off-screen, leaving pyplot and its backend alone."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    chart = Figure(**options)
    FigureCanvasAgg(chart)
    return chart


def _Number(value):
    """A float, or None for anything not finite."""
    value = float(value)
    return value if math.isfinite(value) else None
