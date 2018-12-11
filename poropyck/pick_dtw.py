"""Dynamic time warping"""
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # pylint: disable=unused-import
from scipy.signal import hilbert
from scipy.stats import norm
import uncertainties
from dtw import dtw  # pylint: disable=no-name-in-module


# colours
HIGHLIGHT_COLOR = 'thistle'
ENVELOPE_COLOR = 'lightpink'
ENVELOPE_ALPHA = 0.4
PATH_COLOR = 'black'
SKIP_ROWS_IN_CSV = 21


class DTW:
    """compare using dynamic time warping"""

    def __init__(self, template_path, query_path, length_data,
                 template_color='tan', template_start=None, template_end=None,
                 query_color='blue', query_start=None, query_end=None):
        self.fig = None
        self.ax = None
        self.template_path = template_path
        self.query_path = query_path
        self.template_color = template_color
        self.query_color = query_color
        length_mean = np.mean(length_data)
        length_std = np.std(length_data)
        self.length = (
            length_mean if length_std == 0.0
            else uncertainties.ufloat(length_mean, length_std)
        )
        self.template = Signal(
            np.loadtxt(template_path, delimiter=',',
                       skiprows=SKIP_ROWS_IN_CSV).T[:2],
            self.length,
            color=template_color,
            window_start=template_start,
            window_end=template_end
        )
        self.query = Signal(
            np.loadtxt(query_path, delimiter=',',
                       skiprows=SKIP_ROWS_IN_CSV).T[:2],
            self.length,
            color=query_color,
            window_start=query_start,
            window_end=query_end
        )
        self.indices1 = None
        self.indices1a = None
        self.indices1h = None
        self.indices2 = None
        self.indices2a = None
        self.indices2h = None
        self.summary_xlim = None
        self.summary_ylim = None

    def pick(self):
        """show plots and start the picking process"""
        self.fig = plt.figure()
        self.ax = {
            'template': self.fig.add_axes([0.03, 0.865, 0.9, 0.1]),
            'query': self.fig.add_axes([0.03, 0.7, 0.9, 0.1]),
            'dtw': self.fig.add_axes([0, 0.15, 0.2*1.65, 0.22*1.65], projection='3d'),
            'x': self.fig.add_axes([0.42, 0.48, 0.25, 0.1]),
            'y': self.fig.add_axes([0.35, 0.08, 0.07, 0.4]),
            'template_clicks': self.fig.add_axes([0.72, 0.38, 0.1, 0.2]),
            'template_velocity': self.fig.add_axes([0.84, 0.38, 0.1, 0.2]),
        }
        self.ax['summary'] = self.fig.add_axes(
            [0.42, 0.08, 0.25, 0.4],
            sharex=self.ax['x'],
            sharey=self.ax['y']
        )
        self.ax['query_clicks'] = self.fig.add_axes(
            [0.72, 0.08, 0.1, 0.2],
            sharex=self.ax['template_clicks']
        )
        self.ax['query_velocity'] = self.fig.add_axes(
            [0.84, 0.08, 0.1, 0.2],
            sharex=self.ax['template_velocity']
        )
        self.ax['template'].get_yaxis().set_visible(False)
        self.ax['query'].get_yaxis().set_visible(False)
        self.ax['x'].get_yaxis().set_visible(False)
        self.ax['y'].get_xaxis().set_visible(False)
        self.ax['template_clicks'].get_yaxis().set_visible(False)
        self.ax['query_clicks'].get_yaxis().set_visible(False)
        self.ax['template_velocity'].get_yaxis().set_visible(False)
        self.ax['query_velocity'].get_yaxis().set_visible(False)
        self.ax['summary'].get_yaxis().tick_right()
        self.ax['x'].get_xaxis().tick_top()
        self.ax['template'].set_title(
            'Template signal: select the area of interest')
        self.ax['query'].set_title(
            'Query signal: select the area of interest')

        self.fig.canvas.mpl_connect('button_press_event', self.onpress)
        self.fig.canvas.mpl_connect('button_release_event', self.onrelease)
        self.fig.canvas.mpl_connect('motion_notify_event', self.onmotion)
        self.fig.canvas.mpl_connect('pick_event', self.onpick)

        self.summary_xlim = self.ax['summary'].get_xlim()
        self.summary_ylim = self.ax['summary'].get_ylim()

        self.template.plot(self.ax['template'])
        self.query.plot(self.ax['query'])
        self.run_dtw()
        self.highlight_summary()
        self.plot_results()
        plt.show()
        return {
            'file': self.query_path,
            'window_start': self.query.pick_start,
            'window_end': self.query.pick_end,
            'distance': get_mean(self.length),
            'distance_error': get_std(self.length),
            'time': get_mean(self.query.time),
            'time_error': get_std(self.query.time),
            'velocity': get_mean(self.query.velocity),
            'velocity_error': get_std(self.query.velocity),
            'template': {
                'file': self.template_path,
                'window_start': self.template.pick_start,
                'window_end': self.template.pick_end,
                'time': get_mean(self.template.time),
                'time_error': get_std(self.template.time),
                'velocity': get_mean(self.template.velocity),
                'velocity_error': get_std(self.template.velocity)
            }
        }

    def onpress(self, event):
        """mouse button pressed"""
        if event.inaxes is self.ax['template'] or event.inaxes is self.ax['query']:
            if event.inaxes is self.ax['template']:
                self.template.onpress(event)
                self.template.plot(self.ax['template'])
            if event.inaxes is self.ax['query']:
                self.query.onpress(event)
                self.query.plot(self.ax['query'])
            if event.inaxes is self.ax['template']:
                self.template.pick_start = None
                self.template.pick_end = None
            if event.inaxes is self.ax['query']:
                self.query.pick_start = None
                self.query.pick_end = None
            self.clear_output_axes()
            self.fig.canvas.draw_idle()

    def onrelease(self, event):
        """mouse button released"""
        if self.template.pressed or self.query.pressed:
            if self.template.pressed:
                self.template.onrelease(event)
                self.template.plot(self.ax['template'])
            if self.query.pressed:
                self.query.onrelease(event)
                self.query.plot(self.ax['query'])
            self.run_dtw()
            self.highlight_summary()
            self.plot_results()
            self.fig.canvas.draw_idle()

    def onmotion(self, event):
        """mouse moves"""
        if event.inaxes is self.ax['template'] or event.inaxes is self.ax['query']:
            if self.template.pressed or self.query.pressed:
                if self.template.pressed:
                    self.template.move_line(event)
                    self.template.plot(self.ax['template'])
                if self.query.pressed:
                    self.query.move_line(event)
                    self.query.plot(self.ax['query'])
                self.fig.canvas.draw_idle()

    def onpick(self, event):
        """pick values from the active plot"""
        indices = event.ind
        xdata = event.artist.get_xdata()[indices]
        ydata = event.artist.get_ydata()[indices]
        mouse_x, mouse_y = event.mouseevent.xdata, event.mouseevent.ydata
        dists = np.sqrt((np.array(xdata)-mouse_x)**2 +
                        (np.array(ydata)-mouse_y)**2)
        xpoint = xdata[np.argmin(dists)]
        ypoint = ydata[np.argmin(dists)]

        if self.template.pick_start is None or self.query.pick_start is None:
            self.template.pick_start = ypoint
            self.query.pick_start = xpoint
        elif self.template.pick_end is None or self.query.pick_end is None:
            self.template.pick_end = ypoint
            self.query.pick_end = xpoint
        else:
            dist_to_start = np.sqrt(
                (xpoint - self.query.pick_start)**2 +
                (ypoint - self.template.pick_start)**2
            )
            dist_to_end = np.sqrt(
                (xpoint - self.query.pick_end)**2 +
                (ypoint - self.template.pick_end)**2
            )
            if dist_to_start <= dist_to_end:
                self.template.pick_start = ypoint
                self.query.pick_start = xpoint
            else:
                self.template.pick_end = ypoint
                self.query.pick_end = xpoint
        self.plot_summary(self.ax['x'], self.ax['y'], self.ax['summary'])
        self.highlight_summary()
        self.plot_results()
        self.fig.canvas.draw_idle()

    def onxzoom(self, axes):
        """summary axes is zoomed - x"""
        self.summary_xlim = axes.get_xlim()

    def onyzoom(self, axes):
        """summary axes is zoomed - y"""
        self.summary_ylim = axes.get_ylim()

    def run_dtw(self):
        """run dynamic time warping"""
        self.indices1, self.indices2 = dtw(
            self.query.picked_signal, self.template.picked_signal)[1:3]
        self.indices1h, self.indices2h = dtw(
            self.query.hilbert_abs(), self.template.hilbert_abs())[1:3]
        self.indices1a, self.indices2a = dtw(
            self.query.hilbert_angle(), self.template.hilbert_angle())[1:3]

        self.plot_dtw(self.ax['dtw'])
        self.plot_summary(self.ax['x'], self.ax['y'], self.ax['summary'])

    def plot_dtw(self, ax):
        """plot the 3D DTW data"""
        template_times = self.template.picked_times
        template_signal = self.template.picked_signal
        query_times = self.query.picked_times
        query_signal = self.query.picked_signal
        ax.clear()
        ax.set_title('Dynamic Time Warping Visualization', y=1.15)
        ax.plot(template_times, template_signal, zs=1, c=self.template.color)
        ax.plot(query_times, query_signal, zs=0, c=self.query.color)
        for i in np.arange(0, len(self.indices1) - 1, 20):
            x_start = np.take(template_times, self.indices2[i].astype(int) - 1)
            x_end = np.take(query_times, self.indices1[i].astype(int) - 1)
            y_start = np.take(
                template_signal, self.indices2[i].astype(int) - 1)
            y_end = np.take(query_signal, self.indices1[i].astype(int) - 1)
            ax.plot(
                [x_start, x_end], [y_start, y_end],
                '-', color=HIGHLIGHT_COLOR, lw=0.5, zs=[1, 0]
            )
        self.summary_xlim = (self.query.start, self.query.finish)
        self.summary_ylim = (self.template.start, self.template.finish)

    def plot_summary(self, x_ax, y_ax, summary_ax):
        """plot the time warping summary"""
        query_times = self.query.picked_times
        query_signal = self.query.picked_signal
        queryh = self.query.hilbert_abs()
        x_ax.clear()
        x_ax.set_title(
            ('Use predicted arrival shown or\n' +
             'pick a 95% confidence interval.\n' +
             'Close window when complete.'),
            y=1.25
        )
        x_ax.plot(query_times, query_signal, '-', c=self.query_color, lw=2)
        x_ax.fill_between(query_times, queryh, -queryh,
                          color=ENVELOPE_COLOR, alpha=ENVELOPE_ALPHA)
        try:
            x_ax.set_ylim(1.1 * np.min(query_signal),
                          1.1 * np.max(queryh))
        except TypeError:
            x_ax.set_ylim(-1, 1)
        x_ax.grid(True, axis='x', color='lightgrey')
        x_ax.tick_params(axis='x', which='major')

        template_times = self.template.picked_times
        template_signal = self.template.picked_signal
        templateh = self.template.hilbert_abs()
        y_ax.clear()
        y_ax.plot(template_signal, template_times, c=self.template_color, lw=2)
        y_ax.fill_betweenx(template_times, templateh, -templateh,
                           color=ENVELOPE_COLOR, alpha=ENVELOPE_ALPHA)
        try:
            y_ax.set_xlim(1.1 * np.min(template_signal),
                          1.1 * np.max(templateh))
        except TypeError:
            y_ax.set_xlim(-1, 1)
        y_ax.grid(True, axis='y', color='lightgrey')
        y_ax.tick_params(axis='y', which='major')
        y_ax.invert_xaxis()

        end = len(self.indices1)
        idxto = np.take(template_times, self.indices2[:end].astype(int) - 1)
        idxqo = np.take(query_times, self.indices1.astype(int) - 1)

        end = len(self.indices1h)
        idxtoh = np.take(template_times, self.indices2h[:end].astype(int) - 1)
        idxqoh = np.take(query_times, self.indices1h.astype(int) - 1)

        summary_ax.clear()
        summary_ax.fill_between(idxqoh, idxtoh, idxqoh,
                                color=ENVELOPE_COLOR, alpha=ENVELOPE_ALPHA)
        summary_ax.fill_betweenx(idxqoh, idxtoh, idxqoh,
                                 color=ENVELOPE_COLOR, alpha=ENVELOPE_ALPHA)
        summary_ax.plot(idxqo, idxto, c=self.template_color, picker=5, lw=2)
        summary_ax.plot(idxqo, idxto, c=self.query_color, alpha=0.25, lw=2)

        times = np.linspace(
            0, np.max([self.template.times, self.query.times]), 10)
        summary_ax.plot(times, times, color=HIGHLIGHT_COLOR, lw=1)
        summary_ax.tick_params(axis='both', which='major')
        summary_ax.grid(color='lightgrey')
        summary_ax.set_xlim(self.summary_xlim)
        summary_ax.set_ylim(self.summary_ylim)
        summary_ax.callbacks.connect('xlim_changed', self.onxzoom)
        summary_ax.callbacks.connect('ylim_changed', self.onyzoom)

        x_ax.set_xlim(
            min(summary_ax.get_xlim()[0], summary_ax.get_ylim()[0]),
            max(summary_ax.get_xlim()[1], summary_ax.get_ylim()[1])
        )
        y_ax.set_ylim(
            min(summary_ax.get_xlim()[0], summary_ax.get_ylim()[0]),
            max(summary_ax.get_xlim()[1], summary_ax.get_ylim()[1])
        )

    def highlight_summary(self):
        """highlight summary plot after points are picked"""
        ax = self.ax['summary']
        ax_x = self.ax['x']
        ax_y = self.ax['y']

        mean, std = self.query.time_picks()
        min_ = mean - 2 * std
        max_ = mean + 2 * std
        ax.axvspan(min_, max_, alpha=0.4, color=self.query.color)
        ax.axvline(mean, linewidth=2, color=self.query.color)
        ax_x.axvspan(min_, max_, alpha=0.4, color=self.query.color)
        ax_x.axvline(mean, linewidth=2, color=self.query.color)

        mean, std = self.template.time_picks()
        min_ = mean - 2 * std
        max_ = mean + 2 * std
        ax.axhspan(min_, max_, alpha=0.4, color=self.template.color)
        ax.axhline(mean, linewidth=2, color=self.template.color)
        ax_y.axhspan(min_, max_, alpha=0.4, color=self.template.color)
        ax_y.axhline(mean, linewidth=2, color=self.template.color)

    def plot_results(self):
        """plot the distributions"""
        self.template.plot_time(self.ax['template_clicks'])
        self.template.plot_velocity(self.ax['template_velocity'])
        self.query.plot_time(self.ax['query_clicks'])
        self.query.plot_velocity(self.ax['query_velocity'])

    def clear_output_axes(self):
        """clear all output data"""
        self.ax['dtw'].clear()
        self.ax['summary'].clear()
        self.ax['x'].clear()
        self.ax['y'].clear()
        self.ax['template_clicks'].clear()
        self.ax['query_clicks'].clear()
        self.ax['template_velocity'].clear()
        self.ax['query_velocity'].clear()


def aic(signal):
    """estimate pick point"""
    length = len(signal)
    output = [0]
    with np.errstate(divide='ignore'):
        for k in range(1, length - 1):
            val = (k * np.log(np.var(signal[0:k])) +
                   (length-k-1) * np.log(np.var(signal[k+1:length])))
            if val == -np.inf:
                val = 0
            output.append(val)
    output.append(0)
    return np.argmin(output)


class Signal:
    """one signal to be compared"""

    def __init__(self, data, length, color='blue', window_start=None, window_end=None):
        self.times, self.signal = data
        self.length = length
        self.color = color
        self.start = window_start
        self.finish = window_end
        width = len(self.signal) // 80
        predicted = aic(self.signal)
        if not self.start:
            self.start = self.times[
                max(0, predicted - width)
            ]
        if not self.finish:
            self.finish = self.times[
                min(len(self.signal) - 1, predicted + width)
            ]
        self.pressed = False
        self.pick_start = self.times[predicted]
        self.pick_end = None
        self.picked_times, self.picked_signal = self.get_picked_data()
        self.velocity = None
        self.time = None

    def onpress(self, event):
        """mouse button pressed"""
        self.pressed = True
        self.move_line(event)

    def onrelease(self, event):
        """mouse button released"""
        self.pressed = False
        self.move_line(event)
        self.picked_times, self.picked_signal = self.get_picked_data()
        predicted = aic(self.picked_signal)
        self.pick_start = self.picked_times[predicted]
        self.pick_end = None

    def move_line(self, event):
        """move the nearest line"""
        click = event.xdata
        if abs(click - self.start) < abs(click - self.finish):
            self.start = click
        else:
            self.finish = click

    def get_picked_data(self):
        """return picked data"""
        pick = np.logical_and(self.start <= self.times,
                              self.times <= self.finish)
        times = np.extract(pick, self.times)
        signal = np.extract(pick, self.signal)
        absmax = np.max(np.abs(signal))
        return times, signal / absmax

    def plot(self, ax):
        """plot the signal over time"""
        title = ax.get_title()
        ax.clear()
        ax.set_title(title)
        ax.axvspan(self.start, self.finish, alpha=0.4, color=self.color)
        ax.plot(self.times, self.signal, color=self.color)

    def plot_time(self, ax):
        """plot time distribution"""
        ax.clear()
        ax.set_title('time picks')
        mean, std = self.time_picks()
        ax.set_title('Time\n{:5g}'.format(mean))
        if std != 0.0:
            range_ = np.linspace(mean - 2 * std, mean + 2 * std, 50)
            norm_ = norm.pdf(range_, mean, std)
            ax.plot(range_, norm_, '--', c=self.color, lw=2)
            ax.set_title('Time\n{:5g}±{:5g}'.format(mean, std))
        ax.set_xlabel(r'$\mu$s')

    def plot_velocity(self, ax):
        """plot velocity distribution"""
        time_mean, time_std = self.time_picks()
        self.time = (
            time_mean if time_std == 0.0
            else uncertainties.ufloat(time_mean, time_std)
        )
        self.velocity = (self.length / self.time) * 1e4
        ax.clear()
        plt.sca(ax)
        if isinstance(self.velocity, float):
            ax.set_title('Velocity\n{:5g}'.format(self.velocity))
        else:
            x = np.linspace(
                norm.ppf(0.01, self.velocity.n, self.velocity.s),
                norm.ppf(0.99, self.velocity.n, self.velocity.s),
                1000
            )
            ax.plot(x, norm.pdf(x, self.velocity.n, self.velocity.s),
                    color=self.color, lw=2, ls='dashed')
            ax.set_title('Velocity\n{:5g}±{:5g}'.format(
                self.velocity.n, self.velocity.s))
        ax.set_xlabel('m/s')

    def hilbert_angle(self):
        """get the angle of the hilbert transform"""
        return np.angle(hilbert(self.picked_signal)) / np.pi

    def hilbert_abs(self):
        """get the absolute value of the hilbert transform"""
        return np.abs(hilbert(self.picked_signal))

    def time_picks(self):
        """return picked time data"""
        if self.pick_start:
            if self.pick_end:
                two_std = abs(self.pick_end - self.pick_start) / 2
                mean = min(self.pick_start, self.pick_end) + two_std
                return mean, two_std / 2
            return self.pick_start, 0.0
        return 0.0, 0.0


def get_mean(x):
    if isinstance(x, float):
        return x
    return x.n


def get_std(x):
    if isinstance(x, float):
        return 0.0
    return x.s
