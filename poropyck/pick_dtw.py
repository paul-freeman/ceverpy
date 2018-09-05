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


class DTW:
    """compare using dynamic time warping"""

    def __init__(self, template_path, query_path, length_data, template_color='tan', query_color='blue'):
        self.fig = None
        self.ax = None
        self.template_color = template_color
        self.query_color = query_color
        length_mean = np.mean(length_data)
        length_std = np.std(length_data)
        self.length = (
            length_mean if length_std == 0.0
            else uncertainties.ufloat(length_mean, length_std)
        )
        self.template = Signal(
            np.loadtxt(template_path, delimiter=',', skiprows=21).T[:2],
            self.length,
            color=template_color
        )
        self.query = Signal(
            np.loadtxt(query_path, delimiter=',', skiprows=21).T[:2],
            self.length,
            color=query_color
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
            'template': self.fig.add_axes([0.030, 0.865, 0.90, 0.100]),
            'query': self.fig.add_axes([0.030, 0.700, 0.90, 0.100]),
            'dtw': self.fig.add_axes([0, 0.15, 0.200*1.65, 0.220*1.65], projection='3d'),
            'summary': self.fig.add_axes([0.100+0.32, 0.08, 0.25, 0.4]),
            'template_clicks': self.fig.add_axes([0.400+0.32, 0.38, 0.1, 0.200]),
            'template_velocity': self.fig.add_axes([0.520+0.32, 0.38, 0.1, 0.200]),
            'query_clicks': self.fig.add_axes([0.400+0.32, 0.080, 0.1, 0.200]),
            'query_velocity': self.fig.add_axes([0.520+0.32, 0.080, 0.1, 0.200]),
        }
        self.ax['x'] = self.fig.add_axes(
            [0.100+0.32, 0.48, 0.25, 0.1],
            sharex=self.ax['summary']
        )
        self.ax['y'] = self.fig.add_axes(
            [0.030+0.32, 0.08, 0.07, 0.4],
            sharey=self.ax['summary']
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
        plt.show()
        template_data = {
            'distance': self.length,
            'time': self.template.time,
            'velocity': self.template.velocity
        }
        query_data = {
            'distance': self.length,
            'time': self.query.time,
            'velocity': self.query.velocity
        }
        return template_data, query_data

    def onpress(self, event):
        """mouse button pressed"""
        if event.inaxes is self.ax['template'] or event.inaxes is self.ax['query']:
            if event.inaxes is self.ax['template']:
                self.template.onpress(event)
                self.template.plot(self.ax['template'])
            if event.inaxes is self.ax['query']:
                self.query.onpress(event)
                self.query.plot(self.ax['query'])
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

        self.template.picks.append(ypoint)
        self.query.picks.append(xpoint)
        self.plot_summary(self.ax['x'], self.ax['y'], self.ax['summary'])
        self.highlight_summary()
        self.plot_monte_carlo()
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
        x_ax.set_title('Select points of interest', y=1.25)
        x_ax.plot(query_times, query_signal, '-', c=self.query_color, lw=2)
        x_ax.fill_between(query_times, queryh, -queryh,
                          color=ENVELOPE_COLOR, alpha=ENVELOPE_ALPHA)
        try:
            x_ax.set_ylim(1.1 * np.min(query_signal),
                          1.1 * np.max(queryh))
        except TypeError:
            x_ax.set_ylim(-1, 1)
        x_ax.grid(True, axis='x', color='lightgrey')
        y_ax.tick_params(axis='x', which='major')

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
        summary_ax.axis('equal')
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

        x_ax.set_xlim(summary_ax.get_xlim())
        y_ax.set_ylim(summary_ax.get_ylim())

    def highlight_summary(self):
        """highlight summary plot after points are picked"""
        ax = self.ax['summary']
        ax_x = self.ax['x']
        ax_y = self.ax['y']

        min_, max_, mean, _ = self.query.time_picks()
        ax.axvspan(min_, max_, alpha=0.4, color=self.query.color)
        ax.axvline(mean, linewidth=2, color=self.query.color)
        ax_x.axvspan(min_, max_, alpha=0.4, color=self.query.color)
        ax_x.axvline(mean, linewidth=2, color=self.query.color)

        min_, max_, mean, _ = self.template.time_picks()
        ax.axhspan(min_, max_, alpha=0.4, color=self.template.color)
        ax.axhline(mean, linewidth=2, color=self.template.color)
        ax_y.axhspan(min_, max_, alpha=0.4, color=self.template.color)
        ax_y.axhline(mean, linewidth=2, color=self.template.color)

    def plot_monte_carlo(self):
        """plot the Monte Carlo distributions"""
        self.template.plot_clicks(self.ax['template_clicks'])
        self.template.plot_velocity(self.ax['template_velocity'])
        self.query.plot_clicks(self.ax['query_clicks'])
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


class Signal:
    """one signal to be compared"""

    def __init__(self, data, length, color='blue'):
        secs, self.signal = data
        self.times = secs * 1e6
        self.length = length
        self.color = color
        self.start = self.times[len(self.times) // 2]
        self.finish = self.times[len(self.times) // 2 + 1]
        self.pressed = False
        self.picks = []
        self.picked_times, self.picked_signal = self.get_picked_data()
        self.velocity = None
        self.time = None

    def onpress(self, event):
        """mouse button pressed"""
        self.pressed = True
        self.picks = []
        self.move_line(event)

    def onrelease(self, event):
        """mouse button released"""
        self.pressed = False
        self.move_line(event)
        self.picked_times, self.picked_signal = self.get_picked_data()

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

    def plot_clicks(self, ax):
        """plot click histogram"""
        ax.clear()
        ax.set_title('time picks')
        min_, max_, mean, std = self.time_picks()
        range_ = np.linspace(min_ - 2 * std, max_ + 2 * std, 50)
        norm_ = norm.pdf(range_, mean, std)
        ax.hist(np.array(self.picks), normed=True, color=self.color, alpha=0.6)
        ax.plot(range_, norm_, '--', c=self.color, lw=2)
        ax.set_title('Time ({} clicks)\n{:5g}±{:5g}'.format(
            len(self.picks), mean, std))
        ax.set_xlabel(r'$\mu$s')

    def plot_velocity(self, ax):
        """plot velocity distribution"""
        time_mean, time_std = self.time_picks()[2:]
        distance = self.length
        self.time = time_mean if time_std == 0.0 else uncertainties.ufloat(
            time_mean, time_std)
        self.velocity = (self.length / self.time) * 1e4
        ax.clear()
        plt.sca(ax)
        if isinstance(self.velocity, uncertainties.core.AffineScalarFunc):
            x = np.linspace(
                norm.ppf(0.01, self.velocity.n, self.velocity.s),
                norm.ppf(0.99, self.velocity.n, self.velocity.s),
                1000
            )
            ax.plot(x, norm.pdf(x, self.velocity.n, self.velocity.s),
                    color=self.color, lw=2, ls='dashed')
            ax.set_title('Velocity\n{:5g}±{:5g}'.format(
                self.velocity.n, self.velocity.s))
        else:
            self.velocity.plot(color=self.color, lw=2, ls='dashed')
            self.velocity.plot(hist=True, color=self.color, alpha=0.6)
            ax.set_title('Velocity\n{:5g}±{:5g}'.format(
                self.velocity.mean, self.velocity.std))
        ax.set_xlabel('m/s')

    def hilbert_angle(self):
        """get the angle of the hilbert transform"""
        return np.angle(hilbert(self.picked_signal)) / np.pi

    def hilbert_abs(self):
        """get the absolute value of the hilbert transform"""
        return np.abs(hilbert(self.picked_signal))

    def time_picks(self):
        """return picked time data"""
        if self.picks:
            time_picks = np.array(self.picks)
            if len(time_picks) == 1:
                return time_picks[0], time_picks[0], time_picks[0], 1e-50
            return np.min(time_picks), np.max(time_picks), np.mean(time_picks), np.std(time_picks)
        return -1, 1, 0, 0.25
