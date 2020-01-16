"""
Plot
====

Module providing the plotting capabilities.
"""

from contextlib import suppress
import logging

import matplotlib.pyplot as plt
import numpy as np


logging.basicConfig(level=logging.INFO)


class PlotMixin:
    """Abstract class giving plotting functionalities."""

    # These attributes allow us to proxy the actualy methods extracting header and cycle data.
    # They should be defined in the inherited class
    _get_data_from_headers_name = None
    _get_data_from_cycles_name = None

    # Matplotlib functions to use for plotting (plot, semilogx, ...)
    # Key is a tuple (logx,logy)
    PLOT_FUNC_CHOICES = {
        (False, False): 'plot',
        (True, False): 'semilogx',
        (False, True): 'semilogy',
        (True, True): 'loglog',
    }

    def get_cycle_header(self, name, cycles=None):
        """Return header data from specific cycles.

        ..note:: Must be implemented in the inherited class.
        """
        raise NotImplementedError

    def get_cycle_data(self, name, cycles):
        """Return data from specific cycles.

        ..note:: Must be implemented in the inherited class.
        """
        raise NotImplementedError

    def plot(self, x, ys, cycles=None, logx=False, logy=False,
             xlabel=None, ylabel=None, legend=None, show=True, **kwargs):
        """Generate a 2D plot.

        :param str x: name of the data for the x-axis
        :param ys: name of the data for the y-axis
        :type ys: str or list(str)
        :param cycles: requested cycles
        :type cycles: int or list(int)
        :param bool logx: if True, use logarithmic x-axis
        :param bool logy: if True, use logarithmic y-axis
        :param str xlabel: label for the x-axis
        :param str ylabel: label for the y-axis
        :param tuple(str) legend: legend to be added
        :param bool show: if True, display figure
        :returns: Figure
        :rtype: :py:class:`matplotlib.figure.Figure`
        """

        # Profiles names should be iterable
        ys = [ys] if isinstance(ys, str) else ys

        # Figure
        fig = plt.figure()

        # Distinguish between:
        #  * singleton header data accross cycles
        #  * data for given cycles
        if not cycles:
            x_data = self.get_cycle_header(x)
            for y in ys:
                y_data = self.get_cycle_header(y)
                getattr(plt, self.PLOT_FUNC_CHOICES[(logx, logy)])(x_data, y_data, **kwargs)
        else:
            cycles = [cycles] if isinstance(cycles, int) else cycles
            for c in cycles:
                x_data = self.get_cycle_data(x, c)
                for y in ys:
                    y_data = self.get_cycle_data(y, c)
                    getattr(plt, self.PLOT_FUNC_CHOICES[(logx, logy)])(x_data, y_data, **kwargs)

        # Axes labels
        if xlabel:
            plt.xlabel(xlabel)
        if ylabel:
            plt.ylabel(ylabel)

        # Legend
        if legend:
            plt.legend(legend)

        # Show figure?
        if show:
            fig.show()

        return fig


class StarPlotsMixin(PlotMixin):
    """Class with additional plotting functionalities for stellar data."""

    # Data is accessed through different names depending on its type
    # We do not want to have several obscured if...else statements in the plot functions
    # with several hidden and hard-coded names
    # Instead, we classify these possible names in class attributes and use a utility function
    #
    # The keys is name representing the content of the data, values all the possible names to access it
    LOG_TEFF_NAMES = ('logTeff', 'log_Teff')
    LOG_L_NAMES = ('logL', 'log_L')
    STELLAR_MASS_NAMES = ('mass', 'star_mass', 'total_mass')
    HE_CORE_NAMES = ('h1_boundary_mass', 'he_core_mass')
    C_CORE_NAMES = ('he4_boundary_mass', 'c_core_mass')
    O_CORE_NAMES = ('c12_boundary_mass', 'o_core_mass')

    CORE_NAMES = {
        'He core': HE_CORE_NAMES,
        'C core': C_CORE_NAMES,
        'O core': O_CORE_NAMES,
    }
    MIXING_BOUNDARIES_NAMES = {
        'Mix1 Bot': ('mx1_bot',),
        'Mix1 Top': ('mx1_top',),
        'Mix2 Bot': ('mx2_bot',),
        'Mix2 Top': ('mx2_top',),
    }

    # Plot kwargs for the different types of data.
    PLOT_KWARGS = {
        'He core': {'c': 'k'},
        'C core': {'c': 'g', 'ls': '--'},
        'O core': {'c': 'c', 'ls': ':'},
        'Mix1 Bot': {'c': 'b', 'marker': 'o', 'ls': 'None', 'alpha': 0.3, 'label': 'convection zones'},
        'Mix1 Top': {'c': 'b', 'marker': 'o', 'ls': 'None', 'alpha': 0.3, 'label': None},
        'Mix2 Bot': {'c': 'b', 'marker': 'o', 'ls': 'None', 'alpha': 0.3, 'label': None},
        'Mix2 Top': {'c': 'b', 'marker': 'o', 'ls': 'None', 'alpha': 0.3, 'label': None},
    }

    def _find_data_by_correct_name(self, names):
        """Utility function that finds the correct data name and extracts the data."""

        for name in names:
            with suppress(KeyError, ValueError):
                return name, self.get_cycle_header(name)
        logging.warning(
            'Cannot find data with any of the following names: %s', names)

    def hrd(self, show=True, **kwargs):
        """Hertzsprung-Russel diagramm.
        :param bool show: if True, display figure
        :returns: Figure
        :rtype: :py:class:`matplotlib.figure.Figure`
        """

        _, x_data = self._find_data_by_correct_name(self.LOG_TEFF_NAMES)
        _, y_data = self._find_data_by_correct_name(self.LOG_L_NAMES)

        fig = plt.figure()
        plt.plot(x_data, y_data, **kwargs)

        plt.xlabel(r'$\log(T_{\mathrm{eff}}/T_\odot)$')
        plt.ylabel(r'$\log(L/L_\odot)$')

        # Invert x-axis and return
        fig.gca().invert_xaxis()

        # Show figure?
        if show:
            fig.show()

        return fig

    def kippenhahn(self, x, x0=None, show=True):
        """Kippenhahn diagramm as a function of time or model.

        :param str x: name of the data for the x-axis. Should be either 'star_age' or 'model'.
        :param x0: zero point for the plot. Previous data is not shown. Must be same unit as `x`.
        :type x0: int or float
        :param tuple(str) legend: legend to be added
        :param bool show: if True, display figure
        :returns: Figure
        :rtype: :py:class:`matplotlib.figure.Figure`
        """

        x_data = self.get_cycle_header(x)

        # Data using the utility function
        # We should at least have the mass
        y_data = self._find_data_by_correct_name(self.STELLAR_MASS_NAMES)
        if not y_data:
            raise RuntimeError('Cannot extract stellar mass from the data.')

        # Key for the star mass - will be needed for the mixing boundaries
        star_mass_name = y_data[0]
        y_data = {star_mass_name: y_data[1]}  # change to dict

        # Mixing zones
        for k, zone in self.MIXING_BOUNDARIES_NAMES.items():
            tmp_y = self._find_data_by_correct_name(zone)
            if tmp_y:
                y_data.update([(k, tmp_y[1] * y_data[star_mass_name])])

        # Cores
        for k, core in self.CORE_NAMES.items():
            tmp_y = self._find_data_by_correct_name(core)
            if tmp_y:
                y_data.update([(k, tmp_y[1])])

        # Shift origin and keep only relevant data
        if x0:
            indices = np.where(x_data > x0)
            x_data = x_data[indices] - x0
            y_data = {k: v[indices] for k, v in y_data.items()}

        # Figure
        fig = plt.figure()

        # Plot
        for k, v in y_data.items():
            try:
                if 'label' in self.PLOT_KWARGS[k]:
                    plt.plot(x_data, v, **self.PLOT_KWARGS[k], lw=2)
                else:
                    plt.plot(x_data, v, **self.PLOT_KWARGS[k], label=k, lw=2)
            except KeyError:
                logging.debug('Cannot find specific plotting kwargs for %s', k)
                plt.plot(x_data, v, label=k, lw=2)

        # Labels
        plt.xlabel(x)
        plt.ylabel(r'$m/M_\odot$')

        # Legend
        plt.legend()

        # Show figure?
        if show:
            fig.show()

        return fig
