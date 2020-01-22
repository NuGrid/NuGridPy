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

    def _find_data_by_correct_name(self, names, return_data=True):
        """Utility function that finds the correct data name and extracts the data if requested."""

        for name in names:
            with suppress(KeyError, ValueError):
                data = self.get_cycle_header(name)
                return (name, data) if return_data else name
        logging.warning(
            'Cannot find data with any of the following names: %s', names)

    @staticmethod
    def _slice_array(x, x0=None):
        """Returns indices for which x>x0 if x0, else an ellipsis."""
        if x0:
            return np.where(x > x0)
        return ...

    def plot(self, x, ys, cycles=None, x0=None, logx=False, logy=False,
             xlabel=None, ylabel=None, legend=None, show=True, **kwargs):
        """Generate a 2D plot.

        :param str x: name of the data for the x-axis
        :param ys: name of the data for the y-axis
        :type ys: str or list(str)
        :param cycles: requested cycles
        :type cycles: int or list(int)
        :param x0: if provided, only plots data for x > x0. Must be same unit as `x`.
        :type x0: int or float
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
            indices = self._slice_array(x_data, x0)
            x_data = x_data[indices]
            for y in ys:
                y_data = self.get_cycle_header(y)[indices]
                getattr(plt, self.PLOT_FUNC_CHOICES[(logx, logy)])(x_data, y_data, **kwargs)
        else:
            cycles = [cycles] if isinstance(cycles, int) else cycles
            for c in cycles:
                x_data = self.get_cycle_data(x, c)
                indices = self._slice_array(x_data, x0)
                x_data = x_data[indices]
                for y in ys:
                    y_data = self.get_cycle_data(y, c)[indices]
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
    STELLAR_MASS_NAMES = ('mass', 'star_mass', 'total_mass')
    LOG_TEFF_NAMES = ('logTeff', 'log_Teff')
    LOG_L_NAMES = ('logL', 'log_L')
    LOG_RHOC_NAMES = ('log_center_Rho',)
    LOG_TC_NAMES = ('log_center_T',)

    SURFACE_C12_NAMES = ('surface_c12',)
    SURFACE_O16_NAMES = ('surface_o16',)

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

    def hrd(self, show=True, **kwargs):
        """Hertzsprung-Russel diagramm.

        :param bool show: if True, display figure
        :returns: Figure
        :rtype: :py:class:`matplotlib.figure.Figure`
        """

        teff_name = self._find_data_by_correct_name(self.LOG_TEFF_NAMES, return_data=False)
        l_name = self._find_data_by_correct_name(self.LOG_L_NAMES, return_data=False)
        if (not teff_name) or (not l_name):
            raise RuntimeError("Cannot access temperature or luminosity data.")

        fig = self.plot(teff_name, l_name,
                        xlabel=r'$\log(T_{\mathrm{eff}}/T_\odot)$',
                        ylabel=r'$\log(L/L_\odot)$',
                        show=False, **kwargs)

        # Invert x-axis and return
        fig.gca().invert_xaxis()

        # Show figure?
        if show:
            fig.show()
        return fig

    def tcrhoc(self, show=True, **kwargs):
        """Central temperature against central density plot.

        :param bool show: if True, display figure
        :returns: Figure
        :rtype: :py:class:`matplotlib.figure.Figure`
        """

        rhoc_name = self._find_data_by_correct_name(self.LOG_RHOC_NAMES, return_data=False)
        tc_name = self._find_data_by_correct_name(self.LOG_TC_NAMES, return_data=False)

        if (not rhoc_name) or (not tc_name):
            raise RuntimeError("Cannot access central temperature or density data.")

        fig = self.plot(rhoc_name, tc_name,
                        xlabel=r'$\log(\rho_c/\,[g.cm^{-3}])$',
                        ylabel=r'$\log(T_c/[K])$',
                        show=False, **kwargs)

        # Invert x-axis and return
        fig.gca().invert_xaxis()

        # Show figure?
        if show:
            fig.show()
        return fig

    def kippenhahn(self, x, x0=None, CO_ratio=False, show=True):
        """Kippenhahn diagramm as a function of time or model.

        :param str x: name of the data for the x-axis. Should be either 'star_age' or 'model'.
        :param x0: zero point for the plot. Previous data is not shown. Must be same unit as `x`.
        :type x0: int or float
        :param tuple(str) legend: legend to be added
        :param bool CO_ratio: if True, display C/O ratio
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

        # C/O ratio
        if CO_ratio:
            c12 = self._find_data_by_correct_name(self.SURFACE_C12_NAMES)
            o16 = self._find_data_by_correct_name(self.SURFACE_O16_NAMES)
            if c12 and o16:
                CO_ratio_data = c12[1] / o16[1]
                if x0:
                    CO_ratio_data = CO_ratio_data[indices]
            else:
                CO_ratio = False

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

        # Add y-axis if C/O ratio is plotted
        if CO_ratio:
            plt.twinx()
            plt.plot(x_data, CO_ratio_data, 'k--', lw=2, label='C/O ratio')
            plt.ylabel('C/O ratio')
            plt.legend()

        # Show figure?
        if show:
            fig.show()
        return fig
