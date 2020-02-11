"""
Plot
====

Module providing the plotting capabilities.
"""

from contextlib import suppress
import itertools

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

from .isotopes import get_A_Z
from .config import logger


# Some constants for producing plots
LINESTYLES = ('-', '--', '-.', ':')
COLORS = ('b', 'g', 'r', 'c', 'm', 'k')
MARKERS = ('^', 's', '*', 'h', '+', 'x', 'o', 'v', 'p', 'd', '<')
MARKEVERY = (7, 19, 11, 17, 13)  # prime numbers to define where to put markers


def color_scheme_generator():
    """
    Function for creating unique combinations of style, color, and mark
    that is suitable for color-blind people.

    :returns: Generator of tuples (style, color, marker, markevery)
    :rtype: generator
    """

    # First take combinations for markers and linsetyles
    # which are the important elements for color-blind people
    combinations = np.array([_ for _ in itertools.product(MARKERS, LINESTYLES)])

    # This dictionary contains the following:
    #   - key: combination
    #   - value: penalty given the combinations that have already been used
    #
    # E.g, we do not want to use '^' if it has already been used twice and 's' only once
    # High penalty means we don't want to use this combination
    # The penalty is calculated the following way:
    #   - +1 for each time the marker has already been used
    #   - +1 for each time the linestyle has already been used
    penalties = {tuple(k): 0 for k in combinations}

    sorted_combinations = []
    while penalties:

        # Get combination with minimum penalty and save it
        min_c = min(penalties.keys(), key=lambda k: penalties[k])
        sorted_combinations.append(min_c)

        # Remove it and update penalties for the remaining combinations
        penalties.pop(min_c)
        for c in penalties:
            if c[0] == min_c[0]:
                penalties[c] += 1
            if c[1] == min_c[1]:
                penalties[c] += 1

    # Add colors, cycling through them
    combinations = ((c, combi[0], combi[1])
                    for c, combi in zip(itertools.cycle(COLORS), sorted_combinations))

    # Add markevery element to make sure markers do not overlap
    combinations = ((m,) + combi
                    for m, combi in zip(itertools.cycle(MARKEVERY), combinations))

    # Return a cycle in case we exhaust all combinations - e.g. abundance plot
    return itertools.cycle(combinations)


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

    def get_cattr(self, name, cycles=None):
        """Return header data from specific cycles.

        ..note:: Must be implemented in the inherited class.
        """
        raise NotImplementedError

    def get_dcol(self, name, cycles):
        """Return column data from specific cycles.

        ..note:: Must be implemented in the inherited class.
        """
        raise NotImplementedError

    def _find_data_by_correct_name(self, names, return_data=True):
        """Utility function that finds the correct data name and extracts the data if requested."""

        for name in names:
            with suppress(KeyError, ValueError):
                data = self.get_cattr(name)
                return (name, data) if return_data else name
        logger.warning(
            'Cannot find data with any of the following names: %s', names)

    @staticmethod
    def _slice_array(x, x0=None):
        """Returns indices for which x>x0 if x0, else an ellipsis."""
        if x0:
            return np.where(x > x0)
        return ...

    def plot(self, x, ys, cycles=None, x0=None, logx=False, logy=False,
             xlabel=None, ylabel=None, legend=True, show=True):
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
        :param bool legend: if True, display legend
        :param bool show: if True, display figure
        :returns: Figure
        :rtype: :py:class:`matplotlib.figure.Figure`
        """

        # Plot function lookup
        plot_function = getattr(plt, self.PLOT_FUNC_CHOICES[(logx, logy)])

        # Style generator
        style_gen = color_scheme_generator()

        # Default x label
        xlabel = xlabel or x

        # Profiles names should be iterable
        ys = [ys] if isinstance(ys, str) else ys

        # Figure
        fig = plt.figure()

        # Distinguish between:
        #  * singleton header data accross cycles
        #  * data for given cycles
        if cycles is None:
            x_data = self.get_cattr(x)
            indices = self._slice_array(x_data, x0)
            x_data = x_data[indices]
            for y in ys:
                y_data = self.get_cattr(y)[indices]
                step, color, marker, ls = next(style_gen)
                plot_function(x_data, y_data, color + ls + marker, markevery=step, label=y)
        else:
            cycles = [cycles] if isinstance(cycles, int) else cycles
            # We report the cycle number in the legend if more than one is plotted
            add_cycle_in_legend = len(cycles) > 1
            for c in cycles:
                x_data = self.get_dcol(x, c)
                indices = self._slice_array(x_data, x0)
                x_data = x_data[indices]
                for y in ys:
                    y_data = self.get_dcol(y, c)[indices]
                    step, color, marker, ls = next(style_gen)
                    label = y
                    if add_cycle_in_legend:
                        label += ', cycle {}'.format(c)
                    plot_function(x_data, y_data, color + ls + marker, markevery=step, label=label)

        # Axes labels
        plt.xlabel(xlabel)
        if ylabel:
            plt.ylabel(ylabel)

        # Legend
        if legend:
            plt.legend()

        # Show figure?
        if show:
            fig.show()
        return fig


class MesaPlotMixin(PlotMixin):
    """Class with additional plotting functionalities for MESA data."""

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
    KIPPENHAHN_PLOT_KWARGS = {
        'Mix1 Bot': {'c': 'b', 'marker': 'o', 'ls': 'None', 'alpha': 0.3, 'label': 'convection zones'},
        'Mix1 Top': {'c': 'b', 'marker': 'o', 'ls': 'None', 'alpha': 0.3, 'label': None},
        'Mix2 Bot': {'c': 'b', 'marker': 'o', 'ls': 'None', 'alpha': 0.3, 'label': None},
        'Mix2 Top': {'c': 'b', 'marker': 'o', 'ls': 'None', 'alpha': 0.3, 'label': None},
    }

    def hrd(self, show=True):
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
                        show=False, legend=False)

        # Invert x-axis and return
        fig.gca().invert_xaxis()

        # Show figure?
        if show:
            fig.show()
        return fig

    def tcrhoc(self, show=True):
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
                        show=False, legend=False)

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

        x_data = self.get_cattr(x)

        # Style generator
        style_gen = color_scheme_generator()

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
            x_data = x_data[indices]
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
            if k in self.KIPPENHAHN_PLOT_KWARGS:
                if 'label' in self.KIPPENHAHN_PLOT_KWARGS[k]:
                    plt.plot(x_data, v, **self.KIPPENHAHN_PLOT_KWARGS[k], lw=2)
                else:
                    plt.plot(x_data, v, **self.KIPPENHAHN_PLOT_KWARGS[k], label=k, lw=2)
            else:
                step, color, marker, ls = next(style_gen)
                plt.plot(x_data, v, color + ls + marker, markevery=step, label=k, lw=2)

        # Labels
        plt.xlabel(x)
        plt.ylabel(r'$m/M_\odot$')

        # Legend
        plt.legend()

        # Add y-axis if C/O ratio is plotted
        if CO_ratio:
            plt.twinx()
            step, color, marker, ls = next(style_gen)
            plt.plot(x_data, CO_ratio_data, color + ls + marker,
                     markevery=step, lw=2, label='C/O ratio')
            plt.ylabel('C/O ratio')
            plt.legend()

        # Show figure?
        if show:
            fig.show()
        return fig


class NugridPlotMixin(PlotMixin):
    """Class with additional plotting functionalities for NuGrid data."""

    def abu_profile(self, cycle, isos=None, logx=False, logy=False, show=True):
        """Abundances profiles vs mass coordinates.

        :param int cycle: cycle for which to plot abundances
        :param list(str) isos: list of isotopes to plot. Default: all isotopes
        :param bool show: if True, display figure
        :returns: Figure
        :rtype: :py:class:`matplotlib.figure.Figure`
        """

        isos = isos or self.isotopes

        return self.plot('mass', isos, cycles=cycle,
                         xlabel=r'$m/M_\odot$',
                         ylabel=r'$\log(X)$',
                         logx=logx, logy=logy,
                         show=show, legend=True)

    def iso_abu(self, cycle, a_min=None, a_max=None, threshold=1e-12, show=True):
        """Isotopes abundances plot.

        :param int cycle: cycle for which to plot abundances
        :param int a_min: minimum mass number to consider
        :param int a_max: maximum mass number to consider
        :param float threshold: mass fraction threshold below which element is not plotted
        :param bool show: if True, display figure
        :returns: Figure
        :rtype: :py:class:`matplotlib.figure.Figure`
        """

        a_min = a_min or min(self.A)
        a_max = a_max or max(self.A)

        # Style generator
        style_gen = color_scheme_generator()

        # Data to plot is organized in a dictionary
        # Key is the element name
        # Value is a tuple of lists ([A1, A2,...], [logmassf1, logmassf2...])
        data_to_plot = {}

        for iso in self.isotopes:
            A, _, element = get_A_Z(iso)

            # Filter those not wanted
            if A < a_min or A > a_max:
                continue

            # Get mass fraction - we take the surface value for the moment
            massf = self.get_dcol(iso, cycle)[0]
            if massf < threshold:
                continue
            log_massf = np.log10(massf)

            # Put in the dictrionary
            if element not in data_to_plot:
                data_to_plot[element] = ([A], [log_massf])
            else:
                data_to_plot[element][0].append(A)
                data_to_plot[element][1].append(log_massf)

        # Figure
        fig = plt.figure()

        for element in data_to_plot:
            x = data_to_plot[element][0]
            y = data_to_plot[element][1]

            _, color, marker, ls = next(style_gen)
            plt.plot(x, y, color + ls + marker)

            # Annotate at the largest mass fraction
            x_max, y_max = sorted(zip(x, y), key=lambda e: e[1], reverse=True)[0]
            plt.text(x_max, y_max + 0.1, element)

        plt.xlabel('Mass Number ($A$)')
        plt.ylabel('log(mass fraction)')

        # Show figure?
        if show:
            fig.show()
        return fig

    def abu_chart(self, cycle, n_min=None, n_max=None, z_min=None, z_max=None,
                  threshold=1e-12, show=True):
        """Isotopic abundances chart.

        :param int cycle: cycle for which to plot abundances
        :param int n_min: minimum neutron number to consider
        :param int n_max: maximum neutron number to consider
        :param int z_min: minimum atomic number to consider
        :param int z_max: maximum atomic number to consider
        :param float threshold: mass fraction threshold below which element is not plotted
        :param bool show: if True, display figure
        :returns: Figure
        :rtype: :py:class:`matplotlib.figure.Figure`
        """

        n_min = n_min or min(self.A - self.Z)
        n_max = n_max or max(self.A - self.Z)
        z_min = z_min or min(self.Z)
        z_max = z_max or max(self.Z)

        # Data to plot is organized in a dictionary
        # Key is the element name
        # Value is a list of tuples, where each tuple has the form ((N, Z), log_massf)
        data_to_plot = {}

        for iso in self.isotopes:
            A, Z, element = get_A_Z(iso)
            N = A - Z

            # Filter those not wanted
            if N < n_min or N > n_max or Z < z_min or Z > z_max:
                continue

            # Get mass fraction - we take the surface value for the moment
            massf = self.get_dcol(iso, cycle)[0]
            if massf < threshold:
                log_massf = None
            else:
                log_massf = np.log10(massf)

            # Put in the dictrionary
            coordinates = (N, Z)
            if element not in data_to_plot:
                data_to_plot[element] = [(coordinates, log_massf)]
            else:
                data_to_plot[element].append((coordinates, log_massf))

        # Figure
        fig = plt.figure()

        # Axes
        # -/+1 for the isotpe name and squares to appear properly
        fig.gca().set_xlim(n_min - 1, n_max + 1)
        fig.gca().set_ylim(z_min - 1, z_max + 1)
        plt.xlabel('Neutron number ($A-Z$)')
        plt.ylabel('Atomic number ($Z$)')
        plt.axis('equal')

        # Colormap
        cmap = plt.get_cmap('jet')
        norm = mpl.colors.Normalize(vmin=np.log10(threshold), vmax=0)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm)
        cbar.ax.set_ylabel(r'$\log_{10}(X)$')

        # Plot data
        for element, isos in data_to_plot.items():
            for iso in isos:
                (N, Z), log_massf = iso

                # Draw square and annotate A
                color = cmap(norm(log_massf)) if log_massf else 'none'
                rect = patches.Rectangle(
                    (N - 0.5, Z - 0.5), 1, 1, linewidth=2, edgecolor='k', facecolor=color)
                fig.gca().add_patch(rect)
                plt.annotate(
                    N + Z, (N, Z), horizontalalignment='center', verticalalignment='center')

            # Add isotope name
            N_min_element = min([iso[0][0] for iso in isos])
            plt.annotate(element, (N_min_element - 1, Z),
                         horizontalalignment='center', verticalalignment='center')

        # Title
        plt.title('Isotropic chart for cycle {}'.format(cycle))

        # Show figure?
        if show:
            fig.show()
        return fig
