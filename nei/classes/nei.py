import numpy as np
import matplotlib.pyplot as plt
from typing import Union, Optional, List, Dict, Callable
import astropy.units as u
import plasmapy as pl
import collections
from scipy import interpolate
from .eigenvaluetable import EigenData2
from .ionization_states import IonizationStates

# TODO: Allow this to keep track of velocity and position too, and
# eventually to have density and temperature be able to be functions of
# position.  (and more complicated expressions for density and
# temperature too)

# TODO: Expand Simulation docstring


class NEIError(Exception):
    pass


class Simulation:
    """
    Store results from a non-equilibrium ionization simulation.
    """
    def __init__(self, initial, n_init, T_e_init, max_steps, time_start):

        self._elements = initial.elements
        self._abundances = initial.abundances
        self._max_steps = max_steps

        self._nstates = {elem: pl.atomic.atomic_number(elem) + 1 for elem in self.elements}

        self._ionic_fractions = {
            elem: np.full((self.nstates[elem], max_steps + 1), np.nan, dtype=np.float64)
            for elem in self.elements
        }

        self._number_densities = {
            elem: np.full((self.nstates[elem], max_steps + 1), np.nan, dtype=np.float64) * u.cm ** -3
            for elem in self.elements
        }

        self._n_elem = {elem: np.full(max_steps + 1, np.nan) * u.cm ** -3 for elem in self.elements}
        self._n_e = np.full(max_steps + 1, np.nan) * u.cm ** -3
        self._T_e = np.full(max_steps + 1, np.nan) * u.K
        self._time = np.full(max_steps + 1, np.nan) * u.s

        self._index = 0

        self._assign(
            new_time=time_start,
            new_ionic_fractions=initial.ionic_fractions,
            new_n = n_init,
            new_T_e = T_e_init,
        )

    def _assign(self, new_time, new_ionic_fractions, new_n, new_T_e):

        try:
            self._time[self._index] = new_time
            self._T_e[self._index] = new_T_e

            for elem in self.elements:
                self._ionic_fractions[elem][:, self._index] = new_ionic_fractions[elem][:]

            # Calculate elemental and ionic number densities
            n_elem = {elem: new_n * self.abundances[elem] for elem in self.elements}
            number_densities = {
                elem: n_elem[elem] * new_ionic_fractions[elem]
                for elem in self.elements
            }

            # Calculate the electron number density
            n_e = 0.0 * u.cm ** -3
            for elem in self.elements:
                integer_charges = np.linspace(0, self.nstates[elem]-1, self.nstates[elem])
                n_e += np.sum(number_densities[elem] * integer_charges)

            # Assign densities
            self._n_e[self._index] = n_e
            for elem in self.elements:
                self._n_elem[elem][self._index] = n_elem[elem]
                self._number_densities[elem][:, self._index] = number_densities[elem]

        except Exception as exc:
            raise NEIError(
                f"Unable to assign parameters to Simulation instance "
                f"for index {self._index} at time = {new_time}.  The "
                f"parameters are new_n = {new_n}, new_T_e = {new_T_e}, "
                f"and new_ionic_fractions = {new_ionic_fractions}."
            ) from exc
        finally:
            self._index += 1

    @property
    def max_steps(self):
        return self._max_steps

    @property
    def nstates(self):
        return self._nstates

    @property
    def elements(self):
        return self._elements

    @property
    def abundances(self):
        return self._abundances

    @property
    def ionic_fractions(self):
        return self._ionic_fractions

    @property
    def number_densities(self):
        return self._number_densities

    @property
    def n_elem(self):
        return self._n_elem

    @property
    def n_e(self):
        return self._n_e

    @property
    def T_e(self):
        return self._T_e

    @property
    def time(self):
        return self._time


class NEI:
    r"""
    Perform and analyze a non-equilibrium ionization simulation.

    Parameters
    ----------
    inputs

    T_e: `~astropy.units.Quantity` or `callable`
        The electron temperature, which may be a constant, an array of
        temperatures corresponding to the times in `time_input`, or a
        function that yields the temperature as a function of time.

    n: `~astropy.units.Quantity` or `callable`
        The number density multiplicative factor.  The number density of
        each element will be `n` times the abundance given in
        `abundances`.  For example, if `abundance['H'] = 1`, then this
        will correspond to the number density of hydrogen (including
        neutral hydrogen and protons).  This factor may be a constant,
        an array of number densities over time, or a function that
        yields a number density as a function of time.

    time_input: ~astropy.units.Quantity, optional
        An array containing the times associated with `n` and `T_e` in
        units of time.

    time_start: ~astropy.units.Quantity, optional
        The start time for the simulation.  If density and/or
        temperature are given by arrays, then this argument must be
        greater than `time_input[0]`.  If this argument is not supplied,
        then `time_start` defaults to `time_input[0]` (if given) and
        zero seconds otherwise.

    time_max: ~astropy.units.Quantity
        The maximum time for the simulation.  If density and/or
        temperature are given by arrays, then this argument must be less
        than `time_input[-1]`.

    verbose: bool, optional
        A flag stating whether or not to print out information for every
        time step. Setting `verbose` to `True` is useful for testing.
        Defaults to `False`.

    abundances: dict

    Examples
    --------

    >>> import numpy as np
    >>> import astropy.units as u

    >>> inputs = {'H': [0.9, 0.1], 'He': [0.9, 0.099, 0.001]}
    >>> n = 1e9 * u.m ** -3
    >>> abundances = {'H': 1, 'He': 0.085}
    >>> T_e = np.array([5000, 50000]) * u.K
    >>> time_input = np.array([0, 10]) * u.min

    The initial conditions can be accessed using the initial attribute.

    #>>> nei.initial['H']

    After having inputted all of the necessary information, we can run
    the simulation.

    #>>> nei.simulate()

    The final results can be access with the `final` attribute.

    #>>> nei.final['H']
    #array([0.0, 1.0])
    #>>> nei.final.T_e
    #<Quantity 50000. K>

    """

    def __init__(
            self,
            inputs,
            abundances: Union[Dict, str] = None,
            T_e: Union[Callable, u.Quantity] = None,
            n: Union[Callable, u.Quantity] = None,
            time_input: u.Quantity = None,
            time_start: u.Quantity = None,
            time_max: u.Quantity = None,
            max_steps: Union[int, np.integer] = 1000,
            tol: Union[int, float] = 1e-15,
            dt: u.Quantity = None,
            adapt_dt: bool = None,
            safety_factor: Union[int, float] = 1,
            verbose: bool = False,
    ):

        try:

            self.time_input = time_input
            self.time_start = time_start
            self.time_max = time_max
            self.T_e_input = T_e
            self.n_input = n
            self.max_steps = max_steps
            self.dt_input = dt
            self._dt = dt
            self.adapt_dt = adapt_dt
            self.safety_factor = safety_factor
            self.verbose = verbose

            T_e_init = self.electron_temperature(self.time_start)
            n_init = self.hydrogen_number_density(self.time_start)

            self.initial = IonizationStates(
                inputs=inputs,
                abundances=abundances,
                T_e=T_e_init,
                n_H=n_init,  # TODO: Update n_H in IonizationState(s)
                tol = tol
            )

            self.tol = tol
            self.elements = self.initial.elements

            if 'H' not in self.elements:
                raise NEIError("Must have H in elements")

            self.abundances = self.initial.abundances

            self._EigenDataDict = {element: EigenData2(element) for element in self.elements}

            if self.T_e_input is not None and not isinstance(inputs, dict):
                for element in self.initial.elements:
                    self.initial.ionic_fractions[element] = \
                        self.EigenDataDict[element].equilibrium_state(T_e_init.value)

        except Exception:
            raise NEIError(
                f"Unable to create NEI instance for:\n"
                f"     inputs = {inputs}\n"
                f" abundances = {abundances}\n"
                f"        T_e = {T_e}\n"
                f"          n = {n}\n"
                f" time_input = {time_input}\n"
                f" time_start = {time_start}\n"
                f"   time_max = {time_max}\n"
                f"  max_steps = {max_steps}\n"
            )

    def equil_ionic_fractions(self, T_e=None, time=None) -> dict:
        """
        Return the equilibrium ionic fractions for a temperature or at
        a given time.

        Parameters
        ----------
        T_e: ~astropy.units.Quantity, optional
            The electron temperature in units that can be converted to
            kelvin.

        time: ~astropy.units.Quantity, optional
            The time in units that can be converted to seconds.

        Returns
        -------
        equil_ionfracs: dict
            The equilibrium ionic fractions for the elements contained
            within this class

        Notes
        -----
        Only one of `T_e` and `time` may be included as an argument.  If
        neither `T_e` or `time` is provided and the temperature for the
        simulation is given by a constant, the this method will assume
        that `T_e` is the temperature of the simulation.

        """

        if T_e is not None and time is not None:
            raise NEIError("Only one of T_e and time may be used as an argument.")

        if T_e is None and time is None:
            if self.T_e_input.isscalar:
                T_e = self.T_e_input
            else:
                raise NEIError

        try:
            T_e = T_e.to(u.K) if T_e is not None else None
            time = time.to(u.s) if time is not None else None
        except Exception as exc:
            raise NEIError("Invalid input to equilibrium_ionic_fractions.")

        if time is not None:
            T_e = self.electron_temperature(time)

        if not T_e.isscalar:
            raise NEIError("Need scalar input for equilibrium_ionic_fractions.")

        equil_ionfracs = {}
        for element in self.elements:
            equil_ionfracs[element] = self.EigenDataDict[element].equilibrium_state(T_e.value)

        return equil_ionfracs

    @property
    def elements(self):
        return self._elements

    @elements.setter
    def elements(self, elements):
        # TODO: Update this
        self._elements = elements

    @property
    def abundances(self):
        return self._abundances

    @abundances.setter
    def abundances(self, abund):
        self._abundances = abund

    @property
    def tol(self):
        return self._tol

    @tol.setter
    def tol(self, value):
        self._tol = value

    @property
    def time_input(self):
        return self._time_input

    @time_input.setter
    def time_input(self, times):
        if times is None:
            self._time_input = None
        elif isinstance(times, u.Quantity):
            if times.isscalar:
                raise ValueError("time_input must be an array.")
            try:
                times = times.to(u.s)
            except u.UnitConversionError:
                raise u.UnitsError("time_input must have units of seconds.") from None
            if not np.all(times[1:] > times[:-1]):
                raise ValueError("time_input must monotonically increase.")
            self._time_input = times
        else:
            raise TypeError("Invalid time_input.")

    @property
    def time_start(self):
        return self._time_start

    @time_start.setter
    def time_start(self, time):
        if time is None:
            self._time_start = 0.0 * u.s
        elif isinstance(time, u.Quantity):
            if not time.isscalar:
                raise ValueError("time_start must be a scalar")
            try:
                time = time.to(u.s)
            except u.UnitConversionError:
                raise u.UnitsError("time_start must have units of seconds") from None
            if hasattr(self, '_time_max') \
                    and self._time_max is not None and self._time_max<=time:
                raise ValueError("time_start must be less than time_max")
            if self.time_input is not None and self.time_input.min() > time:
                raise ValueError("time_start must be less than min(time_input)")
            self._time_start = time
        else:
            raise TypeError("Invalid time_start.") from None

    @property
    def time_max(self):
        return self._time_max

    @time_max.setter
    def time_max(self, time):
        if time is None:
            self._time_max = self.time_input[-1] if self.time_input is not None else None
        elif isinstance(time, u.Quantity):
            if not time.isscalar:
                raise ValueError("time_max must be a scalar")
            try:
                time = time.to(u.s)
            except u.UnitConversionError:
                raise u.UnitsError("time_max must have units of seconds") from None
            if hasattr(self, '_time_start') and self._time_start is not None and \
                    self._time_start >= time:
                raise ValueError("time_max must be greater than time_start")
            self._time_max = time
        else:
            raise TypeError("Invalid time_max.") from None

    @property
    def dt_input(self):
        return self._dt

    @dt_input.setter
    def dt_input(self, dt: Optional[u.Quantity]):
        if dt is None:
            self._dt_input = None
            self._dt = None
        elif isinstance(dt, u.Quantity):
            try:
                dt = dt.to(u.s)
                if dt > 0 * u.s:
                    self._dt_input = dt
            except (AttributeError, u.UnitConversionError):
                raise NEIError("Invalid dt.")

    @property
    def adapt_dt(self):
        return self._adapt_dt

    @adapt_dt.setter
    def adapt_dt(self, choice: Optional[bool]):
        if choice is None:
            self._adapt_dt = True if self.dt_input is None else False
        elif choice is True or choice is False:
            self._adapt_dt = choice
        else:
            raise TypeError("Invalid value for adapt_dt")

    @property
    def safety_factor(self):
        return self._safety_factor

    @safety_factor.setter
    def safety_factor(self, value):
        if not isinstance(value, (float, np.float64, np.float32, np.integer, int)):
            raise TypeError
        if 1e-3 <= value <= 1e3:
            self._safety_factor = value
        else:
            raise NEIError("Invalid safety factor.")

    @property
    def verbose(self):
        return self._verbose

    @verbose.setter
    def verbose(self, choice):
        if choice is True or choice is False:
            self._verbose = choice
        else:
            raise TypeError("Invalid choice for verbose.")

    def _check_time(self, time):
        try:
            time_s = time.to(u.s)
        except Exception as exc:
            raise NEIError(f"{time} is not a valid time.")
        if np.isnan(time_s):
            raise NEIError(f"time is not a number.")

        if time_s < self.time_start:
            raise NEIError(
                f"time = {time} is less than time_start = "
                f"{self.time_start}.")
        elif self._time_max is not None and time_s > self.time_max:
            raise NEIError(
                f"time = {time} is greater than time_max = "
                f"{self.time_max}"
            )

    @property
    def max_steps(self):
        return self._max_steps

    @max_steps.setter
    def max_steps(self, n):
        if isinstance(n, (int, np.integer)) and n > 0:
            self._max_steps = n
        else:
            raise TypeError("max_steps must be an integer")

    @property
    def T_e_input(self):
        return self._T_e_input

    @T_e_input.setter
    def T_e_input(self, T_e: Optional[Union[Callable, u.Quantity]]):
        """Set the input electron temperature."""
        if isinstance(T_e, u.Quantity):
            try:
                T_e = T_e.to(u.K, equivalencies=u.temperature_energy())
            except u.UnitConversionError:
                raise u.UnitsError("Invalid electron temperature.") from None
            if T_e.isscalar:
                self._T_e_input = T_e
                self._electron_temperature = lambda time: T_e
            else:
                if self._time_input is None:
                    raise TypeError("Must define time_input prior to T_e for an array.")
                time_input = self.time_input
                if len(time_input) != len(T_e):
                    raise ValueError("len(T_e) is not equal to len(time_input).")
                f = interpolate.interp1d(time_input.value, T_e.value)
                self._electron_temperature = lambda time: f(time.value) * u.K
                self._T_e_input = T_e
        elif callable(T_e):
            if self.time_start is not None:
                try:
                    T_e(self.time_start).to(u.K)
                    T_e(self.time_max).to(u.K)
                except Exception:
                    raise ValueError("Invalid electron temperature function.")
            self._T_e_input = T_e
            self._electron_temperature = T_e
        elif T_e is None:
            self._electron_temperature = lambda: None
        else:
            raise TypeError("Invalid T_e")

    def electron_temperature(self, time):
        try:
            self._check_time(time)
            time = time.to(u.s)
            T_e = self._electron_temperature(time).to(u.K)

            if np.isnan(T_e):
                raise NEIError(f"Finding T_e = {T_e} at time = {time}.")
            elif T_e < 0 * u.K:
                raise
            else:
                return T_e
        except Exception as exc:
            raise NEIError(f"Unable to calculate electron temperature for time {time}")

    @property
    def n_input(self) -> u.Quantity:
        """The hydrogen number density."""
        if 'H' in self.elements:
            return self._n_input
        else:
            raise ValueError

    @n_input.setter
    def n_input(self, n):
        if isinstance(n, u.Quantity):
            try:
                n = n.to(u.cm ** -3)
            except u.UnitConversionError:
                raise u.UnitsError("Invalid hydrogen density.")
            if n.isscalar:
                self._n_input = n
                self.hydrogen_number_density = lambda time: n
            else:
                if self._time_input is None:
                    raise TypeError("Must define time_input prior to n for an array.")
                time_input = self.time_input
                if len(time_input) != len(n):
                    raise ValueError("len(n) is not equal to len(time_input).")
                f = interpolate.interp1d(time_input.value, n.value)
                self._hydrogen_number_density = lambda time: f(time.value) * u.cm ** -3
                self._n_input = n
        elif callable(n):
            if self.time_start is not None:
                try:
                    n(self.time_start).to(u.cm ** -3)
                    n(self.time_max).to(u.cm ** -3)
                except Exception:
                    raise ValueError("Invalid hydrogen number density function.")
            self._n_input = n
            self._hydrogen_number_density = n
        elif n is None:
            self._hydrogen_number_density = lambda: None
        else:
            raise TypeError("Invalid n.")

    def hydrogen_number_density(self, time):
        try:
            time = time.to(u.s)
        except (AttributeError, u.UnitsError):
            raise NEIError("Invalid time in hydrogen_density")
        return self._hydrogen_number_density(time)

    @property
    def EigenDataDict(self):
        return self._EigenDataDict

    @EigenDataDict.setter
    def EigenDataDict(self):
        self._EigenDataDict = {element: EigenData2(element) for element in self.elements}

    @property
    def initial(self):
        """
        The ~plasmapy.atomic.IonizationStates instance representing the
        initial conditions of the simulation.
        """
        return self._initial

    @initial.setter
    def initial(self, initial_states: Optional[IonizationStates]):
        if isinstance(initial_states, IonizationStates):
            self._initial = initial_states
            self._elements = initial_states.elements
        elif initial_states is None:
            self._ionstates = None
        else:
            raise TypeError("Expecting an IonizationStates instance.")

    @property
    def results(self):
        try:
            return self._results
        except Exception:
            raise AttributeError("The simulation has not yet been performed.")

    @results.setter
    def results(self, value):
        if isinstance(value, Simulation):
            self._results = value
        else:
            raise TypeError

    @property
    def final(self):
        try:
            return self._final
        except Exception:
            raise AttributeError("The simulation has not yet been performed.")

    def _initialize_simulation(self):

        self._results = Simulation(
            initial=self.initial,
            n_init=self.hydrogen_number_density(self.time_start),
            T_e_init=self.electron_temperature(self.time_start),
            max_steps=self.max_steps,
            time_start=self.time_start,
        )

    def set_timestep(self, dt=None):
        if dt is not None:
            try:
                dt = dt.to(u.s)
            except Exception:
                raise NEIError("Invalid timestep.")
            finally:
                self._dt = dt
        elif self.adapt_dt:
            raise NotImplementedError(
                "Adaptive time step not yet implemented; set adapt_dt "
                "to False.")
        elif self.dt_input is not None:
            self._dt = self.dt_input
        else:
            raise NEIError("Unable to get set timestep.")

    def simulate(self):
        """
        Perform a non-equilibrium ionization simulation.
        """

        self._initialize_simulation()

        for step in range(self.max_steps):

            try:
                self.set_timestep()
                self.time_advance()
            except Exception as exc:
                raise NEIError(f"Unable to complete simulation.") from exc

        self._finalize_simulation()

    def _finalize_simulation(self):
        ...

    def time_advance(self):

        # TODO: Fully implement units into this.

        step = self.results._index
        T_e = self.results.T_e[step - 1].value
        n_e = self.results.n_e[step - 1].value
        dt = self._dt.value

        if self.verbose:
            print(
                f"step={step}  T_e={T_e}  n_e={n_e}  dt={dt}"
            )

        new_ionic_fractions = {}

        try:
            for elem in self.elements:
                nstates = self.results.nstates[elem]
                f0 = self.results._ionic_fractions[elem][:, self.results._index - 1]
                evals = self.EigenDataDict[elem].eigenvalues(T_e=T_e)
                evect = self.EigenDataDict[elem].eigenvectors(T_e=T_e)
                evect_inverse = self.EigenDataDict[elem].eigenvector_inverses(T_e=T_e)

                diagonal_evals = np.zeros((nstates, nstates), dtype=np.float64)
                for ii in range(0, nstates):
                    diagonal_evals[ii, ii] = np.exp(evals[ii] * dt * n_e)

                matrix_1 = np.dot(diagonal_evals, evect)
                matrix_2 = np.dot(evect_inverse, matrix_1)

                ft = np.dot(f0, matrix_2)

                # Can probably simplify this with np.where
                minconce = 1.0e-15
                for ii in np.arange(0, nstates, dtype=np.int):
                    if (abs(ft[ii]) <= minconce):
                        ft[ii] = 0.0

                new_ionic_fractions[elem] = ft

        except Exception as exc:
            raise NEIError(f"Unable to do time advance for {elem}") from exc
        else:

            new_time = self.results.time[0] + self._dt
            self.results._assign(
                new_time=new_time,
                new_ionic_fractions=new_ionic_fractions,
                new_T_e=self.electron_temperature(new_time),
                new_n=self.hydrogen_number_density(new_time),
            )


        # ------------------------------------------------------------------------------
        # function: Time-Advance solover
        # ------------------------------------------------------------------------------
        # def func_solver_eigenval(natom, te, ne, dt, f0, table):
        #    """
        #        The testing function for performing time_advance calculations.
        #    """

        # !! Change the following to use table.eigen*(T_e=...)

        #    table.temperature = te
        #    evals = table.eigenvalues  # find eigenvalues on the chosen Te node
        #    evect = table.eigenvectors
        #    evect_invers = table.eigenvector_inverses

        #    # define the temperary diagonal matrix
        #    diagona_evals = np.zeros((natom + 1, natom + 1))
        #    for ii in range(0, natom + 1):
        #        diagona_evals[ii, ii] = np.exp(evals[ii] * dt * ne)

        #    # matirx operation
        #    matrix_1 = np.dot(diagona_evals, evect)
        #    matrix_2 = np.dot(evect_invers, matrix_1)

        #    # get ions fraction at (time+dt)
        #    ft = np.dot(f0, matrix_2)

        #    # re-check the smallest value
        #    minconce = 1.0e-15
        #    for ii in np.arange(0, natom + 1, dtype=np.int):
        #        if (abs(ft[ii]) <= minconce):
        #            ft[ii] = 0.0
        #    return ft

    def save(self, filename="nei.h5"):
        ...