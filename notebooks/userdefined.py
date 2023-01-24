# custom wrapper functions to help set up and run chemical kinetic simulations using cantera
# authored by: Adam J. Susa (asusa@stanford.edu)

import cantera as ct
import numpy as np
import os
import matplotlib.pyplot as plt


T0 = 1200  # K
P0 = 101325  # Pa (1 atm)
phi0 = 1.
fuel0 = {'H2': 1.}
oxidizer0 = {'AR': .79, 'O2': .21}  # "airgon"


def setup_simulation(sln, T=T0, P=P0, phi=phi0, f=fuel0, ox=oxidizer0):
    """
    Returns a constant-pressure reactor object and reactor network based upon a ct.Solution object and at optionally
     specified starting conditions, or predefined defaults
    :param sln: ct.Solution object specifying mechanism for given simulation
    :param T: float, optional starting temperature (K)
    :param P: float, optional starting pressure (Pa)
    :param phi: float, optional equivalence ratio (-)
    :param f: dict, optional dict specifying fuel composition
    :param ox: dict, optional dict specifying oxidizer composition
    :return tuple: (ct.ConstPressureReactor, ct.ReactorNet)
    """
    sln.set_equivalence_ratio(phi=phi, fuel=f, oxidizer=ox)
    sln.TP = T, P
    r1 = ct.ConstPressureReactor(sln)
    rnet = ct.ReactorNet([r1])
    return r1, rnet


def new_result_dict(sln: ct.Solution, species_list=None):
    """
    Returns an empty result dictionary used for storing state variables during kinetic simulations
    :param sln: ct.Solution object to be used in the simulation
    :param species_list: list of species for which to record mole fractions; default None to record all species in
     mechansim
    :return dict: result dictionary with empty lists in fields for time, T, P, and all species, along with a list of
     all included species and another of the species indices within the sln.X array
    """
    # if the species list is None, make entries for all the species in the Solution object
    if species_list is None:
        species_list = sln.species_names
    else:  # remove any species in the list not in the Solution to prevent errors
        species_list_2, species_list = list(species_list), []  # store provided species_list as species_list_2, make a
        #   new species list for storing valid species
        for sp in species_list_2:
            if sp in sln.species_names and sp not in species_list:
                species_list += [sp]
            elif sp not in sln.species_names:
                print(f"{sp} not in mechanism; removed from tracked species.")

    # start initializing the dictionary with basic properties - time, temperature (T), pressure (P)
    result_dict = {'time': [],
                   'T': [],
                   'P': []}
    # now add entries for the species
    for sp in species_list:
        result_dict[sp] = []
    # finally, add entries with a list of the species and their indices
    result_dict['species'] = species_list
    result_dict['indices'] = [sln.species_index(sp) for sp in species_list]

    return result_dict


def add_to_result(sln, network, res_dict):
    """
    Adds the current simulation state to the results dictionary
    :param sln: ct.Solution at current simulation state
    :param network: ct.ReactorNet at current simulation state
    :param res_dict: formatted results dictionary, as created by "new_result_dict" function
    :return: None
    """
    res_dict['time'] += [network.time]
    res_dict['T'] += [sln.T]
    res_dict['P'] += [sln.P]
    for sp, x in zip(res_dict['species'], sln.X[res_dict['indices']]):
        res_dict[sp] += [x]


def plot_species_result(result_dict, species=None, logx=True, tmin=1e-6, tmax=None, ylim_x=(1e-6, 1), title=''):
    """

    :param result_dict:
    :param logx:
    :param tmin:
    :param tmax:
    :param ylim_x:
    :param title:
    :return:
    """
    # plot mole fractions vs. time
    plt.figure(figsize=(8, 6))

    if species is None:
        species = result_dict['species']

    ts = result_dict['time']
    for sp in species:
        if logx:
            plt.loglog(ts, result_dict[sp], label=sp)
        else:
            plt.semilogy(ts, result_dict[sp], label=sp)
    plt.legend()
    plt.xlabel('Time (s)')
    plt.ylabel('Mole Fractions (-)')
    plt.ylim(ylim_x)
    plt.xlim((tmin, tmax))
    plt.title(title)


def plot_simulation_result(result_dict, logx=True, tmin=1e-6, tmax=None, ylim_x=(1e-6, 1), title=''):
    """
    plots the simulation results stored in a result dictionary
    :param result_dict: formatted results dictionary, as created by "new_result_dict" and filled using "add_to_result"
    :param logx: bool, True (default) to plot time as log scale, False to plot linear time scale
    :param tmin: float, minimum time to plot (default 1e-6)
    :param tmax: float, maximum time to plot (default None to plot until last simulation time)
    :param ylim_x: tuple, (min, max) y limits of mole-fraction plot
    :param title: string, title to assign to plots
    :return: None
    """
    ts = result_dict['time']
    if tmax is None:
        tmax = ts[-1]

    # plot temperature vs. time
    plt.figure(figsize=(8, 4))
    if logx:
        plt.semilogx(ts, result_dict['T'])
    else:
        plt.plot(ts, result_dict['T'])
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')
    plt.xlim((tmin, tmax))
    plt.title(title)

    # plot mole fractions vs. time
    plt.figure(figsize=(8, 6))
    for sp in result_dict['species']:
        if logx:
            plt.loglog(ts, result_dict[sp], label=sp)
        else:
            plt.semilogy(ts, result_dict[sp], label=sp)
    plt.legend()
    plt.xlabel('Time (s)')
    plt.ylabel('Mole Fractions (-)')
    plt.ylim(ylim_x)
    plt.xlim((tmin, tmax))
    plt.title(title)


def idt_x_max(result_dict, x_key, save_to_dict=True, print_out=False):
    """
    Calculates the ignition delay time based on the maximum value of a tracked variable.
    Recommended for use on transient radical species
    :param result_dict: formatted results dictionary, as created by "new_result_dict" and filled using "add_to_result"
    :param x_key: str, dictionary key of tracked variable by which to define IDT
    :param save_to_dict: bool, default True to add calculated IDT to result_dict
    :param print_out: bool, True to print calculated IDT, default False to suppress output
    :return: None if save_to_dict, else float calculated IDT
    """
    x = result_dict[x_key]

    ix_max = np.argmax(x)  # flips the sign of dx_dt if keyword arg "sign" is negative
    t_ign_x = result_dict['time'][ix_max]

    if print_out:
        print(f'IDT_{x_key+"_max":10s} = {t_ign_x*1e6:4.1f} us')

    if save_to_dict:
        idt_dict = result_dict.get('idt', {})
        idt_dict[x_key] = t_ign_x
        result_dict['idt'] = idt_dict
    else:
        return t_ign_x


def idt_dx_dt(result_dict, x_key, sign=1, save_to_dict=True, print_out=False):
    """
    Calculates the ignition delay time based on the maximum rate of change of a tracked variable.
    Recommended for use on temperature, fuel, oxygen, or major products
    :param result_dict: formatted results dictionary, as created by "new_result_dict" and filled using "add_to_result"
    :param x_key: str, dictionary key of tracked variable by which to define IDT
    :param sign: numeric, positive (e.g. 1, default) to detect rising edge, negative (e.g. -1) to detect falling edge
    :param save_to_dict: bool, default True to add calculated IDT to result_dict
    :param print_out: bool, True to print calculated IDT, default False to suppress output
    :return: None if save_to_dict, else float calculated IDT
    """
    dx = np.diff(result_dict[x_key])
    dt = np.diff(result_dict['time'])

    dx_dt = dx / dt
    t_mean = np.array(result_dict['time'][:-1]) + (dt / 2)

    idx_max = np.argmax(np.sign(sign) * dx_dt)  # flips the sign of dx_dt if keyword arg "sign" is negative
    t_ign_dx = t_mean[idx_max]

    if print_out:
        print(f'IDT_d{x_key+"/dt":9s} = {t_ign_dx*1e6:4.1f} us')

    if save_to_dict:
        idt_dict = result_dict.get('idt', {})
        idt_dict[x_key] = t_ign_dx
        result_dict['idt'] = idt_dict
    else:
        return t_ign_dx


def calculate_idt(result_dict, max_value=('H',), max_rise=('T', 'H2O'), max_fall=None):
    """
    Calculates the IDT using numerous methods and adds them to result_dict
    :param result_dict: formatted results dictionary, as created by "new_result_dict" and filled using "add_to_result"
    :param max_value: iterable of keys for which to calculate IDT as maximum absolute value; e.g. intermediates)
    :param max_rise: iterable of keys for which to calculate IDT as maximum rate of increase; e.g. temperature, products
    :param max_fall: iterable of keys for which to calculate IDT as maximum rate of decrease; e.g. fuel, oxygen
    :return: None
    """
    if max_value is not None:
        for key in max_value:
            idt_x_max(result_dict, key)
    if max_rise is not None:
        for key in max_rise:
            idt_dx_dt(result_dict, key, sign=1)
    if max_fall is not None:
        for key in max_fall:
            idt_dx_dt(result_dict, key, sign=1)


def run_simulation(sln, T=T0, P=P0, phi=phi0, f=fuel0, ox=oxidizer0, species_of_interest=None, t_stop=1e-3):
    """
    Performs kinetic simulation at specified conditions and returns the full results dictionary
    :param sln: ct.Solution object specifying mechanism for given simulation
    :param T: float, optional starting temperature (K)
    :param P: float, optional starting pressure (Pa)
    :param phi: float, optional equivalence ratio (-)
    :param f: dict, optional dict specifying fuel composition
    :param ox: dict, optional dict specifying oxidizer composition
    :param species_of_interest: list of species to track in simulations, default None to track all species
    :param t_stop: float, maximum time to which to run the simulation
    :return: dict of IDTs calculated with various metrics
    """
    # setup the reactor, result objects
    r1, rnet = setup_simulation(sln, T=T, P=P, phi=phi, f=f, ox=ox)
    result_dict = new_result_dict(sln, species_list=species_of_interest)

    # run the simulation
    while rnet.time < t_stop:
        rnet.step()
        add_to_result(sln, rnet, result_dict)

    return result_dict


def simulate_idt(sln, T=T0, P=P0, phi=phi0, f=fuel0, ox=oxidizer0, species_of_interest=None, t_stop=1e-3):
    """
    Performs kinetic simulation at specified conditions and returns a dictionary of ignition delay times
    :param sln: ct.Solution object specifying mechanism for given simulation
    :param T: float, optional starting temperature (K)
    :param P: float, optional starting pressure (Pa)
    :param phi: float, optional equivalence ratio (-)
    :param f: dict, optional dict specifying fuel composition
    :param ox: dict, optional dict specifying oxidizer composition
    :param species_of_interest: list of species to track in simulations, default None to track all species
    :param t_stop: float, maximum time to which to run the simulation
    :return: dict of IDTs calculated with various metrics
    """
    result_dict = run_simulation(sln, T=T, P=P, phi=phi, f=f, ox=ox,
                                 species_of_interest=['H', 'H2O'] + list(f.keys()), t_stop=t_stop)
    calculate_idt(result_dict, max_fall=tuple(f.keys()), max_rise=['H2O'])

    return result_dict['idt']


def calculate_sensitivities(sln, T=T0, P=P0, phi=phi0, f=fuel0, ox=oxidizer0, t_stop=1e-3,
                            plot_result=True, n_sens_max=10, title=None, IDT_def='T'):
    """
    Performs brute-force sensitivity at specified conditions, with option to plot result
    :param sln: ct.Solution object specifying mechanism for given simulation
    :param T: float, optional starting temperature (K)
    :param P: float, optional starting pressure (Pa)
    :param phi: float, optional equivalence ratio (-)
    :param f: dict, optional dict specifying fuel composition
    :param ox: dict, optional dict specifying oxidizer composition
    :param t_stop: float, maximum time to which to run the simulation
    :param plot_result: bool, default True to plot sensitivities, False to suppress plotting
    :param n_sens_max: int, number of most sensitive reactions to include in plot, default 10
    :param title: str, plot title; default None for no title
    :param IDT_def: str, definition of IDT to use for calculating sensitivity; must be among those included in
     "calculate_idt" function (default 'T')
    :return: np.array, sensitivity of IDT to each reaction in mechanism defining sln
    """
    # calculate the base result
    sln.set_multiplier(1.)  # this makes sure all the multipliers are set to 1.
    IDT_base = simulate_idt(sln, T=T, P=P, phi=phi, f=f, ox=ox, t_stop=t_stop)[IDT_def]

    # perturb the mechanism and calculate the new IDT
    f_k = 0.05
    IDT_sens = np.zeros(sln.n_reactions, dtype=float)

    for i in range(sln.n_reactions):
        sln.set_multiplier(1)
        sln.set_multiplier((1+f_k), i)
        IDT_sens[i] = simulate_idt(sln, T=T, P=P, phi=phi, f=f, ox=ox, t_stop=IDT_base*3)[IDT_def]

    # calculate the IDT sensitivity
    dIDT_IDT = (IDT_sens - IDT_base) / IDT_base
    sensitivities = dIDT_IDT / f_k

    if plot_result:
        i_abs_sort = np.argsort(np.abs(sensitivities))[::-1]  # flip the order to sort from high to low
        n_max_sens = sensitivities[i_abs_sort[:n_sens_max]]
        numbered_eqns = [f'{eqn} (r{i})' for i, eqn in enumerate(sln.reaction_equations())]
        n_max_eqns = np.array(numbered_eqns)[i_abs_sort[:n_sens_max]]

        i_sort_nmax = np.argsort(n_max_sens)[::-1]
        n_sort_sens = n_max_sens[i_sort_nmax]
        n_sort_eqns = n_max_eqns[i_sort_nmax]

        plt.figure(figsize=(8, n_sens_max/2))
        plt.barh(range(len(n_sort_sens)), n_sort_sens)
        plt.yticks(range(len(n_sort_sens)), n_sort_eqns)
        plt.xlabel('IDT Sensitivity')
        plt.axvline(0, ls=':', c='k')

        if title is not None:
            plt.title(title)

    return sensitivities