# This code is part of 'The Fifth Harmony Bloch Busters' submission
# for "Track C: Harmonic Oscillator" challenge in Q-Volution 2026
#
# (C) Copyright The Fifth Harmony Bloch Busters, 2026.
#
# This code is licenced by the MIT Licence.

import matplotlib.pyplot as plt
import numpy as np


def plot_expected_vs_actual(t_values:np.array, kinetic_expected:np.array, potential_expected:np.array, kinetic_actual:np.array, potential_actual:np.array, params:dict=None) -> None:
    """ Plots the expected vs actual kinetic and potential energy values

    t_values: time range
    kinetic_expected, potential_expected: arrays of classically computed energies
    kinetic_actual, potential_actual: arrays of energies obtained by quantum algorithm

    params : dictionary containing certain plotting parameters
        params = {'xmin':xmin, 'xmax':xmax, 'ymin':ymin, 'ymax':ymax, 'grid':grid}
    """
    
    plt.figure(figsize=(10, 6))

    # --- Plot data ---
    # Expected values (quantum-algo computed values)
    plt.plot(t_values, kinetic_expected, label="Expected Kinetic Energy", color="blue", linestyle="-")
    plt.plot(t_values, potential_expected, label="Expected Potential Energy", color="green", linestyle="-")    
    # Actual values (classical computed values)
    plt.scatter(t_values, kinetic_actual, label="Actual Kinetic Energy", color="red", marker="o")
    plt.scatter(t_values, potential_actual, label="Actual Potential Energy", color="orange", marker="o")
    
    # --- Plot config ---
    # ranges
    if params is not None:
        try:
            plt.xlim(params['xmin'], params['xmax'])
            plt.ylim(params['ymin'], params['ymax'])
        except KeyError: # key does not exist
            pass
        except NoneType: # dictionary does not exist
            pass
    # text
    plt.xlabel("Time t")
    plt.ylabel("Energy")
    plt.title("Energy Comparison: Expected vs. Actual")
    # aesthetics
    if params is not None:
        try:
            plt.grid(params['grid'])
        except KeyError:
            plt.grid(True)
        except NoneType:
            plt.grid(True)
    
    plt.legend()
    plt.tight_layout()
    
    plt.show()