# This code is part of 'The Fifth Harmony Bloch Busters' submission
# for "Track C: Harmonic Oscillator" challenge in Q-Volution 2026
#
# (C) Copyright The Fifth Harmony Bloch Busters, 2026.
#
# This code is licenced by the MIT Licence.

import numpy as np
from tabulate import tabulate


def print_table(time_range:np.array, kinetic_expected:np.array, potential_expected:np.array, kinetic_actual:np.array, potential_actual:np.array) -> None:
    """ Tables of the expected vs actual kinetic and potential energy values

    t_final: final time
    kinetic_expected, potential_expected: arrays of classically computed energies
    kinetic_actual, potential_actual: arrays of energies obtained by quantum algorithm
    """
    results = [
        [
            t,
            kinetic_expected[i],
            kinetic_actual[i],
            potential_expected[i],
            potential_actual[i],
        ]
        for i, t in enumerate(time_range)
    ]
    table = tabulate(
        results,
        headers=[
            "t",
            "Kinetic expected",
            "Kinetic actual",
            "Potential expected",
            "Potential actual)",
        ],
        numalign="center",
        tablefmt="github",
    )
    print(table)
    
    error_bound_potential = [
        100 * np.abs(a - e) / a for a, e in zip(potential_expected, potential_actual)
    ]
    error_bound_kinetic = [
        100 * np.abs(a - e) / a for a, e in zip(kinetic_expected, kinetic_actual)
    ]
    
    print("\n\nkinetic accuracy (Mean): ", 100 - np.mean(error_bound_kinetic), "%")
    print("potential accuracy (Mean): ", 100 - np.mean(error_bound_potential), "%")