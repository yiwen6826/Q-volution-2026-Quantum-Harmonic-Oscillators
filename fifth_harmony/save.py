# This code is part of 'The Fifth Harmony Bloch Busters' submission
# for "Track C: Harmonic Oscillator" challenge in Q-Volution 2026
#
# (C) Copyright The Fifth Harmony Bloch Busters, 2026.
#
# This code is licenced by the MIT Licence.

import numpy as np

import os


def make_dir(k:int, omega:float, y0:float, vy0:float, b0:float, b1:float) -> str:
    """ Name of directory where data will be saved or loaded from
    
    k : order of the Taylor expansion
    omega : frecuency of the oscillator
    y0 : initial position
    vy0 : initial velocity
    b0 : shift to velocity
    b1 : shift to equilibrium position
    """
    # convention is as follows: k_w_y0_vy0_b0_b1
    # `k` has no decimal places, all other values are rounded to 1 decimal places
    stem = f"{k}_{omega:.1f}_{y0:.1f}_{vy0:.1f}_{b0:.1f}_{b1:.1f}"
    
    return os.path.join("data", stem)


def data_exists(k:int, omega:float, y0:float, vy0:float, b0:float, b1:float) -> bool:
    """ Returns True if all expected and actual data files exist, False otherwise

    k : order of the Taylor expansion
    omega : frecuency of the oscillator
    y0 : initial position
    vy0 : initial velocity
    b0 : shift to velocity
    b1 : shift to equilibrium position
    """
    d = make_dir(k, omega, y0, vy0, b0, b1)
    files = [
        "time.npy",
        "expected_y.npy", "expected_dydt.npy", "expected_KE.npy", "expected_PE.npy",
        "actual_y.npy",   "actual_dydt.npy",   "actual_KE.npy",   "actual_PE.npy",
    ]
    return all(os.path.isfile(os.path.join(d, f)) for f in files)


def save_data(k:int, omega:float, y0:float, vy0:float, b0:float, b1:float, time:np.array, data_expected:dict, data_actual:dict) -> None:
    """ Saves the data for the expected y and dy/dt

    k : order of the Taylor expansion
    omega : frecuency of the oscillator
    y0 : initial position
    vy0 : initial velocity
    b0 : shift to velocity
    b1 : shift to equilibrium position

    time : time data
    data_i : dictionary of the form
        data_i = {"y":data1, "dydt":data2, "KE":data3, "PE":data4}
    """
    d = make_dir(k, omega, y0, vy0, b0, b1)
    os.makedirs(d, exist_ok=True)
    # time values
    np.save(f"{d}/time.npy",  time)
    # expected values
    np.save(f"{d}/expected_y.npy", data_expected["y"])
    np.save(f"{d}/expected_dydt.npy", data_expected["dydt"])
    np.save(f"{d}/expected_KE.npy", data_expected["KE"])
    np.save(f"{d}/expected_PE.npy", data_expected["PE"])
    # actual computed values from quantum algorithms
    np.save(f"{d}/actual_y.npy", data_actual["y"])
    np.save(f"{d}/actual_dydt.npy", data_actual["dydt"])
    np.save(f"{d}/actual_KE.npy", data_actual["KE"])
    np.save(f"{d}/actual_PE.npy", data_actual["PE"])

    print(f"All data has been saved at {d}/")


def load_data(k:int, omega:float, y0:float, vy0:float, b0:float, b1:float, data:str, msg:bool=True) -> tuple[np.array, np.array, np.array]:
    """ Loads the data for the expected y and dy/dt

    k : order of the Taylor expansion
    omega : frecuency of the oscillator
    y0 : initial position
    vy0 : initial velocity
    b0 : shift to velocity
    b1 : shift to equilibrium position

    data : string indicating the data to be loaded
        possible values are 'y', 'dydt', 'KE', 'PE'
    """
    d = make_dir(k, omega, y0, vy0, b0, b1)
    time = np.load(f"{d}/time.npy")
    y_expected = np.load(f"{d}/expected_{data}.npy")
    y_actual = np.load(f"{d}/actual_{data}.npy")

    if msg:
        print(f"Data from {d}/ has been loaded")
    
    return time, y_expected, y_actual