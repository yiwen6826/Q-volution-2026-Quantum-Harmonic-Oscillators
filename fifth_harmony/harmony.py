# This code is part of 'The Fifth Harmony Bloch Busters' submission
# for "Track C: Harmonic Oscillator" challenge in Q-Volution 2026
#
# (C) Copyright The Fifth Harmony Bloch Busters, 2026.
#
# This code is licenced by the MIT Licence.
#
# DISCLAIMER: This Python code is a modification of Classiq's implementation [Ref1] of the paper A Quantum Algorithm for Solving Linear Differential Equations: Theory and Experiment (arXiv:1807.04553) [Ref2].
#
# [Ref1] https://github.com/Classiq/classiq-library/blob/main/community/paper_implementation_project/quantum_algo_for_solving_linear_differential_equations/harmonic_oscillator.ipynb
# [Ref2] https://arxiv.org/pdf/1807.04553

from classiq import *
import numpy as np

authenticate()


def harmonic_oscillator(k:int, omega:float, y0:float, vy0:float, bx:float, by:float, msg:bool=True) -> tuple[np.array, dict, dict]:
    """ Solution for the harmonic oscillator using a modified version of [Ref1]
        which is based on arXiv:1807.04553

    k : order of the Taylor expansion, must be of the form k = 2**m - 1 (m a natural number)
    omega : frequency
    y0 : initial position
    vy0 : initial velocity
    bx : shift to velocity
    by : shift to equilibrium position
    
    msg : to print or not to print messages

    Returns :
        time : time range used for computation
        data_expected : dict of expected values of harmonic oscillator (classical computation)
        data_actual : dict of actual computed values from quantum algorithm
        
        both dictionaries have the following structure
        data_x = {"y":data1, "dydt":data2, "KE":data3, "PE":data4}

    [Ref1] https://github.com/Classiq/classiq-library/blob/main/community/paper_implementation_project/quantum_algo_for_solving_linear_differential_equations/harmonic_oscillator.ipynb
    """

    n = 2  # dimension of vectors x(0) and b
    M = np.array([[0, 1], [-omega**2, 0]])
    
    x0 = np.array([y0, vy0]) # initial vector (pos, vel)
    b = np.array([bx, by]) # shift to (velocity, equilibrium position)
    
    # Norms (used in coefficients)
    x0_norm = np.linalg.norm(x0)
    b_norm = np.linalg.norm(b)
    
    # Note that quantum states |x0> = x0/||x0|| and |b> = b/||b|| are computed at 'encoding'
    
    # Matrices
    # np.linalg.norm(X,2) is the l2 norm of matrix X
    M_norm = np.linalg.norm(M,2)
    A = M/np.linalg.norm(M,2)

    # Shared mutable state (avoids global scope issues inside nested functions)
    state = {
        "vs1": np.eye(k + 1, k + 1),
        "vs2": np.eye(k + 1, k + 1),
        "v": np.eye(2, 2),
        "c": 0.0,
        "d": 0.0,
        "N": 0.0
    }

    # --- VS1 ---    
    def VS1(t):    
        # --- Coefficient C_m, eq. (3) ---
        c = 0
        c_m: np.array = np.zeros(k + 1)
        m_factorial = 1
        for i in range(k + 1):
            c_m[i] = (x0_norm * (pow(t * M_norm, i))) / m_factorial
            c += c_m[i]
            m_factorial *= i + 1
    
        c = np.sqrt(c) # C = sqrt(sum C_m)
        state["c"] = c
    
        # --- Equation (5) ---
        if t == 0:
            state["vs1"] = np.eye(k + 1, k + 1)
            return
    
        # Construct Unitary matrix with the first column as defined above in the markdown
        e = np.zeros(k + 1)
        e[0] = 1
        w = (np.sqrt(c_m) / c) - e
        state["vs1"] = np.identity(k + 1) - np.multiply(2 * (1 / np.inner(w, w)), np.outer(w, w))
    
    
    # --- VS2 ---    
    def VS2(t):
        # --- Coefficient D_m, eq. (3) ---
        d = 0
        d_n: np.array = np.zeros(k + 1)
        n_factorial = 1
        for i in range(1, k + 1):
            d_n[i - 1] = (b_norm * (pow(M_norm * t, i - 1)) * t) / n_factorial
            d += d_n[i - 1]
            n_factorial *= i + 1
        
        d_n[k] = 0 # eq. (5)
        d = np.sqrt(d)  # D = sqrt(sum D_n)
        state["d"] = d
    
        # --- Equation (5) ---
        if d == 0:
            vs2 = np.eye(k + 1, k + 1)
            return
        
        # Construct Unitary matrix with the first column as defined above in the markdown
        e = np.zeros(k + 1)
        e[0] = 1
        w = (np.sqrt(d_n) / d) - e
        state["vs2"] = np.identity(k + 1) - np.multiply(2 * (1 / np.inner(w, w)), np.outer(w, w))
    

    # --- VS ---    
    def V():
        c = state["c"]
        d = state["d"]
        N = np.sqrt(c * c + d * d) # fancy N in the paper (see eq3)
        state["N"] = N
        if N == 0:
            state["v"] = np.eye(2, 2)
        else:
            state["v"] = np.array([[c / N, d / N], [d / N, -c / N]])
    
    
    @qfunc
    def encoding(x: QNum, ancilla: QNum, y: QBit, t: float):
        prob_x0 = [i / x0_norm for i in x0] if x0_norm > 0 else [0.0] * len(x0)
        prob_b  = [i / b_norm for i in b] if b_norm > 0 else [0.0] * len(b)
        
        # inplace_prepare_amplitudes(prob_x0, 0.01, x)
    
        VS1(t)
        VS2(t)
        V()
    
        unitary(state["v"], y)
        
        control(
            y == 0,
            lambda: (inplace_prepare_amplitudes(prob_x0, 0.01, x), unitary(state["vs1"], ancilla)),
            lambda: (inplace_prepare_amplitudes(prob_b,  0.01, x), unitary(state["vs2"], ancilla))
        )
    
    
    @qfunc
    def evolution(x: QNum, ancilla: QNum):
        u_m = np.array([[1, 0], [0, 1]])
    
        for i in range(k + 1):
            U = u_m.copy()
            control(ancilla == i, lambda U=U: unitary(U, x))
            u_m = u_m @ A
    
    
    @qfunc
    def decoding(ancilla: QNum, y: QBit):
        ws1 = state["vs1"].conj().T
        ws2 = state["vs2"].conj().T
        w = state["v"].conj().T
        control(y == 0, lambda: unitary(ws1, ancilla), lambda: unitary(ws2, ancilla))
        unitary(w, y)
    
    
    T = int(np.log2(k + 1))  # no. of ancilla qubits
    dim = int(np.log2(n))  # no. of work qubits
    
    def create_main_for_t(t: float):
        @qfunc
        def main(x: Output[QNum[dim]], ancilla: Output[QNum[T]], y: Output[QBit]):
            # ancilla: Output[QNum[T]] -> log2(k+1)
            allocate(x)
            allocate(ancilla)
            allocate(y)
    
            encoding(x, ancilla, y, t)
            evolution(x, ancilla)
            decoding(ancilla, y)
    
        return main
    
    
    execution_preferences = ExecutionPreferences(
        num_shots=1,
        backend_preferences=ClassiqBackendPreferences(
            backend_name=ClassiqSimulatorBackendNames.SIMULATOR_STATEVECTOR
        ),
    )
    
    
    y = [] # y(t)
    y_dash = [] # y'(t)
    
    # To standarize our results, the initial and final time, as well as the time step will be the below fixed values
    t_initial = 0 # initial time
    t_final = 1 # final time
    t_step = 11 # number of steps to divide [t_initial, t_final]
    time_range = np.linspace(t_initial, t_final, t_step)
    
    for i, t in enumerate(time_range):
        if msg:
            print(f"Time = {t} --- case {i+1}/{t_step}")
        qmod = create_model(create_main_for_t(t))
        qmod = set_execution_preferences(qmod, execution_preferences)
        qprog = synthesize(qmod)
        job = execute(qprog)
        results = job.result_value()
        for j in results.parsed_state_vector:
            if int(j.bitstring[:-dim], 2) == 0:
                if msg:
                    print(j.bitstring, " : ", np.linalg.norm(j.amplitude) * (state["N"] * state["N"]))
                if int(j.bitstring, 2) == 0:
                    y.append(np.linalg.norm(j.amplitude) * (state["N"] * state["N"]))
                if int(j.bitstring, 2) == 1:
                    y_dash.append(np.linalg.norm(j.amplitude) * (state["N"] * state["N"]))
        if msg:
            print(20*"-")

    y_actual = np.array(y)
    y_dash_actual = np.array(y_dash)

    # Expected energies calculated from the formulas
    B = y0 - (by / omega**2)
    C = (vy0 + bx) / omega
    phi = omega*(time_range - t_initial)
    y_expected = B * np.cos(phi) + C * np.sin(phi) + (by / omega**2) # y(t)
    y_dash_expected = omega * C * np.cos(phi) - omega * B * np.sin(phi) # y'(t)
    
    kinetic_expected = np.square(y_dash_expected) / 2 # KE = (1/2) * (y')^2
    potential_expected = np.square(omega) * np.square(y_expected - (by / omega**2)) / 2 # PE = (1/2) * ω^2 * [y - (b_y/ω^2)]^2
    
    # Actual energies calculated from quantum algorithm:
    kinetic_actual = np.square(y_dash_actual + bx)/2
    potential_actual = omega**2 * np.square(y_actual - by/omega**2) /2

    data_expected = {
        "y" : y_expected,
        "dydt" : y_dash_expected,
        "KE": kinetic_expected,
        "PE": potential_expected
    }
    
    data_actual = {
        "y" : y_actual,
        "dydt" : y_dash_actual,
        "KE": kinetic_actual,
        "PE": potential_actual
    }

    return time_range, data_expected, data_actual