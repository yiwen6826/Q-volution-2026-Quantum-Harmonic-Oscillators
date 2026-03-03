# This code is part of 'The Fifth Harmony Bloch Busters' submission
# for "Track C: Harmonic Oscillator" challenge in Q-Volution 2026
#
# (C) Copyright The Fifth Harmony Bloch Busters, 2026.
#
# This code is licenced by the MIT Licence.
#
# DISCLAIMER: This Python code is a modification of Classiq's implementation [Ref1] of the paper
# A Quantum Algorithm for Solving Linear Differential Equations: Theory and Experiment (arXiv:1807.04553) [Ref2].
# Extended to handle nonunitary A (omega != 1) using the alternative decomposition in Appendix D.
#
# [Ref1] https://github.com/Classiq/classiq-library/blob/main/community/paper_implementation_project/quantum_algo_for_solving_linear_differential_equations/harmonic_oscillator.ipynb
# [Ref2] https://arxiv.org/pdf/1807.04553

from classiq import *
from scipy.linalg import sqrtm
from scipy.special import factorial
import numpy as np


authenticate(overwrite=True)


def unitary_harmonic_oscillator(k:int, omega:float, y0:float, vy0:float, b0:float, b1:float, msg:bool=True) -> tuple[np.array, dict, dict]:
    """ Solution for the harmonic oscillator using a modified version of [Ref1]
        which is based on arXiv:1807.04553.

    k     : order of the Taylor expansion, must be of the form k = 2**m - 1 (m a natural number)
    omega : frequency (omega != 1 makes A nonunitary)
    y0    : initial position
    vy0   : initial velocity
    b0    : shift to velocity
    b1    : shift to equilibrium position
    msg   : whether to print progress messages

    Returns:
        time           : time range used for computation
        data_expected  : dict of expected values (classical analytical solution)
        data_actual    : dict of actual values from quantum algorithm

        Both dicts have structure:
            {"y": ..., "dydt": ..., "KE": ..., "PE": ...}

    [Ref1] https://github.com/Classiq/classiq-library/blob/main/community/paper_implementation_project/quantum_algo_for_solving_linear_differential_equations/harmonic_oscillator.ipynb
    """

    # -----------------------------------------------------------------------
    # PROBLEM SETUP
    # -----------------------------------------------------------------------
    n = 2  # dimension of vectors x(0) and b
    M = np.array([[0, 1], [-omega**2, 0]])
    
    x0 = np.array([y0, vy0]) # initial vector (pos, vel)
    b = np.array([b0, b1]) # shift to (velocity, equilibrium position)
    
    # Norms (used in coefficients)
    x0_norm = np.linalg.norm(x0)
    b_norm = np.linalg.norm(b)
    
    # Note that quantum states |x0> = x0/||x0|| and |b> = b/||b|| are computed at 'encoding'
    
    # Matrices
    # np.linalg.norm(X,2) is the l2 norm of matrix X
    M_norm = np.linalg.norm(M,2)
    A = M/np.linalg.norm(M,2)

    # -----------------------------------------------------------------------
    # SHARED MUTABLE STATE
    # -----------------------------------------------------------------------
    # avoids global scope issues inside nested functions
    state = {
        "vs1": np.eye(k + 1, k + 1),
        "vs2": np.eye(k + 1, k + 1),
        "v": np.eye(2, 2),
        "c": 0.0,
        "d": 0.0,
        "N": 0.0
    }

    # -----------------------------------------------------------------------
    # VS1
    # -----------------------------------------------------------------------   
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
    
    
    # -----------------------------------------------------------------------
    # VS2
    # -----------------------------------------------------------------------    
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
    

    # -----------------------------------------------------------------------
    # VS
    # -----------------------------------------------------------------------   
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
    
    # -----------------------------------------------------------------------
    # DECODING — same structure as original
    # -----------------------------------------------------------------------
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
    
    # -----------------------------------------------------------------------
    # EXECUTION PREFERENCES — unchanged
    # -----------------------------------------------------------------------
    execution_preferences = ExecutionPreferences(
        num_shots=1,
        backend_preferences=ClassiqBackendPreferences(
            backend_name=ClassiqSimulatorBackendNames.SIMULATOR_STATEVECTOR
        ),
    )
    
    # -----------------------------------------------------------------------
    # RUN OVER TIME RANGE
    # -----------------------------------------------------------------------
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
    B = y0 - (b1 / omega**2)
    C = (vy0 + b0) / omega
    phi = omega*(time_range - t_initial)
    y_expected = B * np.cos(phi) + C * np.sin(phi) + (b1 / omega**2) # y(t)
    y_dash_expected = omega * C * np.cos(phi) - omega * B * np.sin(phi) # y'(t)
    
    kinetic_expected = np.square(y_dash_expected) / 2 # KE = (1/2) * (y')^2
    potential_expected = np.square(omega) * np.square(y_expected - (b1 / omega**2)) / 2 # PE = (1/2) * ω^2 * [y - (b_y/ω^2)]^2
    
    # Actual energies calculated from quantum algorithm:
    kinetic_actual = np.square(y_dash_actual + b0)/2
    potential_actual = omega**2 * np.square(y_actual - b1/omega**2) /2

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


def non_unitary_harmonic_oscillator(k: int, omega: float, y0: float, vy0: float, b0: float, b1: float, msg: bool = True) -> tuple[np.array, dict, dict]:
    """ Solution for the harmonic oscillator using a modified version of [Ref1]
        which is based on arXiv:1807.04553.
        Extended to handle nonunitary A (arbitrary omega) via Appendix D decomposition.

    k     : order of the Taylor expansion, must be of the form k = 2**m - 1 (m a natural number)
    omega : frequency (omega != 1 makes A nonunitary)
    y0    : initial position
    vy0   : initial velocity
    b0    : shift to velocity
    b1    : shift to equilibrium position
    msg   : whether to print progress messages

    Returns:
        time           : time range used for computation
        data_expected  : dict of expected values (classical analytical solution)
        data_actual    : dict of actual values from quantum algorithm

        Both dicts have structure:
            {"y": ..., "dydt": ..., "KE": ..., "PE": ...}

    [Ref1] https://github.com/Classiq/classiq-library/blob/main/community/paper_implementation_project/quantum_algo_for_solving_linear_differential_equations/harmonic_oscillator.ipynb
    """

    # -----------------------------------------------------------------------
    # PROBLEM SETUP
    # -----------------------------------------------------------------------
    n  = 2
    M  = np.array([[0, 1], [-omega**2, 0]], dtype=complex)
    x0 = np.array([y0, vy0], dtype=complex)
    b  = np.array([b0, b1],  dtype=complex)

    x0_norm = np.linalg.norm(x0)
    b_norm  = np.linalg.norm(b)

    M_norm = np.linalg.norm(M, 2)
    A      = M / M_norm

    # -----------------------------------------------------------------------
    # APPENDIX D — DECOMPOSE A INTO 4 UNITARIES F1..F4
    # A = 0.5 * (F1 + F2 + F3 + F4), Eq. (D2), (D4), (D6)
    # -----------------------------------------------------------------------
    I_mat = np.eye(n, dtype=complex)
    Bh    = 0.5  * (A + A.conj().T)         # Hermitian part
    Ch    = 0.5j * (A.conj().T - A)         # Anti-Hermitian part

    F1 = Bh + 1j * sqrtm(I_mat - Bh @ Bh)
    F2 = Bh - 1j * sqrtm(I_mat - Bh @ Bh)
    F3 = 1j * Ch - sqrtm(I_mat - Ch @ Ch)  # i absorbed into F3, eq. (D6)
    F4 = 1j * Ch + sqrtm(I_mat - Ch @ Ch)  # i absorbed into F4, eq. (D6)

    F_unitaries = [F1, F2, F3, F4]

    assert np.allclose(A, 0.5 * sum(F_unitaries), atol=1e-10), \
        "Decomposition check failed: A != 0.5*(F1+F2+F3+F4)"

    # -----------------------------------------------------------------------
    # ANCILLA DIMENSIONS — Appendix D
    #
    # |x0> branch has K   = (4^(k+1) - 1) / 3  terms
    # |b>  branch has K_b = (4^k     - 1) / 3  terms
    # Second ancilla size is K (the larger of the two).
    # Third ancilla needs 2k qubits to index up to 4^k combinations.
    # -----------------------------------------------------------------------
    K   = (4**(k + 1) - 1) // 3

    # After computing K, find the next power of 2
    K_padded = 2 ** (K).bit_length()   # next power of 2 >= K
    
    # Update T to match
    T  = K_padded.bit_length() - 1   # = log2(K_padded), guaranteed int --- second ancilla qubits
    T3  = 2 * k   # third ancilla qubits
    dim = 1       # log2(2) = 1, hardcode since n=2 always

    if msg:
        print(f"k={k}, K={K}, T={T}, T3={T3}, dim={dim}")

    # -----------------------------------------------------------------------
    # SHARED MUTABLE STATE
    # -----------------------------------------------------------------------
    state = {
        "vs1": np.eye(K_padded, K_padded, dtype=complex),
        "vs2": np.eye(K_padded, K_padded, dtype=complex),
        "v":   np.eye(2,  2, dtype=complex),
        "c":   0.0,
        "d":   0.0,
        "N":   0.0
    }

    # -----------------------------------------------------------------------
    # HOUSEHOLDER HELPER
    # Builds a unitary whose first column equals target_col.
    # -----------------------------------------------------------------------
    def householder_unitary(target_col, padded_size):
        """Builds unitary of size padded_size x padded_size 
           whose first column equals target_col (zero-padded)."""
        col = np.zeros(padded_size, dtype=complex)
        col[:len(target_col)] = target_col   # zero-pad the input column
    
        e    = np.zeros(padded_size, dtype=complex)
        e[0] = 1.0
        w    = col - e
        if np.linalg.norm(w) < 1e-14:
            return np.eye(padded_size, dtype=complex)
        w = w / np.linalg.norm(w)
        return np.eye(padded_size, dtype=complex) - 2.0 * np.outer(w, w.conj())

    # -----------------------------------------------------------------------
    # VS1 — Appendix D, Eq. (D9)
    # Non-zero entries at sparse indices idx = (4^j - 1)//3, j = 0..k
    # C_m[idx] = ||x0|| * (M_norm * t/2)^j / j!
    # -----------------------------------------------------------------------
    def VS1(t):
        c_m = np.zeros(K_padded, dtype=complex)   # K_padded instead of K
        for j in range(k + 1):
            idx      = (4**j - 1) // 3
            c_m[idx] = x0_norm * (M_norm * t / 2)**j / factorial(j)
    
        c          = np.sqrt(np.sum(np.abs(c_m)))
        state["c"] = c
    
        if t == 0 or c < 1e-14:
            state["vs1"] = np.eye(K_padded, dtype=complex)
            return
    
        first_col    = np.sqrt(np.abs(c_m)) / c
        state["vs1"] = householder_unitary(first_col, K_padded)

    # -----------------------------------------------------------------------
    # VS2 — Appendix D, Eq. (D9)
    # Non-zero entries at indices idx = (4^(j-1) - 1)//3, j = 1..k
    # D_n[idx] = ||b|| * (M_norm * t/2)^(j-1) * t / j!
    # Last entry forced to 0 per Eq. (5).
    # -----------------------------------------------------------------------
    def VS2(t):
        d_n = np.zeros(K_padded, dtype=complex)   # K_padded instead of K
        for j in range(1, k + 1):
            idx      = (4**(j - 1) - 1) // 3
            d_n[idx] = b_norm * (M_norm * t / 2)**(j-1) * t / factorial(j)
    
        d_n[K_padded - 1] = 0
        d          = np.sqrt(np.sum(np.abs(d_n)))
        state["d"] = d
    
        if d < 1e-14:
            state["vs2"] = np.eye(K_padded, dtype=complex)
            return
    
        first_col    = np.sqrt(np.abs(d_n)) / d
        state["vs2"] = householder_unitary(first_col, K_padded)

    # -----------------------------------------------------------------------
    # V — unchanged from original, Eq. (4)
    # -----------------------------------------------------------------------
    def V():
        c = state["c"]
        d = state["d"]
        N = np.sqrt(c * c + d * d)
        state["N"] = N
        if N == 0:
            state["v"] = np.eye(2, 2, dtype=complex)
        else:
            state["v"] = np.array([[c / N,  d / N],
                                   [d / N, -c / N]], dtype=complex)

    # -----------------------------------------------------------------------
    # ENCODING — same structure as original
    # -----------------------------------------------------------------------
    @qfunc
    def encoding(x: QNum, ancilla: QNum, y: QBit, t: float):
        prob_x0 = [float(np.real(i / x0_norm)) for i in x0] if x0_norm > 0 else [0.0] * len(x0)
        prob_b  = [float(np.real(i / b_norm))  for i in b]  if b_norm  > 0 else [0.0] * len(b)

        VS1(t)
        VS2(t)
        V()

        unitary(state["v"], y)

        control(
            y == 0,
            lambda: (inplace_prepare_amplitudes(prob_x0, 0.01, x), unitary(state["vs1"], ancilla)),
            lambda: (inplace_prepare_amplitudes(prob_b,  0.01, x), unitary(state["vs2"], ancilla))
        )

    # -----------------------------------------------------------------------
    # EVOLUTION — Appendix D
    #
    # For Taylor order j, there are 4^j combinations of F_unitaries.
    # Each combination sits at second-ancilla index: base = (4^j-1)//3 + combo_idx
    # The third ancilla (ancilla3) labels which combination (0..4^j - 1).
    # The product of F's is built by reading combo_idx in base-4.
    # -----------------------------------------------------------------------
    @qfunc
    def evolution(x: QNum, ancilla: QNum, ancilla3: QNum):
        for j in range(k + 1):
            base_idx = (4**j - 1) // 3
            n_combos = 4**j

            for combo_idx in range(n_combos):
                # Build F product for this combination
                U   = np.eye(n, dtype=complex)
                tmp = combo_idx
                for _ in range(j):
                    U   = F_unitaries[tmp % 4] @ U
                    tmp = tmp // 4

                ancilla_state = base_idx + combo_idx

                control(
                    (ancilla == ancilla_state) & (ancilla3 == combo_idx),
                    lambda U=U: unitary(U, x)
                )

    # -----------------------------------------------------------------------
    # DECODING — same structure as original
    # -----------------------------------------------------------------------
    @qfunc
    def decoding(ancilla: QNum, y: QBit):
        ws1 = state["vs1"].conj().T
        ws2 = state["vs2"].conj().T
        w   = state["v"].conj().T
        control(y == 0, lambda: unitary(ws1, ancilla), lambda: unitary(ws2, ancilla))
        unitary(w, y)

    # -----------------------------------------------------------------------
    # MAIN CIRCUIT FACTORY — third ancilla added
    # -----------------------------------------------------------------------
    def create_main_for_t(t: float):
        @qfunc
        def main(x:        Output[QNum[dim]],
                 ancilla:  Output[QNum[T]],
                 ancilla3: Output[QNum[T3]],
                 y:        Output[QBit]):
            allocate(x)
            allocate(ancilla)
            allocate(ancilla3)
            allocate(y)

            encoding(x, ancilla, y, t)
            evolution(x, ancilla, ancilla3)
            decoding(ancilla, y)

        return main

    # -----------------------------------------------------------------------
    # EXECUTION PREFERENCES — unchanged
    # -----------------------------------------------------------------------
    execution_preferences = ExecutionPreferences(
        num_shots=1,
        backend_preferences=ClassiqBackendPreferences(
            backend_name=ClassiqSimulatorBackendNames.SIMULATOR_STATEVECTOR
        ),
    )

    # -----------------------------------------------------------------------
    # RUN OVER TIME RANGE
    # -----------------------------------------------------------------------
    y_list      = []
    y_dash_list = []

    t_initial  = 0
    t_final    = 1
    t_step     = 11
    time_range = np.linspace(t_initial, t_final, t_step)

    for i, t in enumerate(time_range):
        if msg:
            print(f"Time = {t} --- case {i+1}/{t_step}")

        qmod    = create_model(create_main_for_t(t))
        qmod    = set_execution_preferences(qmod, execution_preferences)
        qprog   = synthesize(qmod)
        job     = execute(qprog)
        results = job.result_value()

        for j in results.parsed_state_vector:
            if int(j.bitstring[:-dim], 2) == 0:
                amplitude_scaled = np.linalg.norm(j.amplitude) * (state["N"] ** 2)
                if msg:
                    print(j.bitstring, " : ", amplitude_scaled)
                if int(j.bitstring, 2) == 0:
                    y_list.append(amplitude_scaled)
                if int(j.bitstring, 2) == 1:
                    y_dash_list.append(amplitude_scaled)

        if msg:
            print(20 * "-")

    y_actual      = np.array(y_list)
    y_dash_actual = np.array(y_dash_list)

    # -----------------------------------------------------------------------
    # ANALYTICAL SOLUTION
    # -----------------------------------------------------------------------
    B_const = y0  - (b1 / omega**2)
    C_const = (vy0 + b0) / omega
    phi     = omega * (time_range - t_initial)

    y_expected      = B_const * np.cos(phi) + C_const * np.sin(phi) + (b1 / omega**2)
    y_dash_expected = omega * C_const * np.cos(phi) - omega * B_const * np.sin(phi)

    kinetic_expected  = np.square(y_dash_expected) / 2
    potential_expected = np.square(omega) * np.square(y_expected - b1 / omega**2) / 2

    kinetic_actual  = np.square(y_dash_actual + b0) / 2
    potential_actual = omega**2 * np.square(y_actual - b1 / omega**2) / 2

    data_expected = {
        "y":    y_expected,
        "dydt": y_dash_expected,
        "KE":   kinetic_expected,
        "PE":   potential_expected
    }

    data_actual = {
        "y":    y_actual,
        "dydt": y_dash_actual,
        "KE":   kinetic_actual,
        "PE":   potential_actual
    }

    return time_range, data_expected, data_actual


def harmonic_oscillator(k:int, omega:float, y0:float, vy0:float, b0:float, b1:float, msg:bool=True):
    """ Compute the harmonic oscillator for the unitary and non-unitary case

    k     : order of the Taylor expansion, must be of the form k = 2**m - 1 (m a natural number)
    omega : frequency (omega != 1 makes A nonunitary)
    y0    : initial position
    vy0   : initial velocity
    b0    : shift to velocity
    b1    : shift to equilibrium position
    msg   : whether to print progress messages

    Returns:
        time           : time range used for computation
        data_expected  : dict of expected values (classical analytical solution)
        data_actual    : dict of actual values from quantum algorithm

        Both dicts have structure:
            {"y": ..., "dydt": ..., "KE": ..., "PE": ...}
    """
    M = np.array([[0, 1], [-omega**2, 0]])
    M_norm = np.linalg.norm(M,2)
    A = M/np.linalg.norm(M,2)

    # check i matrix is unitary
    if np.allclose(A.conj().T @ A, np.eye(A.shape[0])):
        if msg:
            print("--- Unitary case ---")
        time_range, data_expected, data_actual = unitary_harmonic_oscillator(k, omega, y0, vy0, b0, b1, msg=False)
    else:
        if msg:
            print("--- Non-unitary case ---")
        time_range, data_expected, data_actual = non_unitary_harmonic_oscillator(k, omega, y0, vy0, b0, b1, msg=False)

    return time_range, data_expected, data_actual