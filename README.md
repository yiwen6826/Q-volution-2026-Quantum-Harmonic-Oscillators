# Q-volution-Quantum-Harmonic-Oscillators

![logo of Track C: Harmonic Oscillator](imgs/trackC-banner.jpg)

- **Team name:** The Fifth Harmony Bloch Busters.

- **Team members:**
    - [Hari](https://www.linkedin.com/in/hari-krishna-sahoo/)
    - [Julia](https://www.linkedin.com/in/julia-huynh-5117b8273)
    - [Prakriti](https://www.linkedin.com/in/prakriti-shahi/)
    - [Rodrigo](https://www.linkedin.com/in/rodrigo-segura-moreno/)
    - [Yiwen](https://www.linkedin.com/in/yiwen-lu-11a60431b/)

This code is our submission regarding *Track C: Harmonic Oscillator* for the **Girls in Quantum** 2026 Hackathon, [*Q-Volution*](https://www.girlsinquantum.com/hackathon).

The task for this challenge is as follows:

> Implement a quantum algorithm to solve the harmonic oscillator equation, compute kinetic and potential energies, and analyze circuit efficiency. Built on [Classiq](https://www.classiq.io/).

---

## 🎯 The Mission

Apply the algorithm from the paper [A Quantum Algorithm for Solving Linear Differential Equations: Theory and Experiment (arXiv:1807.04553)](https://arxiv.org/pdf/1807.04553) to a simple, well-understood physical system: the **quantum harmonic oscillator**. Specifically, implement a quantum algorithm to solve:

$$ y + \omega^2 y = 0, \quad y(0)=1, \quad y'(0)=1, \quad \omega=1 \,.$$

Once implemented, use the resulting quantum state to evaluate the system's **kinetic and potential energies** as a function of time in the interval $[0,1]$.

Explore how varying algorithmic parameters — such as the bounds used in `inplace_prepare_state()` — affects the accuracy of these energy values. Finally, analyze how resource-efficient your implementation is by comparing circuit depth and width under different optimization settings.

---

## 📋 Deliverables

A notebook containing:
* A quantum program that solves the harmonic oscillator equation using the algorithm from the paper.
* Computed kinetic and potential energy values as functions of time, derived from the simulated output.
* An investigation of how parameter choices (e.g., amplitude bounds in state preparation) affect energy estimations.
* ~~A graphical analysis of circuit depth and width under different optimization settings.~~

## 💡 Repository Summary 

- `index.html`
HTML displaying interactive plots to visualize how energies, positions, and velocities of the harmonic oscillator change with different parameters (initial position `y0`, initial velocity `vy0`, and offset to equilibrium position `b1`) through interactive plots! Visit the GitHub pages website [here](https://yiwen6826.github.io/Q-volution-2026-Quantum-Harmonic-Oscillators/). 

- `1_harmonic_oscillator.ipynb`
Jupyter Notebook explaining the differential equation to solve (harmonic oscillator) and the implementation of the quantum algorithm aforementioned.

- `2_interactive_plots.ipynb`
Jupter Notebook used for creating the interactive plots displayed at our [GitHub pages website](https://yiwen6826.github.io/Q-volution-2026-Quantum-Harmonic-Oscillators/).

- `3_time.ipynb`
Jupyter Notebook exploring time parameters for the harmonic oscillator. At high values of `t`, the expected kinetic and potential energies blow up. We realized that due to symmetry in the periodicity of the solution, we only need to simulate one-fourth of the time period of the oscillator to model energy overall.

- `4_sweep.ipynb`
Jupyter Notebook showcasing the behaviour of the kinetic and potential energies for `k = 3, 5, 9`, `y0 = 0, 1, 10` and `vy0 = 0, 1, 10`. As well as the mean accuracy as a function of the Taylor order expansion `k`.

- `5_run_cases.ipynb`
Jupter Notebook containing a three cell tutorial on how to run the quantum algorithm, save the data and load pre-existing data.

- `sweep_results.csv`
Data generated from `4_sweep.ipynb`

- `fifth_harmony/`
Custom module which handles the data plotting, saving and loading. Includes a modularized version of the code in `1_harmonic_oscillator.ipynb` via the `harmonic_oscillator` function to simplify the generation of datasets.

- `imgs/`
Directory with images used throughout the repository only for illustrative purposes.

- `interactive_plots/`
Directory containing the interactive plots generated with `2_interactive_plots.ipynb` and showcased the repository's [GitHub pages website](https://yiwen6826.github.io/Q-volution-2026-Quantum-Harmonic-Oscillators/).

- `k_plots/`
Directory containing plots showcasing `k` vs the error (mean, kinetic energy and potential energy).

  - `k_plots/optimal_k.ipynb`
Jupyter Notebook allowing for arbitrary `k` parameters for the harmonic oscillator by usage of a padding function. Generally larger values of `k` correspond to smaller error as the Taylor Expansion accuracy increases.
  
 ## How to use this code

1. Create a conda environment with Python v3.12.12
2. `pip install requirements.txt`
3. Enjoy :-)
