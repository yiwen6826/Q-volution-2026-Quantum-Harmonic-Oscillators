# Q-volution-Quantum-Harmonic-Oscillators

<img src="imgs/trackC-banner.jpg">

**Team name:** The Fifth Harmony Bloch Busters.

**Team members:** Hari, Julia, Prakriti, [Rodrigo](https://www.linkedin.com/in/rodrigo-segura-moreno/) and Yiwen.

This code is our submission regarding *Track C: Harmonic Oscillator* for the **Girls in Quantum** 2026 Hackathon, [*Q-Volution*](https://www.girlsinquantum.com/hackathon).

The task for this challenge is as follows:

> Implement a quantum algorithm to solve the harmonic oscillator equation, compute kinetic and potential energies, and analyze circuit efficiency. Built on [Classiq](https://www.classiq.io/).

---

## The Mission

Apply the algorithm from the paper to a simple, well-understood physical system: the **quantum harmonic oscillator**. Specifically, implement a quantum algorithm to solve:

$$ y + \omega^2 y = 0, \quad y(0)=1, \quad y'(0)=1, \quad \omega=1 \,.$$

Once implemented, use the resulting quantum state to evaluate the system's **kinetic and potential energies** as a function of time in the interval $[0,1]$.

Explore how varying algorithmic parameters — such as the bounds used in `inplace_prepare_state()` — affects the accuracy of these energy values. Finally, analyze how resource-efficient your implementation is by comparing circuit depth and width under different optimization settings.

---

## Deliverables

A notebook containing:
* A quantum program that solves the harmonic oscillator equation using the algorithm from the paper.
* Computed kinetic and potential energy values as functions of time, derived from the simulated output.
* An investigation of how parameter choices (e.g., amplitude bounds in state preparation) affect energy estimations.
* A graphical analysis of circuit depth and width under different optimization settings.