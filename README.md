# Q-volution-Quantum-Harmonic-Oscillators

![logo of Track C: Harmonic Oscillator](imgs/trackC-banner.jpg)

- **Team name:** The Fifth Harmony Bloch Busters.

- **Team members:**
    - Hari
    - [Julia](www.linkedin.com/in/julia-huynh-5117b8273)
    - Prakriti
    - [Rodrigo](https://www.linkedin.com/in/rodrigo-segura-moreno/)
    - Yiwen

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
* A graphical analysis of circuit depth and width under different optimization settings.

## 💡Repository Summary 

### **interactive_plots.html**
Visit this website to visualize how energies, positions, and velocities of the harmonic oscillator changes with different parameters (e.g. initial position y0, initial velocity vy0, and offset b1) through interactive plots! 

### **Interactive Plot Codes**
These are the codes associated with the interactive plots for the main html file.

### **Large Data Analysis Code for Varying Parameters** 
Code for the large dataset generated to analyze overall trends over varying parameters. 

### **Longer Time Code** 
Jupyter Notebook exploring time parameters for the harmonic oscillator. At high values of t, the expected kinetic and potential energies blow up. We realized that due to symmetry characteristics, we only need to simulate one-fourth of the time period of the oscillator to model energy overall. 

### **Optimal k (approximation order) Plots**
Jupyter Notebook allowing for arbitrary k parameters for the harmonic oscillator by usage of a padding function. Generally larger k's correspond to smaller error as the Taylor Expansion accuracy increases. 
  
 


