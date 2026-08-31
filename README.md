# 1D-Quantum-Schrodinger-Solver
Finite Difference 1D Time-Independent Schrödinger Equation (TISE) solver in Python with numerical eigenvalue analysis and analytical comparisons.
# 1D Finite Difference Schrödinger Equation Solver

A Python computational physics package for discretizing and solving the 1D Time-Independent Schrödinger Equation (TISE) using finite difference matrix methods.

The solver computes numerical energy eigenvalues $E_n$, eigenstate wavefunctions $\psi_n(x)$, probability densities $|\psi_n(x)|^2$, and quantum mechanical expectation values ($\langle x \rangle$, $\Delta x$, $\Delta p$) across continuous potentials.

---

## Physical Formulation

The 1D Time-Independent Schrödinger Equation is given by:

$$-\frac{\hbar^2}{2m} \frac{d^2 \psi(x)}{dx^2} + V(x)\psi(x) = E \psi(x)$$

### 1. Matrix Discretization
Using a uniform 1D spatial grid with step size $\Delta x = x_{i+1} - x_i$, the second spatial derivative is discretized using a 3-point central finite difference stencil:

$$\frac{d^2 \psi}{dx^2} \approx \frac{\psi_{i+1} - 2\psi_i + \psi_{i-1}}{\Delta x^2}$$

This transforms the continuous differential operator into a symmetric tri-diagonal Hamiltonian matrix $\mathbf{H} = \mathbf{T} + \mathbf{V}$:

$$\mathbf{H}_{ij} = -\frac{\hbar^2}{2m \Delta x^2} \left( \delta_{i, j+1} - 2\delta_{i,j} + \delta_{i, j-1} \right) + V(x_i) \delta_{i,j}$$

### 2. Quantum Expectation Values & Uncertainty
* **Spatial Position Expectation Value:**
  $$\langle x \rangle = \int_{-L}^{L} x \, |\psi(x)|^2 \, dx, \quad \Delta x = \sqrt{\langle x^2 \rangle - \langle x \rangle^2}$$

* **Momentum Uncertainty (Integration by Parts):**
  $$\langle p^2 \rangle = \hbar^2 \int_{-L}^{L} \left| \frac{d\psi}{dx} \right|^2 dx, \quad \Delta p = \sqrt{\langle p^2 \rangle}$$

---

## Supported Potentials

1. **Harmonic Oscillator:**
   $$V(x) = \frac{1}{2} m \omega^2 x^2, \quad E_n^{\text{analytical}} = \hbar \omega \left( n + \frac{1}{2} \right)$$

2. **Infinite Square Well (Box of Width $2L$):**
   $$V(x) = 0 \quad \text{for } x \in [-L, L], \quad E_n^{\text{analytical}} = \frac{(n+1)^2 \pi^2 \hbar^2}{2 m (2L)^2}$$

---

## Requirements & Installation

Dependencies:
* `numpy`
* `scipy`
* `matplotlib`

```bash
pip install numpy scipy matplotlib

How to RunClone or download this repository.Run the main Python script:CS_code Git.py
Enter your parameters in the GUI interface:$L$: Half-length of the spatial domain$N$: Number of grid points (e.g., 400 or 500)Potential: 1 for Harmonic Oscillator, 2 for Infinite Square WellStates: Number of energy levels to calculateClick Run to view numerical vs. analytical energies, percent errors, expectation values, and wavefunctions/probability density plots.
