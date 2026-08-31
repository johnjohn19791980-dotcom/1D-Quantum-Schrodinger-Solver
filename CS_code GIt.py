iimport numpy as np
import scipy.constants
from scipy.linalg import eigh
import matplotlib.pyplot as plt

# ---------------- Physical Constants ----------------
hbar = scipy.constants.hbar
m_e = scipy.constants.m_e
eV_conversion = scipy.constants.eV
omega = 1e15

# ---------------- Grid Setup ----------------
def create_grid(L, N):
    x = np.linspace(-L, L, N)
    dx = x[1] - x[0]
    return x, dx

# ---------------- Potential Functions ----------------
def choose_potential(choice, x, m):
    if choice == 1:    # Harmonic Oscillator
        return 0.5 * m * (omega ** 2) * (x ** 2)
    elif choice == 2:  # Infinite Square Well (-L to L)
        return np.zeros_like(x)
    else:
        raise ValueError("Choice must be 1 (HO) or 2 (Box)")

# ---------------- Hamiltonian Construction ----------------
def build_hamiltonian(V, m, hbar, dx):
    N = len(V)
    T_factor = hbar**2 / (2.0 * m * dx**2)
    
    main_diag = 2.0 * np.ones(N)
    off_diag = -1.0 * np.ones(N - 1)
    
    K = T_factor * (np.diag(main_diag) + np.diag(off_diag, k=1) + np.diag(off_diag, k=-1))
    V_matrix = np.diag(V)
    
    return K + V_matrix

# ---------------- Eigen Solver ----------------
def get_eigenvalues(H, num_states):
    eigvals, eigvecs = eigh(H)
    return eigvals[:num_states], eigvecs[:, :num_states]

# ---------------- Analytical Solution ----------------
def analytical_energy(n, choice, L, m, hbar):
    if choice == 1:
        return hbar * omega * (n + 0.5)
    elif choice == 2:
        a = 2.0 * L  # Box width (from -L to L)
        return ((n + 1)**2 * (np.pi**2) * (hbar**2)) / (2.0 * m * (a**2))

# ---------------- Expectation Values ----------------
def normalize_wavefunction(psi, dx):
    return psi / np.sqrt(np.sum(np.abs(psi)**2) * dx)

def expectation_x(x, psi, dx):
    return np.sum(x * (np.abs(psi)**2)) * dx

def expectation_x2(x, psi, dx):
    return np.sum((x**2) * (np.abs(psi)**2)) * dx

def expectation_p2(psi, dx, hbar):
    dpsi_dx = np.gradient(psi, dx)
    return (hbar**2) * np.sum(np.abs(dpsi_dx)**2) * dx

# ---------------- Main Interactive Routine ----------------
def run_interactive_simulation():
    print("=== 1D Quantum Schrödinger Solver ===")
    
    # User Inputs via console prompt
    try:
        L_input = input("Enter domain half-width L in nanometers [default = 3.0]: ").strip()
        L = float(L_input) * 1e-9 if L_input else 3e-9

        N_input = input("Enter number of grid points N [default = 500]: ").strip()
        N = int(N_input) if N_input else 500

        pot_input = input("Choose potential (1 = Harmonic Oscillator, 2 = Infinite Well) [default = 1]: ").strip()
        choice = int(pot_input) if pot_input else 1

        states_input = input("Enter number of energy states to calculate [default = 4]: ").strip()
        num_states = int(states_input) if states_input else 4
        
    except ValueError as e:
        print("\nInvalid input entered. Please run again with numerical values.")
        return

    # Numerical execution
    x, dx = create_grid(L, N)
    V = choose_potential(choice, x, m_e)

    H = build_hamiltonian(V, m_e, hbar, dx)
    energies, eigvecs = get_eigenvalues(H, num_states)

    # Output Printout Table
    print("\n" + "="*80)
    print(f"{'n':>3} {'E_num (eV)':>12} {'E_ana (eV)':>12} {'%diff':>8} {'<x> (nm)':>10} {'Δx (nm)':>10} {'Δp (kg·m/s)':>14}")
    print("-" * 80)

    for n in range(num_states):
        psi = normalize_wavefunction(eigvecs[:, n], dx)

        E = energies[n]
        E_ana = analytical_energy(n, choice, L, m_e, hbar)
        pdiff = 100.0 * (E - E_ana) / E_ana

        exp_x = expectation_x(x, psi, dx)
        exp_x2 = expectation_x2(x, psi, dx)
        delta_x = np.sqrt(max(0.0, exp_x2 - exp_x**2))

        exp_p2 = expectation_p2(psi, dx, hbar)
        delta_p = np.sqrt(max(0.0, exp_p2))

        print(f"{n:3d} {E/eV_conversion:12.4f} {E_ana/eV_conversion:12.4f} {pdiff:8.2f}% {exp_x*1e9:10.3e} {delta_x*1e9:10.3e} {delta_p:14.3e}")

    print("="*80 + "\n")

    # Plotting Output
    plt.figure(figsize=(12, 5))
    
    # 1. Wavefunctions
    plt.subplot(1, 2, 1)
    for n in range(num_states):
        psi = normalize_wavefunction(eigvecs[:, n], dx)
        plt.plot(x * 1e9, psi, label=rf"$\psi_{{{n}}}(x)$")
    plt.xlabel("Position x (nm)")
    plt.ylabel(r"$\psi(x)$")
    plt.title("Eigenstate Wavefunctions")
    plt.grid(True)
    plt.legend()

    # 2. Probability Densities
    plt.subplot(1, 2, 2)
    for n in range(num_states):
        psi = normalize_wavefunction(eigvecs[:, n], dx)
        plt.plot(x * 1e9, np.abs(psi)**2, label=rf"$|\psi_{{{n}}}(x)|^2$")
    plt.xlabel("Position x (nm)")
    plt.ylabel(r"$|\psi(x)|^2$")
    plt.title("Probability Densities")
    plt.grid(True)
    plt.legend()
    
    plt.tight_layout()
    plt.show()

# Trigger interactive session
run_interactive_simulation()