# Quantum Circuits for the Advection–Diffusion Equation with Boundary Conditions Based on LCHS

Numerical-experiment code for the paper:

> **Quantum circuits for the advection-diffusion equation with boundary conditions based on LCHS**
> arXiv:[2605.17542](https://arxiv.org/abs/2605.17542)

This repository contains the simulation code for the numerical experiments reported in the
paper. Each case is a self-contained C++ program built on the
[QPanda2](https://github.com/OriginQ/QPanda-2) framework and runs on a full state-vector
emulator of a fault-tolerant machine.

---

## Method overview

The solvers implement the framework of the paper end to end:

- **Space discretisation** — the finite volume method (FVM). Cell averages are stored on a
  grid of `N = 2^n` cells per dimension, and the flux at each interface is reconstructed by a
  chosen scheme.
- **Flux construction schemes** — three schemes are implemented, giving the `central` /
  `upwind` / `exp` variants of each equation: a central flux, an upwind flux, and an
  exponential-fitting flux.
- **Boundary conditions** — Robin boundary conditions are handled by the newly constructed
  endpoint-projector circuits (`exp_siama11_gate`, `exp_siama00_gate`), with Dirichlet and
  Neumann obtained as special cases. Periodic boundaries are imposed through the wrap-around
  controlled rotations in `V_gate`.
- **LCHS transformation** — the non-unitary semigroup is represented as a linear combination
  of Hamiltonian simulations. The kernel `kernal_lgh` is evaluated on `M = 2^m` quadrature
  points of the truncated integration interval.
- **Outer LCU** — all physical and boundary-induced inhomogeneous terms are folded into an
  outer linear combination of unitaries, prepared on a separate ancilla register and applied
  through `Select_U_gate`.
- **Time integration** — a second-order Trotter–Suzuki decomposition with `r` Trotter steps of
  length `tau = T / r`.
- **Read-out** — the state is sampled with `probRunDict`, and the physical amplitude is
  recovered from the measured probability by rescaling with the normalisation constants of the
  initial state and of the prepared coefficient vectors.

### Qubit register naming

The code follows the notation of the paper:

| Symbol | Meaning |
| --- | --- |
| `n` / `N` | qubits / number of spatial discretisation points (`N = 2^n`) |
| `m` / `M` | qubits / number of LCHS coefficients (`M = 2^m`) |
| `l` / `L` | qubits / number of outer LCU coefficients (`L = 2^l`) |
| `r` | number of Trotter steps |

---

## Repository layout

```
.
├── 1diff_homo_dirichlet_central/
│   ├── diff_homo_dirichlet_central.cpp
│   ├── <parameters>.plt
│   └── <parameters>error.dat
├── ...
└── 9adv&diff_homo_periodic_exp_2d/
    ├── adv&diff_homo_periodic_exp_2d.cpp
    ├── <parameters>.plt
    └── <parameters>error.dat
```

Every directory holds one `.cpp` solver together with the output files it produced. Directory
names read as `<equation>_<homogeneity>_<boundary condition>_<flux scheme>[_2d]`:

- `diff` — diffusion equation, `adv` — advection equation, `adv&diff` — advection–diffusion equation
- `homo` — homogeneous problem, `inhomo` — problem with an inhomogeneous (source) term
- `dirichlet` / `neumann` / `robin` / `periodic` — boundary condition
- `central` / `upwind` / `exp` — flux construction scheme
- `2d` — two-dimensional problem

The `&` in directory names is part of the name and is interpreted as a literal character by git
and by GitHub; browsers percent-encode it in URLs.

---

## Case index

Each case is validated against the analytic solution of the corresponding model problem.

| Case | Directory | Equation | Boundary condition | Flux scheme | Dim |
| --- | --- | --- | --- | --- | --- |
| 1 | `1diff_homo_dirichlet_central` | diffusion | Dirichlet | central | 1D |
| 2 | `2diff_homo_neumann_central` | diffusion | Neumann | central | 1D |
| 3 | `3adv_homo_periodic_central` | advection | periodic | central | 1D |
| 3 | `3adv_homo_periodic_upwind` | advection | periodic | upwind | 1D |
| 4 | `4adv&diff_homo_dirichlet_central` | advection–diffusion | Dirichlet | central | 1D |
| 4 | `4adv&diff_homo_dirichlet_exp` | advection–diffusion | Dirichlet | exponential | 1D |
| 5 | `5adv&diff_homo_periodic_central` | advection–diffusion | periodic | central | 1D |
| 5 | `5adv&diff_homo_periodic_exp` | advection–diffusion | periodic | exponential | 1D |
| 6 | `6adv_inhomo_periodic_central` | advection with source | periodic | central | 1D |
| 6 | `6adv_inhomo_periodic_upwind` | advection with source | periodic | upwind | 1D |
| 7 | `7adv&diff_inhomo_dirichlet_central` | advection–diffusion with source | Dirichlet | central | 1D |
| 7 | `7adv&diff_inhomo_dirichlet_exp` | advection–diffusion with source | Dirichlet | exponential | 1D |
| 8 | `8diff_homo_robin_central_2d` | diffusion | Robin | central | 2D |
| 9 | `9adv&diff_homo_periodic_central_2d` | advection–diffusion | periodic | central | 2D |
| 9 | `9adv&diff_homo_periodic_exp_2d` | advection–diffusion | periodic | exponential | 2D |

The solvers allocate only the registers they need: the inhomogeneous cases carry an outer LCU
register in addition to the spatial one, and some of them an LCHS register as well, while the
purely advective cases construct the propagator directly.

---

## Building and running

These solvers link against **QPanda 2**, which is not part of this repository. Install it first
by following the official guide:

- QPanda-2 source and build instructions — <https://github.com/OriginQ/QPanda-2>
- Official *download, compile, install and use* tutorial — <https://qcloud.originqc.com.cn/document/qpanda-doc-cn/Tutorial.html>

The official guide covers the prerequisites (a C++14 compiler with OpenMP, CMake, `pkg-config`
and a UUID library), the CMake configuration options, the Windows and Linux/macOS build
sequences, and the default installation directories (`/usr/local` on Linux and macOS,
`C:/QPanda2` on Windows). Please refer to it for anything version-specific.

### Compiling a case

Point a CMake project at your QPanda installation, using the pattern from the official guide:

```cmake
cmake_minimum_required(VERSION 3.1)
project(solver)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${QPANDA_INSTALL_DIR}/lib/cmake")

find_package(QPANDA REQUIRED)
include_directories(${QPANDA_INCLUDE_DIR})

add_executable(solver "4adv&diff_homo_dirichlet_central/adv&diff_homo_dirichlet_central.cpp")
target_link_libraries(solver ${QPANDA_LIBRARIES})
```

```bash
cmake -B build -DQPANDA_INSTALL_DIR=/usr/local    # C:/QPanda2 on Windows
cmake --build build
```

The `Template` directory of the QPanda-2 repository contains further examples that link against
the installed library.

### Running

Each solver writes to a path relative to the working directory — for example
`./4adv&diff_homo_dirichlet_central/` — so **run the binary with the repository root as the
working directory**, and make sure the target directory exists: the solver does not create it.

All physical and numerical parameters are declared near the top of `main()` and are meant to be
edited in place; the 2D solvers take a per-dimension vector of spatial qubit counts. Output file
names encode the parameters of the run that produced them, so a directory may contain results
from several parameter sweeps alongside the current source defaults. The 2D cases and the
larger grid sizes are correspondingly more expensive to emulate.

---

## Output

Each run writes two files into its own directory.

**`<parameters>.plt`** — Tecplot point data, one row per cell centre, with the coordinates, the
exact solution, the LCHS solution, and the pointwise absolute error.

**`<parameters>error.dat`** — scalar diagnostics, one per line: the relative and absolute errors
in the `L-infinity`, `L1` and `L2` norms, the normalisation constants used to rescale the
measured probabilities, the post-selection success rate and its square root, and the wall-clock
runtime.

---

## Notes

- The emulator is QPanda's state-vector CPU backend. Memory and runtime scale exponentially in
  the qubit count, so the reported runtimes are emulator timings, not estimates of
  fault-tolerant quantum cost.
- Every solver is standalone: the gate helpers (`U_gate`, `W_gate`, `V_gate`,
  `exp_siama11_gate`, `exp_siama00_gate`) are duplicated in each file rather than shared, so a
  case can be read and modified in isolation.

## Citation

```bibtex
@article{advectiondiffusionLCHS,
  title  = {Quantum circuits for the advection-diffusion equation with
            boundary conditions based on LCHS},
  year   = {2026},
  eprint = {2605.17542},
  archivePrefix = {arXiv},
  primaryClass  = {quant-ph}
}
```

---

## License

Released under the MIT License. See [LICENSE](LICENSE).
