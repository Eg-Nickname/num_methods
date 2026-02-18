# Numerical methods | Jagiellonian university

**winter semester 2025/2026**

## Tech

- **Language** C++23
- **Linear Algebra Library** Eigen 5
- **Visualisations** Python (Pandas & Numpy & Matplotlib)
- **Typesetting** Typst (Reports)

## Project Structure

Codebase is documented in english but the formal reports are written in Polish. Each task (`N0X/`) follows a layout:

```text
.
├── src/                # C++ source code
├── scripts/            # Python scripts for data processing and chart generation
├── data/               # Raw output files (CSV) from simulations
├── figures/            # Generated plots
├── N0X_raport.pdf      # Compiled report (in Polish)
└── N0X_raport.typ      # Results analysis and theoretical background typestetting (in Polish)

```

## How to run?

In order to build specific task (N0X), navigate to its directory and use provided instructions in its README.md .

# Tasks

[N01] Eigen decompositions and Thomas Algorithm

Comparison and benchmarks of different methods for solving systems of equations with general Eigen methods and tailored solution for thridiagonal matrix. Program generates data for different matrix sizes and methods.

![Avg Speed Plot](./N01/figures/all_combined_10k.jpg)

[N02] Eigen decompositions and Thomas Algorithm with Sherman-Morrison Formula

Comparison and benchmarks of different methods for solving systems of equations with general Eigen methods and tailored solution for thridiagonal matrix and Sherman-Morrison formula. Program generates data for different matrix sizes and methods.

![Avg Speed Plot](./N02/figures/all_combined_10k.jpg)

[N03] Eigenvalue Problem

Program for finding biggest and smallest eigenvalues of a matrix. Comparison between Power method and Rayleigh method.

![Convergence Plot](./N03/figures/biggest_ev_cmp_plot.jpg)

[N04] Polynomial Interpolation & Runge Phenomenon

Program for finding lagrange interpolation of a given function. Shows comparison between equidistant nodes and Chebyshev nodes. Compares direct evaluation with evaluation from polynomial factors.

![Comparison Plot](./N04/figures/czebyszew_plot.jpg)

[N05] Numerical Integration

Numerical calculation of closed Newton-Cotes quadratures with diffrent interpolation polynomials.

![Interpolation polynomial](./N05/figures/trapezoid_plot.jpg)

[N06] Polynomial Root finding

Program for finding all roots of polynomial with Laguerre method and smoothing.

[N07] Solving Non-Linear Systems

Program solving multivariable systems of equations using damped multidimensional Newton method.

# Author

[Jakub Kurek] Jagiellonian University, Faculty of Physics, Astronomy and Applied Computer Science
