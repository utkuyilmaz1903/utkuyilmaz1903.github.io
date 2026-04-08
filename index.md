+++
title = "The Illusion of Uniformity"
tags = ["scientific-computing", "julia", "sciml", "pde"]
hasmath = true
hascode = true
+++

# The Illusion of Uniformity: The Need for Non-Uniform Grids and Dispatch Architecture in High-Order PDE Schemes

Finite difference discretization is a cornerstone of numerical PDE solving. For many standard physical models, uniform grids yield excellent and highly performant results. However, in scenarios where physical phenomena occur in very narrow regions with steep gradients—such as reaction interfaces in battery electrodes or sharp shockwaves—computing the entire domain with the same density leads to a massive waste of memory and CPU resources. 

In such models, to increase solution accuracy and prevent computational waste, we need non-uniform grids where node points are clustered specifically towards the boundaries where the physical event occurs.

![](/assets/grid_comparison.png)

## 1. Clustered Grids and the Stability Penalty

Clustering a grid into a specific region is not just a geometric operation; it directly impacts the stability limits of the differential equation. When we cluster points, the local grid spacing ($\Delta x$) suddenly shrinks. 

Considering the standard CFL (Courant-Friedrichs-Lewy) stability condition for diffusion equations:
\[ \Delta t \le \frac{(\Delta x)^2}{2\alpha} \]

This equation reveals a harsh reality: as the $\Delta x$ value decreases, the time step ($\Delta t$) must drop quadratically to microscopic levels for the system to remain stable. The fact that explicit solvers get choked at these microscopic time steps is the greatest proof of why we directly need the powerful implicit solvers and flexible discretization infrastructures present in the SciML ecosystem.

## 2. The Current Architectural Bottleneck in MethodOfLines.jl

When using `MethodOfLines.jl` (MOL) to solve these challenging problems, we can observe that basic advection schemes provide a foundational level of support for non-uniform grids via Fornberg algorithms. However, when we attempt to invoke high-order, shock-capturing schemes like `WENOScheme` on non-uniform grids to dampen the numerical oscillations generated at high-gradient boundaries, the system produces the following error:

`ERROR: LoadError: Scheme WENO applied to n(t, u)... is not defined.`

The root cause of this situation is that the existing mathematical operators for high-order schemes are defined assuming a scalar $\Delta x$ spacing. When the system is fed an `AbstractVector` (a grid array with variable steps), the current `DiscreteSpace` generator cannot handle this state, and the operators remain fundamentally undefined.

![](/assets/weno_necessity_real.png)

## 3. Comprehensive Vision and WENO Mathematics

Our goal is to overcome this architectural bottleneck and bring comprehensive, end-to-end non-uniform grid support to the MOL ecosystem. This vision includes refactoring basic schemes to handle variable step sizes flawlessly, mathematically deriving local stencils for WENO to build them as an isolated engine, and laying the groundwork for cell-volume data structures required for future Finite Volume Method (FVM) integration.

Specifically, in the most challenging part of the system—the WENO integration—we will implement the following mathematical pipeline:

**Step A: Dynamic Smoothness Indicators ($\beta$)**
Classic WENO uses constant fractions like $\frac{13}{12}$ optimized for uniform spacing. On non-uniform grids, these constants are invalid. As a solution, drawing from established literature (Shu, 1998), we will dynamically derive the smoothness indicators on-the-fly using Lagrange interpolation over local geometric distances ($h$):
\[ \beta = \int_{x_L}^{x_R} (P'(x))^2 dx + \int_{x_L}^{x_R} (P''(x))^2 dx \]

**Step B: Negative Weight Regularization (Shi-Hu-Shu Splitting)**
On highly stretched grids, the Fornberg algorithm can produce negative linear weights ($d_k$) within the stencil. This causes WENO to risk division by zero and collapses the Partition of Unity rule.
To prevent this, we apply the Shi-Hu-Shu (2002) approach, splitting the weights into two strictly positive sub-groups:
\[ \tilde{d}_k^+ = \frac{1}{2} (d_k + 3|d_k|) \quad \text{and} \quad \tilde{d}_k^- = \tilde{d}_k^+ - d_k \]
Thus, both groups remain safely positive ($\sigma^+ = \sum \tilde{d}_k^+$ and $\sigma^- = \sum \tilde{d}_k^-$).

**Step C: Final Non-Linear Weights ($\omega$)**
We run the standard WENO formula for these two split groups and normalize them:
\[ \alpha_k^\pm = \frac{\tilde{d}_k^\pm}{(\epsilon + \beta_k)^2} \quad \implies \quad \omega_k^\pm = \frac{\alpha_k^\pm}{\sum \alpha_j^\pm} \]
To unify the system back into a single final derivative coefficient, we stitch the non-linear weights together:
\[ \omega_k^{final} = \sigma^+ \cdot \omega_k^+ - \sigma^- \cdot \omega_k^- \]
This exact mathematical flow zeroes out oscillations at shockwaves while guaranteeing the stability of the system.

## 4. Architectural Routing via Multiple Dispatch

To avoid disrupting the library's existing high-performance structure, the integration will be entirely built upon Julia's Multiple Dispatch capability. 

When the system detects a `StepRangeLen` (legacy uniform grids), it will route directly to the existing high-performance operators, ensuring zero breaking changes. However, when it detects an `AbstractVector` (non-uniform grid), it will bypass the current generic generator and directly trigger the dynamic weight calculators formulated above, which are fully isolated within the `WENO` module.

*(Conclusion)*
Building this infrastructure for MethodOfLines.jl will open the door to much more realistic physical simulations on irregular geometries. I will explore deeper topics, such as the formulaic derivation of dynamic smoothness indicators and how these 1D schemes will be embedded into the Multi-Dimensional architecture, in the technical blog posts I will publish in the upcoming period. 

Thanks for reading.