+++
title = "Handling Non-Uniform Boundaries for WENO-5"
description = "Architecting zero-allocation, ghost-node-free boundary conditions for non-uniform WENO-5 using compile-time dispatch."
tags = ["gsoc", "julia", "sciml", "pde", "weno", "math"]
hasmath = true
hascode = true
date = Date(2026, 7, 13)
rss_title = "Non-Uniform Boundaries for WENO-5"
rss_description = "A deep dive into the mathematical and architectural implementation of non-uniform WENO-5 boundary conditions in MethodOfLines.jl."
+++

@@blog-post

# Handling Non-Uniform Boundaries for WENO-5

@@article-header
*GSoC 2026 · MethodOfLines.jl · Technical DevLog*

July 2026 · Boundary architecture
@@

\toc

In our previous explorations of the WENO-5 scheme for `MethodOfLines.jl`, we tackled the challenge of deriving dynamic smoothness indicators for non-uniform grids. But there is a glaring physical reality we must confront: What happens when our grid hits a wall?

When discretizing advection terms in Partial Differential Equations (PDEs), interior nodes are relatively forgiving. A standard 5-point stencil ($x_1$ to $x_5$) places the target node perfectly in the center at $x_3$. However, as we approach the physical boundaries of our domain, this symmetry shatters. We simply run out of spatial points.

Many traditional solvers handle this by extrapolating "ghost nodes" outside the physical domain. In our implementation (`nonuniform_weno.jl`), we took a different route. We strictly avoid ghost nodes (setting `extent = 0` for non-uniform grids) and compute the exact derivative $\partial u/\partial x$ using strictly one-sided, interior points. 

In this post, I will walk you through the mathematics and the software architecture behind our node-centered Lagrange WENO-5 boundary reconstruction, and how we engineered a single, allocation-free core to handle both the interior and the boundaries.

## 1. The Boundary Problem and The `Val{Target}` Architecture

To calculate a 4th-order accurate first derivative on a non-uniform grid, we always need a 5-point stencil. Let's visualize our local subset of the grid:

```text
x1 ---- x2 ---- x3 ---- x4 ---- x5
u1      u2      u3      u4      u5
```

For a standard interior node, we want to find the derivative at the center. Our target is $x_3$. But what if we need the derivative exactly at the left boundary wall? The physical points available to us are still $x_1$ through $x_5$, but our target node is now $x_1$. 

This means the sub-stencils we use to reconstruct the polynomial don't change, but the **point of evaluation** and the **Simpson integration cell** must shift dramatically.

If we were to write this naively, we would end up with massive `if-else` blocks checking whether the current node is a boundary, a near-boundary, or an interior node. Inside the hot loop of a PDE solver, these runtime branches destroy performance.

To solve this, we implemented a unified core function (`_weno_f_nonuniform_core`) driven entirely by Julia's multiple dispatch system using `Val{Target}`. 

We categorized every possible evaluation point in our 5-point stencil:
*   `Val(1)`: Left wall (target $x_1$)
*   `Val(2)`: Second node (target $x_2$)
*   `Val(3)`: Interior center node (target $x_3$)
*   `Val(4)`: Inner-right node (target $x_4$)
*   `Val(5)`: Right wall (target $x_5$)

By leveraging an explicit `Val{T}` parameter, we instruct the Julia compiler to generate highly optimized, target-specific machine code for each boundary condition. The scheme remains completely unified—the geometry, the Fornberg weights, and the ideal weight formulas automatically align themselves at compile-time based on the specific target node.

Let's open the black box and see exactly how this geometric shift affects our mathematical engine.

## 2. Redefining Geometry: Simpson Cells and Fornberg Weights

In the classic Jiang-Shu WENO scheme, we measure how "smooth" each sub-stencil is by calculating a $\beta_k$ indicator. This involves integrating the squared derivatives of our polynomial over a specific cell. Most modern implementations use **Simpson's Quadrature** to solve this integral numerically.

But here we face a critical problem: Where should the boundaries of this integration cell be?

For a standard interior node (`Val(3)`), our cell is perfectly symmetric around the target: $[x_{i-1/2}, x_{i+1/2}]$. But what if we are at the left wall (`Val(1)`) and our target is $x_1$? Physically, there is no $x_{1-1/2}$; that point falls into the void outside our domain. If we try to integrate over that range, our scheme bleeds out of bounds.

This is where our `Val{Target}` architecture shines. Instead of creating ghost nodes and extending the domain outward, **we shrink and shift the Simpson cell inward, strictly inside the physical domain.**

In our code, we handle this entirely at compile-time with the `_weno_target_geometry` function:

```julia
@inline _weno_target_geometry(::Val{1}, x1, x2, x3, x4, x5) = (x1, x1, (x1 + x2) / 2)
@inline _weno_target_geometry(::Val{3}, x1, x2, x3, x4, x5) = (x3, (x2 + x3) / 2, (x3 + x4) / 2)
@inline _weno_target_geometry(::Val{5}, x1, x2, x3, x4, x5) = (x5, (x4 + x5) / 2, x5)
```

As you can see, for an interior node `Val(3)`, our integration cell spans `[(x2+x3)/2, (x3+x4)/2]` (a full cell). But when we hit the left wall `Val(1)`, the cell doesn't spill over. It strictly starts at $x_1$ and shrinks to `[x1, (x1+x2)/2]` (a half cell). This guarantees we always stay within physical bounds.

### Dynamic Derivatives with Fornberg Weights

We have secured the left boundary ($x_L$), the right boundary ($x_{ph}$), and the midpoint ($x_M$) of our integration cell. To run the Simpson quadrature, we need to compute the first and second derivatives at these exact points.

Because we are on a non-uniform grid, the $\Delta x$ distances are constantly changing. We cannot use fixed, pre-calculated coefficients for our derivatives. Instead, we must generate them dynamically on the fly. To do this, we bring in Bengt Fornberg's (1988) algorithm.

We pass our 3-point sub-stencil and our current evaluation point ($x_t$) into the `_fornberg3_weights` function. It gives us the exact 0th, 1st, and 2nd derivative operator weights aligned perfectly for that specific point.

To avoid any runtime allocations, we bake these Fornberg calculations directly into the CPU cache using the `@inline` macro and `StaticArrays` (`SVector{3}`).

Now we have our custom-shifted cell and the dynamic Fornberg weights to compute derivatives at its boundaries. In the next step, we will use these tools to solve one of the biggest headaches of non-uniform grids: **Ideal Weights**.

## 3. The Ideal Weights and the Negative Weight Crisis

Once we have our three 3-point sub-stencil derivatives, we need to combine them. In a perfectly smooth region, WENO should combine these three 3-point stencils to perfectly recreate the accuracy of a full 5-point Lagrange derivative ($P'_{5pt}$). 

To achieve this, we use **ideal weights** ($d_0, d_1, d_2$). Mathematically, they must sum to 1 ($\sum d_k = 1$), representing a convex partition.

Because our target node shifts depending on our `Val{Target}`, the closed-form formulas for these ideal weights change for every single boundary condition. For example, our interior node `Val(3)` evaluates to:

```julia
d0 = ((x3 - x4) * (x3 - x5)) / ((x1 - x4) * (x1 - x5))
d2 = ((x3 - x1) * (x3 - x2)) / ((x5 - x1) * (x5 - x2))
```
*(Note: $d_1$ is simply $1 - d_0 - d_2$)*

If we perform a sanity check and assume a perfectly uniform grid ($\Delta x_i = \text{constant}$), this formula beautifully collapses back into the classic magic fractions: $(1/6, 2/3, 1/6)$. 

However, boundary nodes like `Val(1)` require massive, asymmetric algebraic formulas. But the real problem isn't the size of the algebra; the real problem is the physics of non-uniform grids.

**The Crisis:** On a non-uniform grid, there is absolutely no mathematical guarantee that $d_k$ will remain positive. When the grid stretches too rapidly, these ideal weights can become negative.

If an ideal weight becomes negative ($d_k < 0$), the foundational logic of WENO shatters. The combination is no longer convex, which leads to severe numerical instability and wild oscillations at the boundaries—exactly what WENO is designed to prevent.

### The Solution: Shi, Hu & Shu Weight Splitting

We cannot just throw away negative weights, because that would destroy our 4th-order accuracy. Instead, we implement a brilliant mathematical trick introduced by Shi, Hu & Shu (2002): **Positive/Negative Weight Splitting**.

Inside `_weno_f_nonuniform_core`, we define a splitting parameter $\theta = 3$. We then split every single ideal weight into a strictly positive component ($d_k^+$) and a negative counterpart ($d_k^-$):

$$d_k^+ = \frac{1}{2} (d_k + \theta |d_k|)$$
$$d_k^- = d_k^+ - d_k$$

Now that we have successfully isolated the negative weights, we compute *two entirely separate sets* of nonlinear weights. We calculate a positive set ($\omega_p$) using the $d_k^+$ weights, and a negative set ($\omega_m$) using the $d_k^-$ weights. Both sets are regularized by our smoothness indicators ($\beta_k$) and $\epsilon$:

```julia
ap0 = (dp0 / σp) / (ε + β0)^2
am0 = (dm0 / σm) / (ε + β0)^2
# ... normalized into ωp0 and ωm0
```

Finally, we construct two separate reconstructions: one for the positive weights ($R_p$) and one for the negative weights ($R_m$). In the final return statement, the negative oscillations are perfectly subtracted out:

```julia
return σp * Rp - σm * Rm
```
This guarantees that our scheme remains completely stable, highly accurate, and mathematically valid, even when the non-uniform grid tries to force negative weights upon us.

## 4. Engineering for Performance: Zero Allocations

Paper mathematics is beautiful, but mapping this heavy algebra into CPU instructions without melting the solver is an entirely different engineering challenge. 

Every component of this implementation—from the Fornberg weights to the Shi-Hu-Shu splitting—happens inside the hot loop of the PDE solver. If we allocate memory here, the garbage collector will paralyze the simulation.

To prevent this, the entire `_weno_f_nonuniform_core` is architected for zero allocations. 
*   **Compile-Time Dispatch:** By passing `Val{1}` through `Val{5}`, the Julia compiler hardcodes the specific geometry and ideal weight formulas for that exact boundary. The compiler eliminates the code for the boundaries we aren't using.
*   **Static Structures:** We strictly use `NTuple` and `SVector` from `StaticArrays.jl`. 
*   **Branchless Logic:** We replace standard `max` and conditional checks with `IfElse.ifelse` to keep the code branchless and GPU-friendly.
*   **Type Stability & AD Compatibility:** Due to the nature of the SciML ecosystem, it is highly critical that our code works seamlessly with Automatic Differentiation (AD). To ensure users can safely run their optimization loops or machine learning models, we took care to design the core to be fully type-stable. By ensuring all types are resolved and propagated correctly, our boundary scheme natively supports types like `Float32`, `Float64`, `ForwardDiff.Dual`, and `Symbolics.Num` without any issues.

The result is a unified, 4th-order accurate, non-uniform WENO-5 boundary reconstruction scheme that allocates exactly 0 bytes during the solver loop. No ghost nodes, no mathematical gaps, and no performance compromises.

## 5. Zooming Out: Multi-Dimensional Integration

Up to this point, we have focused entirely on a 1D stencil. But `MethodOfLines.jl` is designed to solve complex, multi-dimensional PDEs. How does our 1D boundary core handle 2D or 3D domains?

The beauty of the `MethodOfLines.jl` architecture is that our WENO core doesn't need to know about multiple dimensions. The solver tackles the domain direction by direction. 

During the discretization phase, an internal routing mechanism (`get_f_and_taps`) inspects the current node's index. If it detects that the node is approaching a boundary along the $x$-axis, it intercepts the standard interior call. Instead, it selects the exact boundary functor we defined (for example, `WENONonUniformBoundary{2}`) and passes a 1D slice (a `view`) of the $u$ and $x$ arrays directly into our core.

Because our core ignores uniform $dx$ assumptions and relies strictly on the raw $x$ coordinates, it blindly and efficiently computes the one-sided derivative for that specific slice. This separation of concerns means we didn't have to write messy "corner-case" logic for 2D or 3D boundaries. The 1D core handles the raw math, and the library handles the spatial routing.

## 6. Verification: How We Know It Works

In scientific computing, it is not enough for heavy algebra to work on paper; we must verify it through tests. To ensure our pipeline (Fornberg $\rightarrow$ Simpson $\beta$ $\rightarrow$ Shifted Ideal Weights $\rightarrow$ Shi-Hu-Shu Splitting) behaves as expected, we check the following fundamental criteria in our test suite:

1.  **Convexity Check:** We test that our dynamically generated ideal weights always sum to 1 ($\sum d_k = 1$) across all targets (`Val{1}` through `Val{5}`), even on variably spaced grids.
2.  **Polynomial Exactness:** We verify that the scheme can flawlessly reconstruct polynomials up to degree 2. This ensures that our positive/negative weight splitting operation does not destroy the linear backbone of the scheme.
3.  **4th-Order Convergence:** As a natural result of our formulation and the verified polynomial exactness, we observe in our base integrations that the scheme achieves its theoretically expected 4th-order convergence, even on non-uniform grids.

## Wrapping Up: Try It Yourself

By shrinking the Simpson cell, dynamically calculating Fornberg weights, and leveraging compile-time dispatch (`Val{Target}`), we successfully built a zero-allocation, ghost-node-free boundary scheme for non-uniform WENO-5. 

If you want to see this core engine in action without spinning up a full PDE solver, you can call it directly in the Julia REPL. Here is a minimal example showing how the same grid generates different derivatives simply by shifting the target:

```julia
using MethodOfLines

# Setup a non-uniform grid and some function values
const ε = (1.0e-6,)
xs = [0.0, 0.3, 0.9, 1.7, 2.2]
u  = sin.(xs)
dx = diff(xs)

# 1. Standard Interior Derivative (Defaults to Val{3}, target is x=0.9)
D_interior = MethodOfLines.weno_f_nonuniform(u, ε, 0.0, xs, dx)

# 2. Left Wall Derivative (Explicit Val{1}, target is x=0.0)
D_left = MethodOfLines.weno_f_nonuniform(u, ε, 0.0, xs, dx, Val(1))

# 3. Right Wall Derivative (Explicit Val{5}, target is x=2.2)
D_right = MethodOfLines.weno_f_nonuniform(u, ε, 0.0, xs, dx, Val(5))
```

This dynamic approach to boundaries ensures that `MethodOfLines.jl` can handle the variable-grid physics demanded by modern scientific simulations, all while keeping the hot loops incredibly fast. 

*This work is part of the Google Summer of Code 2026 project for SciML. You can explore the full implementation in [`nonuniform_weno.jl`](https://github.com/SciML/MethodOfLines.jl/blob/master/src/discretization/schemes/WENO/nonuniform_weno.jl).*

@@series-nav

**GSoC 2026 series**

- [Project Proposal](/gsoc-proposal/)
- [Kickoff](/kickoff/)
- [WENO Smoothness Indicators](/weno-smoothness/)
- **Non-Uniform Boundaries** *(this post)*
- [Showcase](/weno-showcase/)

@@

@@
