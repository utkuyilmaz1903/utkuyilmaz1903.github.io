+++
title = "GSoC 2026 Final Work Product: Non-Uniform Grid Support for MethodOfLines.jl"
description = "Official GSoC 2026 final work product: native non-uniform grid support for WENO-5 in MethodOfLines.jl."
tags = ["gsoc", "julia", "sciml", "pde", "weno"]
hasmath = true
hascode = true
date = Date(2026, 8, 17)
rss_title = "GSoC 2026 Final Work Product"
rss_description = "Official GSoC 2026 final work product: native non-uniform grid support for WENO-5 in MethodOfLines.jl."
+++

@@blog-post
@@work-product

# GSoC 2026 Final Work Product: Non-Uniform Grid Support for MethodOfLines.jl

@@article-header
*GSoC 2026 · MethodOfLines.jl · Final Work Product*

August 2026 · Official submission
@@

@@work-meta
**Contributor:** Utku Yılmaz ([@utkuyilmaz1903](https://github.com/utkuyilmaz1903))

**Organization:** [NumFOCUS](https://numfocus.org/) / [SciML](https://sciml.ai/)

**Project:** [MethodOfLines.jl](https://github.com/SciML/MethodOfLines.jl)

**Mentors:** Chris Rackauckas
@@

@@lead
This is my final work product submission for Google Summer of Code 2026.
@@

\toc

## Project Goals

MethodOfLines.jl discretizes symbolic PDE systems into ODE problems using finite differences. When I started, the high-order WENO advection scheme only worked on uniform grids, and the non-uniform code paths in general had little test coverage. My goal for the summer was to bring native non-uniform grid support to the WENO-5 scheme — from the core stencil math to boundaries, interfaces, tests, benchmarks, and documentation — without changing anything for existing uniform-grid users.

## What I Did

Roughly in order:

@@work-timeline

- **Started with tests.** Added a regression suite for the existing upwind scheme on non-uniform grids ([#562](https://github.com/SciML/MethodOfLines.jl/pull/562)) so I had a safe baseline before touching any scheme code.
- **Set up the dispatch architecture** ([#542](https://github.com/SciML/MethodOfLines.jl/pull/542)): uniform grids keep taking the existing path, vector grids route to the new one.
- **Wrote the non-uniform WENO-5 core** ([#581](https://github.com/SciML/MethodOfLines.jl/pull/581)): stencil weights and smoothness indicators computed from the actual grid coordinates instead of assuming constant spacing. The kernel is zero-allocation and checked against exact solutions.
- **Added boundary handling** ([#582](https://github.com/SciML/MethodOfLines.jl/pull/582)) with one-sided stencils near the domain edges, and **wired everything into the discretizer** ([#583](https://github.com/SciML/MethodOfLines.jl/pull/583)) so it works through the normal `discretize` workflow.
- **Added interface and periodic boundary support**, both for WENO ([#602](https://github.com/SciML/MethodOfLines.jl/pull/602)) and for the upwind scheme ([#572](https://github.com/SciML/MethodOfLines.jl/pull/572)). This includes joining two subdomains that use independently generated non-uniform grids.
- **Validated and measured it**: an MMS convergence test suite ([#591](https://github.com/SciML/MethodOfLines.jl/pull/591)) and a benchmark suite that now runs on every PR ([#609](https://github.com/SciML/MethodOfLines.jl/pull/609)).
- **Documented it**: a tutorial ([#615](https://github.com/SciML/MethodOfLines.jl/pull/615)) and a showcase page with two worked examples against exact solutions ([#616](https://github.com/SciML/MethodOfLines.jl/pull/616)), also published as a [blog post](https://utkuyilmaz1903.github.io/weno-showcase/).
- Along the way I helped keep CI green with some smaller fixes that weren't strictly part of the project but were blocking reviews ([#584](https://github.com/SciML/MethodOfLines.jl/pull/584), [#590](https://github.com/SciML/MethodOfLines.jl/pull/590), [#599](https://github.com/SciML/MethodOfLines.jl/pull/599), [#612](https://github.com/SciML/MethodOfLines.jl/pull/612), [#620](https://github.com/SciML/MethodOfLines.jl/pull/620), [#621](https://github.com/SciML/MethodOfLines.jl/pull/621)).

@@

A result I'm happy with, from the showcase: on a viscous shock problem, WENO on a grid clustered around the layer reached about two orders of magnitude lower error than a uniform grid with twice the points. Concentrating resolution where the physics needs it really pays off.

## Current State

Everything above is merged and released. A user can now pass any monotone grid vector and use WENO end to end:

```julia
disc = MOLFiniteDifference([x => xgrid], t; advection_scheme = WENOScheme())
prob = discretize(pdesys, disc)
```

Some limitations are documented in the tutorial: the formulation is non-conservative, mismatched-grid interfaces currently support first-order derivatives only, and the advection velocity needs to be continuous across an interface.

## What's Left to Do

The GSoC scope itself is complete. The natural next step is **grid adaptivity** — placing points well currently requires knowing where the steep features are in advance. I opened [issue #610](https://github.com/SciML/MethodOfLines.jl/issues/610) with a roadmap for this, and it's something I'd like to keep working on with the SciML community after GSoC. The limitations listed above are also good candidates for follow-up contributions.

## Merged PRs

All work was merged into [SciML/MethodOfLines.jl](https://github.com/SciML/MethodOfLines.jl):

@@pr-table

| PR | Description |
|---|---|
| [#562](https://github.com/SciML/MethodOfLines.jl/pull/562) | Non-uniform regression tests for the upwind scheme |
| [#542](https://github.com/SciML/MethodOfLines.jl/pull/542) | Dispatch architecture for non-uniform WENO |
| [#581](https://github.com/SciML/MethodOfLines.jl/pull/581) | Non-uniform WENO-5 core |
| [#582](https://github.com/SciML/MethodOfLines.jl/pull/582) | Boundary conditions for non-uniform WENO |
| [#583](https://github.com/SciML/MethodOfLines.jl/pull/583) | Discretizer integration |
| [#591](https://github.com/SciML/MethodOfLines.jl/pull/591) | MMS convergence tests |
| [#602](https://github.com/SciML/MethodOfLines.jl/pull/602) | Interface and periodic boundaries for non-uniform WENO |
| [#572](https://github.com/SciML/MethodOfLines.jl/pull/572) | Interface and periodic boundaries for non-uniform upwind |
| [#609](https://github.com/SciML/MethodOfLines.jl/pull/609) | Benchmark suite (runs on every PR) |
| [#615](https://github.com/SciML/MethodOfLines.jl/pull/615) | Tutorial |
| [#616](https://github.com/SciML/MethodOfLines.jl/pull/616) | Quantitative showcase |
| [#584](https://github.com/SciML/MethodOfLines.jl/pull/584), [#590](https://github.com/SciML/MethodOfLines.jl/pull/590), [#599](https://github.com/SciML/MethodOfLines.jl/pull/599), [#612](https://github.com/SciML/MethodOfLines.jl/pull/612), [#620](https://github.com/SciML/MethodOfLines.jl/pull/620), [#621](https://github.com/SciML/MethodOfLines.jl/pull/621) | Various CI and compatibility fixes |

@@

## Challenges and Learnings

@@learnings

- **The math came before the code.** The standard WENO formulas quietly assume uniform spacing, so I spent a good part of the summer deriving and verifying the non-uniform versions before writing any integration code. Testing against exact solutions from day one saved me many times.
- **Small, stacked PRs work.** Splitting the project into reviewable pieces made feedback faster and kept the uniform path verifiably untouched at every step.
- **Keeping CI green is part of the job.** When upstream changes broke tests mid-project, fixing them wasn't a distraction from the project — it was what made reviewing the project possible. That gave me a bit of a maintainer's perspective I didn't have before.
- **Measure honestly.** One of the more useful benchmark results was a negative one: a stretched grid on its own doesn't help unless it's actually adapted to the solution. That observation is what shaped the adaptivity roadmap.

@@

@@ack
Many thanks to my mentors Chris Rackauckas for their guidance throughout the summer, and to the SciML community for the reviews and discussions.
@@

@@series-nav

**GSoC 2026 series**

- [Project Proposal](/gsoc-proposal/)
- [Kickoff](/kickoff/)
- [WENO Smoothness Indicators](/weno-smoothness/)
- [Non-Uniform Boundaries](/nonuniform-weno-boundaries/)
- [Showcase](/weno-showcase/)
- **Final Work Product** *(this post)*

@@

@@
@@
