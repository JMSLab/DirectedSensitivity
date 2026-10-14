# Directed sensitivity workflow

## What this does

Given a filled [`TEMPLATE.md`](TEMPLATE.md), report how far the object of interest $c$ can move, in each direction, across the nesting model $\mathcal{N}$, among primitives that fit the data at least as well as the working model.
Fit is the criterion $\hat{q}$; the budget is the working model's fit, $\bar{q} = \min_{\alpha\in \mathcal{A}}\hat{q}(\alpha) = \hat{q}(\hat{\alpha})$.
The feasible set is $\lbrace \nu\in\mathcal{N} : \hat{q}(\nu) \le \bar{q} \rbrace$, and a *certificate* is a feasible primitive realizing — or approaching — the extreme value of $c$ in one direction.

## How to run

Before starting, read the worked examples in [`examples/`](examples) for context.
Run all steps end to end without stopping.
Where a choice is required, such as setting a tolerance, make it and record it; do not pause for clearance.
Collect every such judgment call, and anything you could not verify, for the final [report](#outputs).

Verification of the nesting model may be approximate if necessary.
Shape and smoothness restrictions can be checked on a grid or up to a stated tolerance, and that is acceptable.
State the proxy you used; where a condition cannot be checked at all, prefer a construction that satisfies it by design over a check that cannot see it.

Every reported certificate must be a primitive the problem accepts as feasible (`IsFeasible`), reported with its criterion slack and its value of $c$.
Stay within the evaluation limits from [Computational constraints](TEMPLATE.md#5-computational-constraints-optional) and report the spend from `trace.Snapshot()`.

## Steps

### 0. Set up the problem in code

Read the filled template and build, with `SensitivityProblem` and `EvaluationTrace(limits=...)`:

- the object $c$ (`target`) and criterion $\hat{q}$ (`criterion`), from the code in [Location of data, code, and estimates](TEMPLATE.md#6-location-of-data-code-and-estimates) if provided, or from scratch if not;
- the nesting-model restrictions as named constraint callables, each nonnegative on $\mathcal{N}$;
- the baseline estimate $\hat{\alpha}$ (`estimate`) from [Location of data, code, and estimates](TEMPLATE.md#6-location-of-data-code-and-estimates), giving the budget $\bar{q} = \hat{q}(\hat{\alpha})$.

Confirm the baseline is admissible and that the working model nests (`NestsBaseline`).
If at any step you find a point in $\mathcal{A}$ with a lower criterion than $\hat{\alpha}$, make that point the baseline, recompute the budget, and record the change for the report.
Record any interpretation of a prose restriction, and any approximation, for the report.


### 1. Analytic exploration

Attempt to provide an **exact solution** (the supremum and infimum of $c$ over the feasible set) to the best of your ability using an analytic argument.
This is the most valuable deliverable from this step, and is worth substantial effort.

If you are unable to provide an exact solution, attempt a **problem reduction,** (e.g. from the infinite-dimensional $\mathcal{N}$ to a finite-dimensional problem) and/or proof of **informative bounds**, by analytic argument.

If [Workflow configuration](TEMPLATE.md#8-workflow-configuration) requests it, formalize the non-trivial mathematical claims in Lean, above all any claim that a value is the infimum or supremum.


### 2. Basis model

Propose and implement a basis model $\mathcal{B}$: a low-dimensional family in the nesting model that nests the working model.
Define it according to what [Analytic exploration](#1-analytic-exploration) produced.

- **Exact solution.** With $\hat{\nu}_L, \hat{\nu}_U$ the exact extremal primitives, deform between them and the working model,

  $$\mathcal{B} = \left\lbrace \beta = w \hat{\nu} + (1-w) \alpha : w \in [0,1], \alpha \in \mathcal{A}, \hat{\nu} \in \lbrace \hat{\nu}_L, \hat{\nu}_U \rbrace, \beta \in \mathcal{N} \right\rbrace,$$


  a one-dimensional extension of the working model in $w$; its certificate recovers the exact answer (up to attainment).
- **Problem reduction.** Replace the fixed anchor with the reduced family $\hat{\nu}(\theta)$, where $\theta$ ranges over the finite-dimensional parameter space $\Theta$ from the reduction,

  $$\mathcal{B} = \left\lbrace \beta = w \hat{\nu}(\theta) + (1-w) \alpha : w \in [0,1], \theta \in \Theta, \alpha \in \mathcal{A}, \beta \in \mathcal{N} \right\rbrace.$$


- **Bounds.** With $\hat{\nu}_L, \hat{\nu}_U$ the certificates achieving the bounds, take the union of the exact-solution deformation with a further finite-dimensional family $\mathcal{B}'\subseteq \mathcal{N}$ that may improve upon them,

  $$\mathcal{B} = \left\lbrace \beta = w \hat{\nu} + (1-w) \alpha : w \in [0,1], \alpha \in \mathcal{A}, \hat{\nu} \in \lbrace \hat{\nu}_L, \hat{\nu}_U \rbrace, \beta \in \mathcal{N} \right\rbrace \cup \mathcal{B}'.$$


- **None of the above.** Use the mathematical problem structure learned in [Analytic exploration](#1-analytic-exploration) to define a family $\mathcal{B}$ such that $\mathcal{A} \subset \mathcal{B} \subset \mathcal{N}$ and it is feasible to find informative certificates in $\mathcal{B}$ while staying within the evaluation limits.

Find the certificates: feasible primitives realizing the extreme $c$ in each direction, within the evaluation limits.
If it would induce a meaningful speedup, utilize symbolic or automatic differentiation instead of numerical differentiation.
The provided solvers may be used but are not required.
Include analytic certificates derived in [Analytic exploration](#1-analytic-exploration), as well as the output of numerical solvers.
Verify each certificate with `IsFeasible` and report its slack and $c$.
If there is an analytic certificate, the solver should return a certificate in the basis model that realizes or approaches it.


## Outputs

Write all of the following to the run directory.
Use fields and headings exactly as they appear in the templates.
Run `python -m DirectedSensitivity.MarkdownStyle` on every markdown file you write and fix what it reports; [`AGENTS.md`](../AGENTS.md) gives the conventions.

**`report.md`** — fill in the shell at [`REPORT.md`](REPORT.md).

**`certificates.json`** — Application-specific fields may be added alongside these, which are required:

```json
{
  "baseline": {"target": 0.0, "criterion": 0.0, "primitive": {}},
  "lower": {"target": 0.0, "criterion": 0.0, "slack": 0.0, "deviation": 0.0,
            "feasible": true, "primitive": {}, "source": "the solver or the anchor"},
  "upper": {"target": 0.0, "criterion": 0.0, "slack": 0.0, "deviation": 0.0,
            "feasible": true, "primitive": {}, "source": "the solver or the anchor"},
  "analytic": {"lower": 0.0, "upper": 0.0},
  "spend": "trace.Snapshot()",
  "limits": {"criterion": 0, "target": 0}
}
```

An `analytic` value is `null` on a side where [Analytic exploration](#1-analytic-exploration) produced none.
A `primitive` holds whatever $c$ and $\hat{q}$ read of it, in the parameters of the basis model.

**`animation.gif`** and **`certificates.png`** — a looping animation of the deformation path and a static figure of its ends, both from `DeformationPath` and one `WriteAnimation` call with `figure` set.
The loop runs out to one certificate and back, and carries no write-up: it is the plot and nothing else.
Embed it in [Results](REPORT.md#1-results) so the report reads on github.com.
If it would not be visually distracting, pass the observed data as `markers`.
Pass `DeformationPath` a feasible path: a sequence of primitives running from the lower certificate of [Results](REPORT.md#1-results), through the baseline estimate, to the upper one.
Every point on it must pass `IsFeasible`, and the ends must be the certificates the table reports.
