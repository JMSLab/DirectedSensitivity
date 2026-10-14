# DirectedSensitivity

`DirectedSensitivity` is an LLM workflow for functional form sensitivity analysis.

The workflow helps a user to evaluate the sensitivity of some object of interest (e.g., consumer surplus) to the form of some unknown function (e.g., a demand curve), when the object of interest has been estimated on some data under an assumed form (e.g., linear).

The workflow returns a lower bound on how much it is possible to increase or decrease the object of interest, without diminishing fit to the data, by selecting a different form for the function from within a (possibly large) class of functions.

Inputs to the workflow are given by filling out a template in markdown format; the LLM works on its own and returns its output in another markdown with an accompanying graphic.

Details on the method are given in our working paper:
Isaiah Andrews, Simon Essig Aberg, and Jesse M. Shapiro, "An LLM Workflow for Functional Form Sensitivity Analysis."
If you find this workflow useful in your work, please [cite that paper and this package](CITATION.cff).

## Applications

- The welfare loss from adverse selection in employer-provided health insurance a la [Einav, Finkelstein, and Cullen (2010)](https://doi.org/10.1162/qjec.2010.125.3.877); see [filled template](applications/EFC2010/filled_template.md) and [report](applications/EFC2010/output/report.md).
- The welfare effect of a fund-raising visit by a charity a la 
  [DellaVigna, List, and Malmendier (2012)](https://doi.org/10.1093/qje/qjr050); see [filled template](applications/DLM2012/filled_template.md) and [report](applications/DLM2012/output/report.md).
- The average effect of lobbying dollars on the probability of enacting an energy policy a la [Kang (2016)](https://doi.org/10.1093/restud/rdv029); see [filled template](applications/K2016/filled_template.md) and [report](applications/K2016/output/report.md).
- The discount factor
   - estimated from consumption data a la [Stock and Wright (2000)](https://doi.org/10.1111/1468-0262.00151); see [filled template](applications/SW2000/filled_template.md) and [report](applications/SW2000/output/report.md).
   - estimated from tax filing data a la [Martinez, Meier, and Sprenger (2023)](https://doi.org/10.1093/jeea/jvac067); see [filled template](applications/MMS2023/filled_template.md) and [report](applications/MMS2023/output/report.md).
 
All reports (and other application outputs) are LLM-generated.

## To use the workflow in your own application

1. Clone the repository and, from its root, install the package with `pip install .`.
2. Copy [`workflow/TEMPLATE.md`](workflow/TEMPLATE.md), fill it out for your problem, and have an LLM run its consistency and completeness checklist.
3. Ask a coding agent to follow [`workflow/WORKFLOW.md`](workflow/WORKFLOW.md) on the filled template.
   We recommend an agent with capabilities at least at the level of Claude Opus 5.x or GPT Sol 5.x.
4. Read the `report.md` written by the LLM.

(Optional): For automated proof verification, install [Lean 4](https://lean-lang.org/install/); see [`lean/`](lean) for the toolchain and Mathlib.


## To rerun the workflow on one of our applications:
1. Complete the above steps.
2. Install the applications with `pip install ".[applications]"`.
3. Download the original article and the files listed in the application's `README.md` to `temp/<APP>/`.
4. Navigate to the relevant subfolder of [`applications/`](applications), delete the files under `./output`, and ask a coding agent to follow [`workflow/WORKFLOW.md`](workflow/WORKFLOW.md) on the filled template.

## Repository

- [`workflow/`](workflow): the template, the workflow, the report shell, and the stylized examples.
- [`DirectedSensitivity/`](DirectedSensitivity): solvers and evaluation tracking.
- [`lean/`](lean): the Lean framework for the optional formalization.
- [`applications/`](applications): the applications; files under `/output` are LLM-generated.
- [`tests/`](tests): tests of the package.
- [`AGENTS.md`](AGENTS.md): conventions for coding agents working here or following the workflow.
