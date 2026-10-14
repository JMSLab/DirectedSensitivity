# Directed sensitivity: {APPLICATION}

This is the shell for the `report.md` that [`WORKFLOW.md`](WORKFLOW.md) requires.
Keep every heading and every table column below, in this order.
Fill each `{...}` slot and delete the braces.
**DO NOT WRITE ANYTHING OUTSIDE OF THE PROVIDED BRACES AND BE CONCISE.**
Report only numbers that code you ran produced.
Where a step failed, was skipped, or did not converge, say so in the slot where its result belongs.
Round reported numbers to around three significant figures for readability.
Delete this preface.

{ONE PARAGRAPH: the object of interest, the criterion, the working model, and the nesting model, one clause each; the path of the filled template; the paths of the code and outputs of this run; the coding agent and model identifier that executed this run; and the date.}

## 1. Results

| Direction | Object $c$ | Criterion $\hat q$ | Verified feasible? | Origin | Exact inf/sup? |
|---|---|---|---|---|---|
| Baseline | {} | {} | {`IsFeasible()`} | N/A | N/A |
| Lower | {} | {} | {`IsFeasible()`} | {Analytic exploration or Basis model} | {Yes/no} |
| Upper | {} | {} | {`IsFeasible()`} | {Analytic exploration or Basis model} | {Yes/no} |

![{ALT TEXT: what the loop shows, in one clause}](animation.gif)

## 2. Analytic exploration

{WHAT [Analytic exploration](WORKFLOW.md#1-analytic-exploration) ESTABLISHED, AS NUMBERED CLAIMS: an exact solution, a problem reduction, informative bounds, or nothing.
Proofs go in the [Appendix: proofs](#appendix-proofs), referenced by number.
If a Lean formalization was requested, report what compiled, what was left as an axiom or `sorry`, and the path of the Lean source.}

## 3. Basis model

{THE DEFINITION OF $\mathcal{B}$, MATHEMATICALLY, IN THE CONCISE FORM THE TEMPLATE USES FOR THE WORKING MODEL.
BE CONCISE, PRECISE, AND MATHEMATICAL.}

## 4. Numerical decisions

{NUMERICAL CHOICES THAT SHAPE THE REPORTED NUMBERS, ONE SHORT (1-2 sentences) BULLET EACH.
ONLY THE ESSENTIAL DECISIONS.
FEW OR NONE ARE EXPECTED.
ONLY INCLUDE NUMERICAL CHOICES THAT BEAR ON THE ANALYTIC OR NUMERICAL FINDINGS IN THE REPORT, OMITTING THOSE THAT BEAR ON THE ANIMATION.}

## 5. Judgement calls

{EACH CHOICE THE FILLED TEMPLATE DID NOT DETERMINE, ONE SHORT (1-2 sentences) BULLET EACH.
ONLY ESSENTIAL DECISIONS.
FEW OR NONE ARE EXPECTED.
ONLY INCLUDE JUDGMENT CALLS THAT BEAR ON THE ANALYTIC OR NUMERICAL FINDINGS IN THE REPORT, OMITTING THOSE THAT BEAR ON THE ANIMATION.}

## 6. Evaluation spend

{THE SPEND FROM `trace.Snapshot()`, AGAINST THE LIMITS.
NO DISCUSSION OR INTERPRETATION.}

## Appendix: proofs

{THE PROOF OF EACH CLAIM IN [Analytic exploration](#2-analytic-exploration), IN ORDER, EACH LABELLED AS IT WAS THERE.
BE CONCISE.}
