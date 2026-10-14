# Problem specification

Fill out this Markdown file.
After completing the verification steps below, hand the filled Markdown file to an LLM instructed to follow the workflow in [`WORKFLOW.md`](../../workflow/WORKFLOW.md).
Provide precise mathematical definitions in each section unless instructed otherwise.

The _object of interest_ $c(\nu) \in \mathbb{R}$ depends on an (unknown) parameter $\nu \in \mathcal{N}$, where the space $\mathcal{N}$ is the _nesting model_ and may be infinite-dimensional.
The data-dependent _criterion_ $\hat{q}(\nu) \in \mathbb{R}$ describes the fit of the nesting model to the data.
The _working model_ is a finite-dimensional subset $\mathcal{A}$ of $\mathcal{N}$.
The working model's fit is $\min_{\alpha\in \mathcal{A}}\hat q(\alpha)$, and the baseline estimate $\hat{\alpha} \in \mathcal{A}$ is the minimizer.

The objectives of this template are:
1. Precisely define the following objects:
   - The nesting model $\mathcal{N}$,
   - The working model $\mathcal{A}$,
   - The object of interest $c(\cdot)$, and
   - The criterion $\hat{q}(\cdot)$.
2. Specify the location of data, code, and estimates.

## 1. The nesting model

Define the nesting model $\mathcal{N}$, which encodes economic restrictions on unknown parameters, including unknown functions.
Define the domain and range of all unknown functions formally.
Specify any restrictions on the function class, such as continuity and shape restrictions.
The nesting model may also include a finite set of other known or unknown parameters.
Define the domain or values of these parameters.

**Example**. The nesting model $\mathcal{N}$ is the set of all demand functions $D: [P_0, P_1] \rightarrow [0, 1]$ that are:
- Continuous, and
- Strictly decreasing.

The bounds $P_0$ and $P_1$ are given in [Location of data, code, and estimates](#6-location-of-data-code-and-estimates).

| Your version: |
| --- |

Throughout, the primitive $\nu$ is the demand curve $D(p)$.
The cost parameters $(\gamma, \delta)$ are fixed.
They give the average cost curve $AC(p) = \gamma + \delta p$ and the marginal cost curve

$$MC(p) = AC(p) + \delta \frac{D(p)}{D'(p)}.$$

$\mathcal{N}$ is the set of demand curves $D : [0, 800] \to [0, 1]$ that satisfy the following properties:

- Continuous and differentiable.
- Strictly decreasing: $D'(p) < 0$ for all $p \in [0, 800]$.
- Bounded: $0 \le D(p) \le 1$ for all $p \in [0, 800]$.
- Single crossing of marginal cost and demand:
if $MC(p)>p$ for some $p\in [0, 800]$, then $MC(p')>p'$ for all $p'\in [0, p)$.
- There exists $\bar p \in [0, 800]$ with $D(\bar p)>0$ and $MC(p)<p$ for all $p\in (\bar p, 800]$.
- Nonnegative marginal cost: $MC(p) \ge 0$ for all $p \in [0, 800]$.

The cost parameters $(\gamma, \delta)$ are given in [Location of data, code, and estimates](#6-location-of-data-code-and-estimates).

<br>

## 2. The working model

Define the working model $\mathcal{A} \subset \mathcal{N}$.
Give its parameters, their domain, and the map from parameters to primitives.

**Example**. The working model is the isoelastic family

$$\mathcal{A} = \left\lbrace D_{a,b}(p)=\exp(a-b\log(p)): b > 0, a \le b \log P_0 \right\rbrace.$$

Parameters $a$ and $b$ map to the primitive $\alpha = D_{a,b} \in \mathcal{A}$.

| Your version: |
| --- |

$\mathcal{A}$ is the linear-demand family

$$\mathcal{A} = \left\lbrace D_{a,b}(p) = a + b p : D_{a,b} \in \mathcal{N} \right\rbrace,$$


with parameters $(a, b)$ mapping to the primitive $D_{a,b}$.

<br>

## 3. The object of interest

Define the object of interest $c(\cdot): \mathcal{N} \to \mathbb{R}$.
It must be well-defined and finite for every $\nu \in \mathcal{N}$.

**Example**. Writing $\nu = D$, the object is the loss in consumer surplus from a price increase $P_0 \to P_1$:

$$c(D) = \int_{P_0}^{P_1} D(p) \mathrm{d}p.$$

Since $D \in [0,1]$ on the bounded interval $[P_0, P_1]$, the integral lies in $[0, P_1 - P_0]$ and is finite for every $D \in \mathcal{N}$.

| Your version: |
| --- |

The competitive price $p_C$ is defined such that $p_C = AC(p_C)$, and the efficient price is the smallest price at which marginal cost is met from below:

$$p_E = \min\Big(\lbrace p\in[0,800] : p = MC(p) \rbrace\Big).$$

$c(D)$ is the deadweight-loss area between demand and marginal cost, from the competitive price $p_C$ to the efficient price $p_E$,

$$c(D) = \int_{D(p_C)}^{D(p_E)} \big[D^{-1}(Q) - MC\big(D^{-1}(Q)\big)\big]dQ.$$


Inverse demand $D^{-1}(Q)$ is well-defined for all $D \in \mathcal{N}$ because $D$ is strictly decreasing.
Inverse demand and marginal cost are bounded, so the integral exists and is finite for all $D \in \mathcal{N}$.

<br>

## 4. The criterion

Define the criterion $\hat{q}(\nu)$.
The criterion is a data-dependent measure of fit, where smaller values indicate a better fit.
The criterion must be finite and well-defined for every $\nu \in \mathcal{N}$.

**Example**. With observed shares $s_i$ at prices $p_i \in [P_0, P_1]$, for $i = 1, \dots, n$,

$$\hat{q}(D) = \frac{1}{n} \sum_{i=1}^{n} \big( s_i - D(p_i) \big)^2.$$

This evaluates $D\in[0,1]$ at finitely many points, so it is finite for every $D \in \mathcal{N}$.
The shares $s_i$ and prices $p_i$ are given under [Location of data, code, and estimates](#6-location-of-data-code-and-estimates).


| Your version: |
| --- |

$\hat{q}(D)$ is the employee-weighted sum of squared share residuals at the six observed prices,

$$\hat{q}(D) = \sum_{i=1}^{6} n_i \big( s_i - D(p_i) \big)^2,$$


where $(p_i, s_i, n_i)$ are the observed price, share, and employee count.

<br>

## 5. Computational constraints (optional)

State the largest number of evaluations of $c(\cdot)$ and $\hat{q}(\cdot)$ you are willing to spend, for when either is expensive.
Fill in the blanks below.

At most 100,000 evaluations of $c(\cdot)$ and at most 100,000 evaluations of $\hat{q}(\cdot)$.


## 6. Location of data, code, and estimates

Provide the location of data inputs.
Provide the location of code that implements the object of interest and criterion (or numerical approximations of these functions).
Provide the values of any parameters referenced in the prior sections.
Provide the location or values of the parameters estimated under the working model.

**Example**.
- $\lbrace (p_i, s_i) \rbrace_{i=1}^{n}$ — `demand.csv`, columns `price` and `share`.
- $c(D)$ — implemented in `consumer_surplus.py`.
- $\hat{q}(D)$ — implemented in `estimator.py`.
- $P_0, P_1$ — equal 100, 200.
- $\hat{a}, \hat{b}$ parameter estimates - `estimates.csv`, columns `intercept` and `slope`.

| Your version: |
| --- |

- $\lbrace (p_i, s_i, n_i) \rbrace_{i=1}^{6}$ and the average incremental costs — `EFC2010_table2.csv`, columns `relative_price`, `fraction_chose_contract_H`, `number_of_employees`, and `average_incremental_cost_contract_H`.
  Keyed in from Table II of Einav, Finkelstein, and Cullen (2010), https://doi.org/10.1162/qjec.2010.125.3.877.
- $(\gamma, \delta)$ — `specification.py:{AVG_COST_INTERCEPT,AVG_COST_SLOPE}`, the least-squares fit of average incremental cost on price with weights $n_i s_i$.
- $c(D)$ and $\hat{q}(D)$ — `specification.py:{Target,Criterion}`. `specification.py:LinearDemand` implements $D_{a,b}$.
- $\hat\alpha = D_{\hat a, \hat b}$ — `output/estimate.json`, produced by `estimate.py`, and read back by `specification.py:Baseline`.
  The least-squares fit of $s_i$ on $p_i$ with weights $n_i$ minimizes $\hat{q}$ over $\mathcal{A}$.
- Replication of Table III of the article — `output/replicate.txt`, produced by `replicate.py`, which also reads `EFC2010_table3.xlsx`, keyed in from Table III of the article.

<br>

## 7. Consistency and completeness checklist

After you complete the template, ask your LLM to verify the following prior to launching the workflow:

- [x] [The nesting model](#1-the-nesting-model), [The working model](#2-the-working-model), [The object of interest](#3-the-object-of-interest), and [The criterion](#4-the-criterion) are each filled with precise mathematical definitions.
- [x] The nesting model $\mathcal{N}$ is non-empty.
- [x] The working model is a genuine subset of the nesting model, $\mathcal{A} \subset \mathcal{N}$, and is finite-dimensional.
- [x] The object of interest $c(\nu)$ is well-defined and finite for every $\nu \in \mathcal{N}$.
- [x] The criterion $\hat{q}(\nu)$ is well-defined and finite for every $\nu \in \mathcal{N}$.
- [x] Every data or code symbol used in [The nesting model](#1-the-nesting-model) through [Computational constraints](#5-computational-constraints-optional) is resolved in [Location of data, code, and estimates](#6-location-of-data-code-and-estimates).

## 8. Workflow configuration

Should the LLM attempt a Lean formalization of mathematical claims?
Formalization may add time and token cost, and is not required to verify numerical certificates.

- [x] Yes
- [ ] No
