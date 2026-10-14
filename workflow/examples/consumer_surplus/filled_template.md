# Problem specification

Fill out this Markdown file.
After completing the verification steps below, hand the filled Markdown file to an LLM instructed to follow the workflow in [`WORKFLOW.md`](../../WORKFLOW.md).
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

The nesting model $\mathcal{N}$ is the set of all demand functions $D: [P_0, P_1] \rightarrow [0, 1]$ that are:
- Continuous, and
- Strictly decreasing.

The bounds $P_0$ and $P_1$ are given in [Location of data, code, and estimates](#6-location-of-data-code-and-estimates).

<br>

## 2. The working model

Define the working model $\mathcal{A} \subset \mathcal{N}$.
Give its parameters, their domain, and the map from parameters to primitives.

**Example**. The working model is the isoelastic family

$$\mathcal{A} = \left\lbrace D_{a,b}(p)=\exp(a-b\log(p)): b > 0, a \le b \log P_0 \right\rbrace.$$

Parameters $a$ and $b$ map to the primitive $\alpha = D_{a,b} \in \mathcal{A}$.

| Your version: |
| --- |

The working model is the isoelastic family

$$\mathcal{A} = \left\lbrace D_{a,b}(p)=\exp(a-b\log(p)): b > 0, a \le b \log P_0 \right\rbrace.$$

Parameters $a$ and $b$ map to the primitive $\alpha = D_{a,b} \in \mathcal{A}$.

<br>

## 3. The object of interest

Define the object of interest $c(\cdot): \mathcal{N} \to \mathbb{R}$.
It must be well-defined and finite for every $\nu \in \mathcal{N}$.

**Example**. Writing $\nu = D$, the object is the loss in consumer surplus from a price increase $P_0 \to P_1$:

$$c(D) = \int_{P_0}^{P_1} D(p) \mathrm{d}p.$$

Since $D \in [0,1]$ on the bounded interval $[P_0, P_1]$, the integral lies in $[0, P_1 - P_0]$ and is finite for every $D \in \mathcal{N}$.

| Your version: |
| --- |

Writing $\nu = D$, the object is the loss in consumer surplus from a price increase $P_0 \to P_1$:

$$c(D) = \int_{P_0}^{P_1} D(p) \mathrm{d}p.$$

Since $D \in [0,1]$ on the bounded interval $[P_0, P_1]$, the integral lies in $[0, P_1 - P_0]$ and is finite for every $D \in \mathcal{N}$.

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

With observed shares $s_i$ at prices $p_i \in [P_0, P_1]$, for $i = 1, \dots, n$,

$$\hat{q}(D) = \frac{1}{n} \sum_{i=1}^{n} \big( s_i - D(p_i) \big)^2.$$

This evaluates $D\in[0,1]$ at finitely many points, so it is finite for every $D \in \mathcal{N}$.
The shares $s_i$ and prices $p_i$ are given under [Location of data, code, and estimates](#6-location-of-data-code-and-estimates).

<br>

## 5. Computational constraints (optional)

State the largest number of evaluations of $c(\cdot)$ and $\hat{q}(\cdot)$ you are willing to spend, for when either is expensive.
Fill in the blanks below.

At most ___ evaluations of $c(\cdot)$ and at most ___ evaluations of $\hat{q}(\cdot)$.


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

- $\lbrace (p_i, s_i) \rbrace_{i=1}^{n}$ — `demand.csv`, columns `price` and `share`.
- $c(D)$ — implemented in `consumer_surplus.py`.
- $\hat{q}(D)$ — implemented in `estimator.py`.
- $P_0, P_1$ — equal 100, 200.
- $\hat{a}, \hat{b}$ parameter estimates - `estimates.csv`, columns `intercept` and `slope`.

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

- [ ] Yes
- [x] No
