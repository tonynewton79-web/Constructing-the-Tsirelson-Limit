# Constructing the Tsirelson Limit

## Finite Weyl Observables and Double-Cone Modular Geometry for Bosonic Bell–CHSH Violations

**Tony Newton**  
Newton Astro Labs*  
London, UK  
tony.newton79@gmail.com  

**S. P. Sorella**  
Instituto de Física Armando Dias Tavares  
Universidade do Estado do Rio de Janeiro  
Rua São Francisco Xavier 524  
20550-900 Maracanã, Rio de Janeiro, Brazil  
silvio.sorella@fis.uerj.br  

> *Newton Astro Labs is the trading name under which Tony Newton conducts independent computational research in the United Kingdom; it is not a limited company.*

---

## Overview

This repository contains the **reproducibility and certification suite** for the paper

> **Constructing the Tsirelson Limit: Finite Weyl Observables and Double-Cone Modular Geometry for Bosonic Bell–CHSH Violations**

The project studies finite-band bosonic Bell–CHSH constructions built from Weyl observables and modular/double-cone geometry. The computational objective is not merely to find large numerical Bell values, but to attach **proof-carrying certificates** to the optimization: exact rational reductions, polynomial boundedness checks, branch and shell certification, LP/Farkas certificates, high-precision independent regressions, and endpoint analysis.

The repository is intentionally compact. The complete public replay suite is consolidated into **three Python files** so that the principal claims can be checked without a large external framework.

---

## Main certified results

### 1. Degree-15 global optimum at \(\sigma=0.02\)

For the declared degree-15 finite-Weyl optimization problem at

\[
\sigma=0.02,
\]

the end-to-end certified Bell–CHSH optimum satisfies

\[
\boxed{
2.36461694035383
<
\mathcal B_{\max}^{(15)}(0.02)
<
2.36461694036
}
\]

with release status

```text
DEGREE15_GLOBAL_OPTIMUM_END_TO_END_CERTIFIED
```

The certification chain includes exact rational and symbolic components rather than relying only on floating-point optimization.

---

### 2. Smaller bandwidth improves the certified construction

The \(\sigma=0.01\) replay gives a rigorous degree-15 witness that **beats the certified \(\sigma=0.02\) value**.

This establishes that the observed improvement toward the zero-bandwidth regime is not merely a floating-point optimizer artifact.

---

### 3. Zero-bandwidth endpoint

For the degree-15 zero-bandwidth problem,

\[
\sigma=0,
\]

the current globally certified bracket is

\[
\boxed{
2.56173695411020
<
\mathcal B_{\max}^{(15)}(0)
<
2.56173836501593
}
\]

and the isolated KKT stationary candidate is

\[
\boxed{
\mathcal B_{\mathrm{KKT}}
=
2.561736954110666788\ldots
}
\]

The candidate has been uniquely isolated within the certified stationary analysis, but the repository **does not claim that the exact global optimum has been proved equal to that candidate**.

Accordingly, the stronger status

```text
ZERO_BANDWIDTH_GLOBAL_OPTIMUM_CERTIFIED
```

is **not** asserted at this stage.

The publication-safe conclusion is:

> the zero-bandwidth degree-15 optimum is globally enclosed in the certified interval above, with a unique isolated KKT stationary candidate extremely close to the certified lower edge; equality between the global optimum and that candidate remains a separate proof obligation.

---

## Relation to the Tsirelson bound

The quantum Bell–CHSH Tsirelson bound is

\[
2\sqrt 2
\approx
2.828427124746190\ldots
\]

The present repository does **not** claim saturation of \(2\sqrt2\), nor does it claim a universal optimization over every possible quantum observable.

Instead, it gives a proof-carrying construction and optimization program inside a specified **finite Weyl / degree-15 observable grammar**, together with a controlled bandwidth continuation toward the \(\sigma=0\) endpoint.

The distinction is important:

```text
finite-grammar global certification
            ≠
unrestricted Tsirelson saturation
```

The research target is to understand how close a rigorously controlled finite construction can move toward the Tsirelson regime and which structural mechanisms control that approach.

---

# Repository contents

The current repository consists of three consolidated replay files.

| File | Purpose |
|---|---|
| `test_01_sigma002_end_to_end.py` | Exact end-to-end certification of the \(\sigma=0.02\), degree-15 result |
| `test_02_zero_bandwidth_global.py` | Global zero-bandwidth certification and endpoint analysis |
| `test_03_bandwidth_reflection_regressions.py` | Independent high-precision regressions, reflection branch tests, \(\sigma=0.01\) witness, and deterministic bandwidth continuation |

No auxiliary JSON certificate is required for the first replay: its exact exterior ledger is compressed and embedded directly in the Python source.

---

# Test 1 — \(\sigma=0.02\) end-to-end certificate

Run

```bash
python test_01_sigma002_end_to_end.py
```

This file consolidates the publication tests for:

- exact rational Gaussian-moment reconstruction;
- canonical covariance construction;
- Bell-kernel enclosure;
- exact rational Sturm boundedness of the degree-15 witness;
- exact local contact-defect shell replay;
- exact exterior LP/Farkas certification;
- replay of the embedded exact exterior ledger;
- integrity checking of the embedded ledger by SHA-256.

The source records the embedded-ledger hash

```text
983f184f4a1d8eca8ddb61534bc4a24f374c52a1c382967d7f97f1b0b6439796
```

and the exterior certificate contains

```text
401 leaves
368 PRUNE
33 EMPTY
```

The file is designed to run using the Python standard library.

### Expected scientific conclusion

A successful replay supports

```text
DEGREE15_GLOBAL_OPTIMUM_END_TO_END_CERTIFIED
```

for the declared \(\sigma=0.02\), degree-15 optimization problem, with

```text
2.36461694035383 < B_max^(15)(0.02) < 2.36461694036
```

---

# Test 2 — zero-bandwidth global certificate

Run

```bash
python test_02_zero_bandwidth_global.py
```

This is the large endpoint/global-certificate replay.

Its role is to certify the degree-15 \(\sigma=0\) problem through the global branch/orthant decomposition and associated exact certificate machinery, and to separate two logically different statements:

1. **global enclosure of the optimum**, and
2. **identification of the exact optimizer/value**.

The first has been certified to the bracket

```text
2.56173695411020 < B_max^(15)(0) < 2.56173836501593
```

The KKT analysis isolates the stationary candidate

```text
2.561736954110666788...
```

but the second statement — exact equality of the global optimum with that candidate — is deliberately not promoted without the remaining equality proof.

This distinction is part of the certificate protocol and should not be removed from downstream summaries of the result.

---

# Test 3 — bandwidth and reflection regressions

Run

```bash
python test_03_bandwidth_reflection_regressions.py
```

This file consolidates several independent numerical and exact regression layers:

- independent high-precision reconstruction of the \(\sigma=0.02\) covariance;
- \(L=8\) Fejér normalization checks;
- finite-Weyl Bell values;
- reflection correlation;
- Fejér-reflection benchmarks;
- exact rational Sturm boundedness for the optimized hybrid axis;
- numerical evaluation of the optimized hybrid construction;
- end-to-end exact \(\sigma=0.01\) witness certification;
- deterministic degree-15 bandwidth continuation toward

\[
K(0)=I_8.
\]

Representative regression values encoded in the test include

\[
\nu_A
\approx
45.00841733645926200048,
\]

\[
\nu_B
\approx
45.03061897035379099827,
\]

\[
\kappa_q
\approx
45.00840364042463202550,
\]

\[
\kappa_p
\approx
44.98621295944687856821,
\]

and the \(L=8\) Fejér normalization

\[
M_8
\approx
0.96036378670045260128.
\]

The regression suite also checks representative Bell quantities including

\[
B_{\rm raw}
\approx
2.0949957173756400870,
\]

\[
B_{\rm Fej}
\approx
2.2714939727791660463,
\]

and a hybrid benchmark near

\[
B_{\rm hyb}
\approx
2.55990385767281.
\]

These regression numbers are intermediate consistency checks. They should not be confused with the globally certified degree-15 optimum values stated above.

---

# Requirements

Recommended environment:

```text
Python >= 3.10
```

Test 1 is written against the Python standard library.

For the complete three-file suite, install:

```bash
python -m pip install mpmath numpy scipy
```

Optional test-runner support:

```bash
python -m pip install pytest
```

A minimal setup is therefore

```bash
python -m venv .venv
```

Linux/macOS:

```bash
source .venv/bin/activate
```

Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

then

```bash
python -m pip install --upgrade pip
python -m pip install mpmath numpy scipy pytest
```

---

# Running the complete replay

The canonical direct-execution sequence is

```bash
python test_01_sigma002_end_to_end.py
python test_02_zero_bandwidth_global.py
python test_03_bandwidth_reflection_regressions.py
```

Run the files from the repository root.

Where pytest discovery is supported by the individual file, it may also be used, for example:

```bash
pytest -q test_01_sigma002_end_to_end.py
pytest -q test_03_bandwidth_reflection_regressions.py
```

For publication/release verification, direct execution of all three files is preferred because it exercises the complete scripted replay path.

---

# What is being certified?

The repository separates **discovery** from **release certification**.

A large numerical Bell value by itself is not treated as sufficient evidence of a theorem-level result. The release chain uses several complementary mechanisms.

## Exact rational arithmetic

Where possible, coefficients, polynomial identities, interval endpoints, and certificate objects are represented exactly using rational arithmetic.

## Sturm certification

Polynomial boundedness and absence/presence of roots on declared intervals are checked using exact Sturm-sequence logic rather than sampled plotting.

## LP/Farkas certificates

Exterior regions are discharged by linear-programming/Farkas-type certificates. The purpose is to turn a numerical exclusion into a replayable proof object.

## Shell and contact-defect replay

Near candidate contact regions, the optimization is subdivided and checked with a dedicated local certificate rather than assuming that a coarse exterior argument remains valid.

## Independent high-precision regression

`test_03_bandwidth_reflection_regressions.py` reconstructs central quantities independently at high precision so that the exact machinery is not the only line of defense against transcription or implementation errors.

## Endpoint/KKT isolation

At zero bandwidth, stationary-point analysis is kept logically separate from the global enclosure. An isolated KKT point is evidence for the optimizer but is not automatically promoted to equality with the global optimum.

---

# Certificate logic

The intended logical structure is

```text
physical / modular construction
          |
          v
finite Weyl observable grammar
          |
          v
exact or high-precision moment reconstruction
          |
          v
degree-15 Bell functional
          |
          v
local + exterior domain decomposition
          |
          +--> exact Sturm checks
          |
          +--> LP / Farkas certificates
          |
          +--> shell/contact replay
          |
          v
global interval enclosure
          |
          +--> independent numerical regression
          |
          +--> bandwidth continuation
          |
          v
publication-level certified statement
```

At \(\sigma=0\), the final step deliberately branches:

```text
global interval enclosure  --------> CERTIFIED

unique isolated KKT candidate ------> CERTIFIED AS A CANDIDATE

global optimum = KKT candidate -----> NOT YET CLAIMED
```

---

# Reproducibility policy

A result should be regarded as reproduced only when the corresponding script terminates without assertion failure.

The scripts use assertions as certificate gates. A failed assertion should therefore be treated as a failed replay rather than silently bypassed.

For a clean independent reproduction:

1. clone the repository;
2. create a fresh Python environment;
3. install only the dependencies listed above;
4. run the three scripts without modifying coefficients or tolerances;
5. retain the terminal output;
6. record the Python and package versions;
7. compare any reported hashes with those embedded in the source.

For stronger archival reproducibility, record

```bash
python --version
python -m pip freeze
```

alongside the test output.

---

# Suggested repository layout

The current three test files are sufficient to expose the computational core. For a fuller archival release, the repository can use the following layout:

```text
Constructing-the-Tsirelson-Limit/
|
|-- README.md
|-- LICENSE
|-- CITATION.cff
|-- requirements.txt
|
|-- test_01_sigma002_end_to_end.py
|-- test_02_zero_bandwidth_global.py
|-- test_03_bandwidth_reflection_regressions.py
|
|-- paper/
|   |-- Constructing_the_Tsirelson_Limit.tex
|   |-- Constructing_the_Tsirelson_Limit.pdf
|   `-- references.bib
|
`-- results/
    |-- replay_environment.txt
    `-- replay_output.txt
```

Only the three Python files are required for the present computational replay. The other files are recommended for a complete archival/publication package.

---

# Recommended `requirements.txt`

A simple repository-wide dependency file is

```text
mpmath
numpy
scipy
pytest
```

If exact version pinning is desired for an archival release, generate it from the verified environment after the final replay.

---

# Recommended GitHub Actions check

A future continuous-integration workflow can run the compact tests on each push.

Because the zero-bandwidth certificate file is substantially larger and may be computationally heavier, it can either be included in the normal workflow or placed in a separate manual/release job.

The important scientific rule is that CI should **replay the certificates**, not replace them with looser surrogate tests.

---

# Scope of the claims

This repository supports claims only for the mathematical and computational model actually encoded in the scripts.

In particular, it does not by itself prove:

- universal saturation of the Tsirelson bound;
- optimality over all bounded quantum observables;
- optimality over every bosonic state or every modular construction;
- experimental realization of the proposed observables;
- device-independent experimental certification;
- equality of the zero-bandwidth global optimum with the isolated KKT candidate.

The strongest current result is a **global proof-carrying optimization inside the declared finite degree-15 Weyl grammar**, together with rigorous bandwidth-improvement and zero-bandwidth endpoint analysis.

---

# Why the zero-bandwidth distinction matters

It is tempting to identify

\[
2.561736954110666788\ldots
\]

with the exact degree-15 endpoint optimum because it lies extremely close to the lower edge of the certified global bracket.

That inference is not made here.

A unique isolated stationary candidate and a global enclosing interval are two different proof objects. The repository keeps them separate until the final equality step has its own certificate.

This is why the project currently reports

```text
ZERO-BANDWIDTH GLOBAL BRACKET: CERTIFIED
KKT CANDIDATE ISOLATION:       CERTIFIED
EXACT GLOBAL EQUALITY:         OPEN
```

rather than prematurely promoting the stronger terminal label.

---

# Connection between the three tests

The files are complementary rather than redundant.

```text
test_01
  |
  |  proves the principal finite-band global result at sigma=0.02
  v
certified finite-band baseline
  |
  |  test_03 checks independent reconstruction,
  |  reflection structure, sigma=0.01 improvement,
  |  and continuation toward zero bandwidth
  v
bandwidth trajectory
  |
  |  test_02 attacks the limiting sigma=0 problem globally
  v
certified zero-bandwidth bracket + isolated KKT candidate
```

Together they provide a compact evidence chain from a finite nonzero-bandwidth theorem to the controlled endpoint problem.

---

# Scientific interpretation

The computational evidence supports three increasingly strong conclusions.

**First**, a nontrivial Bell–CHSH violation can be constructed and globally optimized within the declared finite-Weyl degree-15 grammar.

**Second**, decreasing the bandwidth from \(\sigma=0.02\) to \(\sigma=0.01\) produces a rigorously improved witness, supporting a genuine bandwidth-dependent mechanism rather than an isolated optimizer accident.

**Third**, the zero-bandwidth endpoint admits a substantially larger globally certified Bell value, with an isolated stationary candidate near \(2.561736954110666788\ldots\).

What remains open is whether the current finite grammar can be analytically pushed further toward \(2\sqrt2\), whether increasing degree changes the endpoint structure, and whether the isolated zero-bandwidth KKT candidate can be promoted to the exact global optimizer.

---

# Open proof targets

The most immediate mathematical targets are:

1. close the remaining equality argument between the zero-bandwidth global optimum and the isolated KKT candidate;
2. extend the exact global machinery beyond degree 15;
3. determine the degree dependence of the certified optimum;
4. determine whether the finite-Weyl sequence converges and, if so, identify its limiting value;
5. test whether the limiting construction reaches \(2\sqrt2\) or a strictly smaller grammar-dependent ceiling;
6. derive analytic structure explaining the numerically observed bandwidth improvement;
7. preserve exact/replayable certification as the grammar and degree are enlarged.

These are research targets, not claims already established by this repository.

---

# Paper

**Title**

> *Constructing the Tsirelson Limit: Finite Weyl Observables and Double-Cone Modular Geometry for Bosonic Bell–CHSH Violations*

**Authors**

Tony Newton and S. P. Sorella

The manuscript is intended to be read together with this repository: the paper provides the mathematical and physical derivation, while the repository provides the compact computational replay/certificate layer.

When an arXiv identifier or DOI is assigned, it should be added here together with the final bibliographic citation.

---

# Citation

Until a DOI/arXiv identifier is available, cite the work by title and authors:

```text
Tony Newton and S. P. Sorella,
"Constructing the Tsirelson Limit: Finite Weyl Observables and
Double-Cone Modular Geometry for Bosonic Bell–CHSH Violations."
```

A `CITATION.cff` file can be added once the final public version and persistent identifier are frozen.

---

# Integrity and verification

The computational package is designed around independently replayable proof objects rather than trust in a single optimizer run.

Important numerical constants should therefore be interpreted together with their corresponding certificate type:

| Result | Status |
|---|---|
| \(\sigma=0.02\), degree-15 global optimum bracket | **Certified** |
| End-to-end \(\sigma=0.02\) certificate chain | **Certified** |
| \(\sigma=0.01\) improvement over \(\sigma=0.02\) | **Certified** |
| Zero-bandwidth global interval | **Certified** |
| Zero-bandwidth KKT candidate isolation | **Certified** |
| Equality of zero-bandwidth optimum and KKT candidate | **Not yet certified** |
| Universal Tsirelson saturation | **Not claimed** |

---

# License

No software license is implied by the existence of a public GitHub repository.

If the repository does not yet contain a `LICENSE` file, normal copyright restrictions apply. A specific open-source license should be added explicitly if reuse, redistribution, or modification rights are intended.

---

# Contact

### Tony Newton

Newton Astro Labs  
London, UK  
tony.newton79@gmail.com

### S. P. Sorella

Instituto de Física Armando Dias Tavares  
Universidade do Estado do Rio de Janeiro  
Rua São Francisco Xavier 524  
20550-900 Maracanã  
Rio de Janeiro, Brazil  
silvio.sorella@fis.uerj.br

---

# Computational Verification Statement

All computational algorithms, proof-search procedures, and certificate criteria were developed by the author(s). Large Language Models (LLMs) such as ChatGPT and Claude were utilized strictly as computational assistants for debugging code and cross-checking routines. No mathematical claims rely on LLM outputs; all results are grounded exclusively in independently reproducible, rigorously verified certificates (exact, symbolic, rational, interval, exhaustive, or controlled high-precision).
