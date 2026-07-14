# GRIFFIN

GRIFFIN (Gauge-invariant Resonance In Four-Fermion INteractions) is an
object-oriented C++ library for electroweak radiative corrections. It provides
a modular description of 2-to-2 fermion scattering, with particular attention
to a gauge-invariant treatment of gauge-boson resonances.

The released version 1.1 provides Standard Model predictions with full NNLO
and selected leading higher-order contributions on the Z resonance,
NLO + O(αα<sub>s</sub>) off the Z resonance, and NLO Bhabha scattering.

This development branch additionally contains a linear tree-level
dimension-six SMEFT implementation and a new prediction engine that separates
the SM order, model order, combination rule, and output components. Its
current process coverage and release limitations are documented explicitly
below.

## Download and installation

The released version 1.1 is available
[on GitHub](https://github.com/lisongc/GRIFFIN/releases/tag/v1.1).

Build the static library out of source:

```bash
cmake -S . -B build
cmake --build build
```

Build the public examples, including the order-combination example, with:

```bash
cmake -S . -B build -DGRIFFIN_BUILD_EXAMPLES=ON
cmake --build build
./build/examples/testprediction
```

Run the automated tests with:

```bash
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

The convenience script performs a clean configured build, runs the complete
test suite, and executes the prediction-engine example:

```bash
./scripts/run_prediction_architecture_checks.sh
```

## Development documentation

- [SMEFT implementation](SMEFT_IMPLEMENTATION.md)
- [Prediction engine, class structure, and call flow](PREDICTION_ENGINE.md)
- [Numerical comparison with Jackson's implementation](SMEFT_JACKSON_COMPARISON.md)
- [Known issues and next development steps](KNOWN_ISSUES.md)

## Manual

The version 1.1 manual is available in the
[GRIFFIN manual repository](https://github.com/lisongc/GRIFFIN_manual).

## Using GRIFFIN in research

Please cite the [GRIFFIN paper](https://arxiv.org/pdf/2211.16272.pdf) when
using the library in research.

## License

The authors of this public repository consent to external users using,
reproducing, forking, and distributing the content with appropriate citation.

## Contacts

- Ayres Freitas: afreitas@pitt.edu
- Lisong Chen: lisong.chen@kit.edu
