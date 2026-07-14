#!/bin/sh

set -eu

ROOT=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
BUILD_DIR=${1:-"$ROOT/build-prediction-checks"}
JOBS=${JOBS:-4}

cmake -S "$ROOT" -B "$BUILD_DIR" \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON \
  -DGRIFFIN_BUILD_EXAMPLES=ON
cmake --build "$BUILD_DIR" -j"$JOBS"
ctest --test-dir "$BUILD_DIR" --output-on-failure

"$BUILD_DIR/examples/testprediction"
