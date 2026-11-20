#!/usr/bin/env bash
set -euo pipefail
# Portable build: RUSTFLAGS from the environment overrides any target-cpu=native in ~/.cargo/config.toml
export RUSTFLAGS="-C target-cpu=x86-64-v3"
cargo install --locked --force --root "$CONDA_PREFIX" \
    --git https://github.com/rickbeeloo/barbell.git --rev 671cab1a barbell