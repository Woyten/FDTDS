#!/usr/bin/env bash
set -eu

echo "=== Compiling with Rust ==="
cargo build --release

echo "=== Compiling with C ==="
g++ -O3 *.c

echo "=== Running with Rust ==="
time target/release/fdtds -n out-rust -o bin3d -g 200

echo "=== Running with C ==="
time ./a.out -n out-c -o bin3d -g 200
