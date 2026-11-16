#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if [[ -n "${SAGE_PYTHON:-}" ]]; then
    sage_python=("$SAGE_PYTHON")
else
    sage_python=(sage -python)
fi
mkdir -p results ../build
"${CXX:-g++}" -O3 -std=c++17 enumerate_curves.cpp -o ../build/ppas_enumerate_curves
for p in 7 11 13; do ../build/ppas_enumerate_curves "$p" "results/p$p"; done
"${sage_python[@]}" run_experiments.py verify
"${sage_python[@]}" run_experiments.py products
"${sage_python[@]}" restrictions.py
for p in 7 11 13; do "${sage_python[@]}" run_experiments.py richelot "$p"; done
"${sage_python[@]}" check_degrees.py
"${sage_python[@]}" local431.py
"${sage_python[@]}" plot_graphs.py
bash run_maximal.sh
