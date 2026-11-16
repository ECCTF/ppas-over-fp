#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if [[ -n "${SAGE_PYTHON:-}" ]]; then
    sage_python=("$SAGE_PYTHON")
else
    sage_python=(sage -python)
fi
mkdir -p results ../build
"${CXX:-g++}" -O3 -std=c++17 enumerate_maximal.cpp -o ../build/ppas_maximal
for p in 7 11 13 31; do
    ../build/ppas_maximal "$p" "results/p$p"
done
"${sage_python[@]}" maximal_graphs.py
"${sage_python[@]}" maximal_three.py
"${sage_python[@]}" verify_count_formulas.py
"${sage_python[@]}" verify_class_group_counts.py
"${sage_python[@]}" verify_automorphisms_duality.py
"${sage_python[@]}" verify_restriction_descent.py
"${sage_python[@]}" audit_current_results.py
"${sage_python[@]}" verify_connectivity.py
"${sage_python[@]}" verify_dyadic_seven.py
"${sage_python[@]}" verify_local_graph_formulas.py
"${sage_python[@]}" verify_rational_loops.py
"${sage_python[@]}" plot_maximal.py
"${sage_python[@]}" plot_complete_p13.py
