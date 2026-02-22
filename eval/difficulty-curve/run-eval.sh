#!/usr/bin/env bash
# A/B Evaluation Runner for Difficulty Curve (Suite 3)
#
# Usage:
#   ./eval/difficulty-curve/run-eval.sh <level>
#   ./eval/difficulty-curve/run-eval.sh D1         # Run level D1, both conditions
#   ./eval/difficulty-curve/run-eval.sh D5         # Run level D5, both conditions
#   ./eval/difficulty-curve/run-eval.sh all        # Run all 8 levels
#
# This script runs each prompt through two Claude sessions:
#   A) WITH skills (all three pipeline skills loaded normally)
#   B) WITHOUT skills (--disable-slash-commands blocks all skills)
#
# Outputs are saved to eval/difficulty-curve/outputs/{with,without}-skill/prompt-DN.md

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Ensure CLAUDECODE is unset so we can launch nested sessions
unset CLAUDECODE 2>/dev/null || true

run_prompt() {
    local num=$1
    local prompt_file="$SCRIPT_DIR/prompt-${num}.txt"

    if [ ! -f "$prompt_file" ]; then
        echo "ERROR: $prompt_file not found."
        return 1
    fi

    echo "=========================================="
    echo "LEVEL $num"
    echo "=========================================="
    cat "$prompt_file"
    echo ""

    # --- Session A: WITH skills ---
    echo "--- Running WITH skills ---"
    local out_a="$SCRIPT_DIR/outputs/with-skill/prompt-${num}.md"
    mkdir -p "$(dirname "$out_a")"
    cd "$PROJECT_DIR"
    claude --print --dangerously-skip-permissions < "$prompt_file" > "$out_a" 2>/dev/null || true
    echo "  Saved to: $out_a ($(wc -c < "$out_a") bytes)"

    # --- Session B: WITHOUT skills (clean temp directory) ---
    echo "--- Running WITHOUT skills (clean temp dir) ---"
    local out_b="$SCRIPT_DIR/outputs/without-skill/prompt-${num}.md"
    mkdir -p "$(dirname "$out_b")"

    # Create isolated environment with no LASER-specific knowledge
    local tmpdir
    tmpdir=$(mktemp -d)
    cp "$prompt_file" "$tmpdir/"
    cat > "$tmpdir/CLAUDE.md" << 'CTXEOF'
# Project Context
This is a public health epidemiological modeling project using the open-source
LASER framework (https://laser.idmod.org/) from the Institute for Disease
Modeling. The project builds spatial disease transmission models for public
health planning.
CTXEOF

    cd "$tmpdir"
    claude --print --dangerously-skip-permissions --disable-slash-commands \
        < "prompt-${num}.txt" > "$out_b" 2>/dev/null || true
    echo "  Saved to: $out_b ($(wc -c < "$out_b") bytes)"
    rm -rf "$tmpdir"

    echo ""
}

# Main
if [ "${1:-}" = "all" ]; then
    for level in D1 D2 D3 D4 D5 D6 D7 D8; do
        run_prompt "$level"
    done
elif [ -n "${1:-}" ]; then
    run_prompt "$1"
else
    echo "Usage: $0 <level|all>"
    echo "  $0 D1      Run level D1 (both conditions)"
    echo "  $0 D5      Run level D5 (both conditions)"
    echo "  $0 all     Run all 8 levels (D1-D8)"
    exit 1
fi

echo "Done. Review outputs in eval/difficulty-curve/outputs/"
echo "Score using the rubric in eval/difficulty-curve/rubric.md"
