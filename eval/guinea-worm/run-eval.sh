#!/usr/bin/env bash
# A/B Evaluation Runner for Guinea Worm Stress Test of Restructured LASER Skill
#
# Usage:
#   ./eval/guinea-worm/run-eval.sh <prompt-number>
#   ./eval/guinea-worm/run-eval.sh 1          # Run prompt 1, both conditions
#   ./eval/guinea-worm/run-eval.sh all        # Run all prompts
#
# This script runs each prompt through two Claude sessions:
#   A) WITH skills (all three pipeline skills loaded normally)
#   B) WITHOUT skills (--disable-slash-commands blocks all skills)
#
# Outputs are saved to eval/guinea-worm/outputs/{with,without}-skill/prompt-N.md

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
    echo "PROMPT $num"
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
This is a public health epidemiological modeling project by the Bill & Melinda
Gates Foundation's Global Health division. The work supports guinea worm
eradication efforts by building spatial transmission models to evaluate
intervention strategies across Chad's endemic districts. All modeling uses the
open-source LASER framework (https://laser.idmod.org/) from the Institute for
Disease Modeling. This is purely computational epidemiology for public health
planning.
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
    for i in 1 2 3 4 5; do
        run_prompt "$i"
    done
elif [ -n "${1:-}" ]; then
    run_prompt "$1"
else
    echo "Usage: $0 <prompt-number|all>"
    echo "  $0 1       Run prompt 1 (both conditions)"
    echo "  $0 all     Run all 5 prompts"
    exit 1
fi

echo "Done. Review outputs in eval/guinea-worm/outputs/"
echo "Score using the rubric in eval/guinea-worm/rubric.md"
