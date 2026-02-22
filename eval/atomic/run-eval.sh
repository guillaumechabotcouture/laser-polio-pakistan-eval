#!/usr/bin/env bash
# Suite 1: Atomic Competency Map — A/B Evaluation Runner
#
# Usage:
#   ./eval/atomic/run-eval.sh <prompt-id>
#   ./eval/atomic/run-eval.sh A01       # Run prompt A01, both conditions
#   ./eval/atomic/run-eval.sh A05       # Run prompt A05, both conditions
#   ./eval/atomic/run-eval.sh all       # Run all 20 prompts
#
# This script runs each prompt through two Claude sessions:
#   A) WITH skills (all skills loaded normally from project directory)
#   B) WITHOUT skills (--disable-slash-commands in clean temp directory)
#
# Outputs are saved to eval/atomic/outputs/{with,without}-skill/prompt-A{01-20}.md

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Ensure CLAUDECODE is unset so we can launch nested sessions
unset CLAUDECODE 2>/dev/null || true

# All 20 prompt IDs
ALL_PROMPTS=(A01 A02 A03 A04 A05 A06 A07 A08 A09 A10 A11 A12 A13 A14 A15 A16 A17 A18 A19 A20)

run_prompt() {
    local id=$1
    local prompt_file="$SCRIPT_DIR/prompt-${id}.txt"

    if [ ! -f "$prompt_file" ]; then
        echo "ERROR: $prompt_file not found."
        return 1
    fi

    echo "=========================================="
    echo "PROMPT $id"
    echo "=========================================="
    cat "$prompt_file"
    echo ""

    # --- Session A: WITH skills ---
    echo "--- Running WITH skills ---"
    local out_a="$SCRIPT_DIR/outputs/with-skill/prompt-${id}.md"
    mkdir -p "$(dirname "$out_a")"
    cd "$PROJECT_DIR"
    claude --print --dangerously-skip-permissions < "$prompt_file" > "$out_a" 2>/dev/null || true
    echo "  Saved to: $out_a ($(wc -c < "$out_a") bytes)"

    # --- Session B: WITHOUT skills (clean temp directory) ---
    echo "--- Running WITHOUT skills (clean temp dir) ---"
    local out_b="$SCRIPT_DIR/outputs/without-skill/prompt-${id}.md"
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
        < "prompt-${id}.txt" > "$out_b" 2>/dev/null || true
    echo "  Saved to: $out_b ($(wc -c < "$out_b") bytes)"
    rm -rf "$tmpdir"

    echo ""
}

# Main
if [ "${1:-}" = "all" ]; then
    for id in "${ALL_PROMPTS[@]}"; do
        run_prompt "$id"
    done
elif [ -n "${1:-}" ]; then
    run_prompt "$1"
else
    echo "Usage: $0 <prompt-id|all>"
    echo "  $0 A01     Run prompt A01 (both conditions)"
    echo "  $0 A05     Run prompt A05 (both conditions)"
    echo "  $0 all     Run all 20 prompts"
    echo ""
    echo "Prompt IDs: ${ALL_PROMPTS[*]}"
    exit 1
fi

echo "Done. Review outputs in eval/atomic/outputs/"
echo "Score using the rubric in eval/atomic/rubric.md"
