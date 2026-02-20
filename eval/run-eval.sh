#!/usr/bin/env bash
# A/B Evaluation Runner for LASER Skill
#
# Usage:
#   ./eval/run-eval.sh <prompt-number>
#   ./eval/run-eval.sh 1          # Run prompt 1, both conditions
#   ./eval/run-eval.sh all        # Run all prompts
#
# This script runs each prompt through two Claude sessions:
#   A) WITH the LASER skill (project .claude/skills/ loaded normally)
#   B) WITHOUT the skill (--deny flag blocks it)
#
# Outputs are saved to eval/outputs/{with,without}-skill/prompt-N.md

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
PROMPTS_FILE="$SCRIPT_DIR/prompts.md"

# Extract a specific prompt by number from prompts.md
extract_prompt() {
    local num=$1
    # Extract text between ```  blocks under "### Prompt N"
    awk "/^### Prompt ${num}:/,/^\`\`\`$/" "$PROMPTS_FILE" \
        | sed -n '/^```$/,/^```$/p' \
        | sed '1d;$d'
}

run_prompt() {
    local num=$1
    local prompt
    prompt=$(extract_prompt "$num")

    if [ -z "$prompt" ]; then
        echo "ERROR: Could not extract prompt $num"
        return 1
    fi

    echo "=========================================="
    echo "PROMPT $num"
    echo "=========================================="
    echo "$prompt"
    echo ""

    # --- Session A: WITH skill ---
    echo "--- Running WITH skill ---"
    local out_a="$SCRIPT_DIR/outputs/with-skill/prompt-${num}.md"
    cd "$PROJECT_DIR"
    echo "$prompt" | claude --print 2>/dev/null > "$out_a" || true
    echo "  Saved to: $out_a"

    # --- Session B: WITHOUT skill ---
    echo "--- Running WITHOUT skill ---"
    local out_b="$SCRIPT_DIR/outputs/without-skill/prompt-${num}.md"
    cd "$PROJECT_DIR"
    echo "$prompt" | claude --print --deny "Skill(laser-spatial-disease-modeling)" 2>/dev/null > "$out_b" || true
    echo "  Saved to: $out_b"

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

echo "Done. Review outputs in eval/outputs/"
echo "Score using the rubric in eval/rubric.md"
