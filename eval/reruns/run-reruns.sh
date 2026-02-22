#!/usr/bin/env bash
# Suite 2: Statistical Re-Runs — A/B Evaluation Runner
#
# Re-runs the polio and guinea worm A/B tests 3 times each to establish
# score distributions and confidence intervals.
#
# Usage:
#   ./eval/reruns/run-reruns.sh all              # All diseases, all 3 runs (60 sessions)
#   ./eval/reruns/run-reruns.sh polio             # Polio only, all 3 runs (30 sessions)
#   ./eval/reruns/run-reruns.sh guinea-worm       # Guinea worm only, all 3 runs (30 sessions)
#   ./eval/reruns/run-reruns.sh polio 2           # Polio only, run 2 only (10 sessions)
#   ./eval/reruns/run-reruns.sh all 1             # All diseases, run 1 only (20 sessions)
#
# Outputs are saved to:
#   eval/reruns/outputs/{polio,guinea-worm}/run-{1,2,3}/{with,without}-skill/prompt-{1-5}.md
#
# Total sessions: up to 60 (5 prompts x 2 conditions x 2 diseases x 3 runs)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
EVAL_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
OUTPUT_BASE="$SCRIPT_DIR/outputs"

# Ensure CLAUDECODE is unset so we can launch nested Claude sessions
unset CLAUDECODE 2>/dev/null || true

# ---------------------------------------------------------------------------
# Prompt and CLAUDE.md source paths
# ---------------------------------------------------------------------------
POLIO_PROMPT_DIR="$EVAL_DIR"
GUINEA_PROMPT_DIR="$EVAL_DIR/guinea-worm"

# ---------------------------------------------------------------------------
# CLAUDE.md content for WITHOUT condition (clean tmpdir)
# ---------------------------------------------------------------------------
POLIO_CLAUDE_MD='# Project Context
This is a public health epidemiological modeling project by the Bill & Melinda
Gates Foundation'"'"'s Global Health division. The work supports polio eradication
efforts by building spatial transmission models to evaluate vaccination
strategies across Pakistan'"'"'s districts. All modeling uses the open-source LASER
framework (https://laser.idmod.org/) from the Institute for Disease Modeling.
This is purely computational epidemiology for public health planning.'

GUINEA_CLAUDE_MD='# Project Context
This is a public health epidemiological modeling project by the Bill & Melinda
Gates Foundation'"'"'s Global Health division. The work supports guinea worm
eradication efforts by building spatial transmission models to evaluate
intervention strategies across Chad'"'"'s endemic districts. All modeling uses the
open-source LASER framework (https://laser.idmod.org/) from the Institute for
Disease Modeling. This is purely computational epidemiology for public health
planning.'

# ---------------------------------------------------------------------------
# Counters for progress reporting
# ---------------------------------------------------------------------------
TOTAL_SESSIONS=0
COMPLETED_SESSIONS=0

count_sessions() {
    local diseases=("$@")
    local runs=("${RUN_LIST[@]}")
    TOTAL_SESSIONS=$(( ${#diseases[@]} * ${#runs[@]} * 5 * 2 ))
}

# ---------------------------------------------------------------------------
# Core: run a single prompt for a single disease/run/condition
# ---------------------------------------------------------------------------
run_single_prompt() {
    local disease=$1
    local run_num=$2
    local prompt_num=$3
    local condition=$4  # "with-skill" or "without-skill"

    # Resolve prompt file
    local prompt_dir
    if [ "$disease" = "polio" ]; then
        prompt_dir="$POLIO_PROMPT_DIR"
    else
        prompt_dir="$GUINEA_PROMPT_DIR"
    fi
    local prompt_file="$prompt_dir/prompt-${prompt_num}.txt"

    if [ ! -f "$prompt_file" ]; then
        echo "  ERROR: $prompt_file not found. Skipping."
        return 1
    fi

    # Output path
    local out_file="$OUTPUT_BASE/$disease/run-${run_num}/${condition}/prompt-${prompt_num}.md"
    mkdir -p "$(dirname "$out_file")"

    # Skip if output already exists (resume support)
    if [ -f "$out_file" ] && [ -s "$out_file" ]; then
        echo "  SKIP (already exists): $out_file"
        return 0
    fi

    if [ "$condition" = "with-skill" ]; then
        # --- WITH skills: run in project directory (skills available) ---
        cd "$PROJECT_DIR"
        claude --print --dangerously-skip-permissions \
            < "$prompt_file" > "$out_file" 2>/dev/null || true
    else
        # --- WITHOUT skills: run in clean tmpdir (no skills) ---
        local tmpdir
        tmpdir=$(mktemp -d)

        # Copy prompt file
        cp "$prompt_file" "$tmpdir/prompt-${prompt_num}.txt"

        # Write disease-appropriate CLAUDE.md
        if [ "$disease" = "polio" ]; then
            printf '%s\n' "$POLIO_CLAUDE_MD" > "$tmpdir/CLAUDE.md"
        else
            printf '%s\n' "$GUINEA_CLAUDE_MD" > "$tmpdir/CLAUDE.md"
        fi

        cd "$tmpdir"
        claude --print --dangerously-skip-permissions --disable-slash-commands \
            < "prompt-${prompt_num}.txt" > "$out_file" 2>/dev/null || true
        rm -rf "$tmpdir"
    fi

    local bytes
    bytes=$(wc -c < "$out_file" | tr -d ' ')
    echo "  Saved: $out_file ($bytes bytes)"
}

# ---------------------------------------------------------------------------
# Run all prompts for a given disease and run number
# ---------------------------------------------------------------------------
run_disease_run() {
    local disease=$1
    local run_num=$2

    for prompt_num in 1 2 3 4 5; do
        for condition in with-skill without-skill; do
            COMPLETED_SESSIONS=$((COMPLETED_SESSIONS + 1))
            echo ""
            echo "=========================================="
            echo "  Run $run_num/3 | Disease: $disease | Prompt $prompt_num/5 | Condition: $condition"
            echo "  Session $COMPLETED_SESSIONS/$TOTAL_SESSIONS"
            echo "=========================================="

            run_single_prompt "$disease" "$run_num" "$prompt_num" "$condition"
        done
    done
}

# ---------------------------------------------------------------------------
# Timing utilities
# ---------------------------------------------------------------------------
format_duration() {
    local secs=$1
    local h=$((secs / 3600))
    local m=$(( (secs % 3600) / 60 ))
    local s=$((secs % 60))
    if [ $h -gt 0 ]; then
        printf '%dh %dm %ds' $h $m $s
    elif [ $m -gt 0 ]; then
        printf '%dm %ds' $m $s
    else
        printf '%ds' $s
    fi
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
DISEASE_ARG="${1:-}"
RUN_ARG="${2:-all}"

# Validate disease argument
case "$DISEASE_ARG" in
    polio)
        DISEASE_LIST=(polio)
        ;;
    guinea-worm)
        DISEASE_LIST=(guinea-worm)
        ;;
    all)
        DISEASE_LIST=(polio guinea-worm)
        ;;
    "")
        echo "Usage: $0 {polio|guinea-worm|all} [run-number|all]"
        echo ""
        echo "Arguments:"
        echo "  Disease:  polio | guinea-worm | all"
        echo "  Run:      1 | 2 | 3 | all (default: all)"
        echo ""
        echo "Examples:"
        echo "  $0 all              # All diseases, all 3 runs (60 sessions)"
        echo "  $0 polio            # Polio only, all 3 runs (30 sessions)"
        echo "  $0 guinea-worm 2    # Guinea worm only, run 2 (10 sessions)"
        echo "  $0 all 1            # All diseases, run 1 only (20 sessions)"
        echo ""
        echo "Outputs: eval/reruns/outputs/{disease}/run-{N}/{with,without}-skill/prompt-{1-5}.md"
        echo ""
        echo "Note: Existing outputs are skipped (resume support). Delete outputs to re-run."
        exit 1
        ;;
    *)
        echo "ERROR: Unknown disease '$DISEASE_ARG'. Use: polio | guinea-worm | all"
        exit 1
        ;;
esac

# Validate run argument
case "$RUN_ARG" in
    1|2|3)
        RUN_LIST=("$RUN_ARG")
        ;;
    all)
        RUN_LIST=(1 2 3)
        ;;
    *)
        echo "ERROR: Unknown run number '$RUN_ARG'. Use: 1 | 2 | 3 | all"
        exit 1
        ;;
esac

# Count total sessions for progress
count_sessions "${DISEASE_LIST[@]}"

echo "============================================================"
echo "  Suite 2: Statistical Re-Runs"
echo "============================================================"
echo "  Diseases:  ${DISEASE_LIST[*]}"
echo "  Runs:      ${RUN_LIST[*]}"
echo "  Prompts:   5 per disease"
echo "  Conditions: with-skill, without-skill"
echo "  Total sessions: $TOTAL_SESSIONS"
echo "  Output base: $OUTPUT_BASE"
echo "============================================================"
echo ""

START_TIME=$(date +%s)

for run_num in "${RUN_LIST[@]}"; do
    for disease in "${DISEASE_LIST[@]}"; do
        run_disease_run "$disease" "$run_num"
    done
done

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "============================================================"
echo "  COMPLETE"
echo "============================================================"
echo "  Sessions completed: $COMPLETED_SESSIONS/$TOTAL_SESSIONS"
echo "  Wall time: $(format_duration $ELAPSED)"
echo "  Output base: $OUTPUT_BASE"
echo ""
echo "  Next steps:"
echo "    1. Score all outputs using the rubrics:"
echo "       - Polio:       eval/rubric.md"
echo "       - Guinea worm: eval/guinea-worm/rubric.md"
echo "    2. Record scores in: eval/reruns/analysis.md"
echo "    3. Compute confidence intervals and variance analysis"
echo "============================================================"
