#!/usr/bin/env bash
# Unified containerized eval runner for all suites
#
# Uses Docker to isolate WITH and WITHOUT conditions, preventing the WITHOUT
# condition from reading laser source code via --dangerously-skip-permissions.
#
# Usage:
#   ./eval/docker/run-all-containerized.sh <suite> [-j N]
#
#   suite: atomic | difficulty | reruns | multi-turn | all
#   -j N:  max parallel jobs (default 4)
#
# Prerequisites:
#   - Docker Desktop running
#   - ~/.claude/ directory with valid OAuth credentials
#
# Outputs go to outputs-containerized/ dirs (preserving original results):
#   eval/atomic/outputs-containerized/{with,without}-skill/
#   eval/difficulty-curve/outputs-containerized/{with,without}-skill/
#   eval/reruns/outputs-containerized/{polio,guinea-worm}/run-{1,2,3}/{with,without}-skill/
#   eval/multi-turn/outputs-containerized/{with,without}-skill/

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
EVAL_DIR="$PROJECT_DIR/eval"

# Docker image names
IMG_WITH="laser-eval-with"
IMG_WITHOUT="laser-eval-without"

# ---------------------------------------------------------------------------
# CLAUDE.md content for WITHOUT condition (per-disease variants)
# ---------------------------------------------------------------------------
CLAUDE_MD_GENERIC='# Project Context
This is a public health epidemiological modeling project using the open-source
LASER framework (https://laser.idmod.org/) from the Institute for Disease
Modeling. The project builds spatial disease transmission models for public
health planning.'

CLAUDE_MD_POLIO='# Project Context
This is a public health epidemiological modeling project by the Bill & Melinda
Gates Foundation'"'"'s Global Health division. The work supports polio eradication
efforts by building spatial transmission models to evaluate vaccination
strategies across Pakistan'"'"'s districts. All modeling uses the open-source LASER
framework (https://laser.idmod.org/) from the Institute for Disease Modeling.
This is purely computational epidemiology for public health planning.'

CLAUDE_MD_GUINEA='# Project Context
This is a public health epidemiological modeling project by the Bill & Melinda
Gates Foundation'"'"'s Global Health division. The work supports guinea worm
eradication efforts by building spatial transmission models to evaluate
intervention strategies across Chad'"'"'s endemic districts. All modeling uses the
open-source LASER framework (https://laser.idmod.org/) from the Institute for
Disease Modeling. This is purely computational epidemiology for public health
planning.'

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
SUITE="${1:-}"
MAX_JOBS=4

shift || true
while [ $# -gt 0 ]; do
    case "$1" in
        -j)
            MAX_JOBS="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

if [ -z "$SUITE" ]; then
    echo "Usage: $0 <suite> [-j N]"
    echo ""
    echo "  suite:  atomic | difficulty | reruns | multi-turn | all"
    echo "  -j N:   max parallel jobs (default 4)"
    echo ""
    echo "Examples:"
    echo "  $0 all                # Run all 4 suites"
    echo "  $0 atomic -j 8       # Run atomic suite with 8 parallel jobs"
    echo "  $0 reruns             # Run reruns suite"
    echo "  $0 multi-turn         # Run multi-turn suite"
    exit 1
fi

# ---------------------------------------------------------------------------
# Preflight checks
# ---------------------------------------------------------------------------
if ! docker info >/dev/null 2>&1; then
    echo "ERROR: Docker is not running. Start Docker Desktop and try again."
    exit 1
fi

# Auth: read OAuth token from file or environment.
# Generate with: claude setup-token > ~/.claude/oauth-token
TOKEN_FILE="$HOME/.claude/oauth-token"
if [ -n "${CLAUDE_CODE_OAUTH_TOKEN:-}" ]; then
    DOCKER_AUTH_TOKEN="$CLAUDE_CODE_OAUTH_TOKEN"
elif [ -f "$TOKEN_FILE" ]; then
    DOCKER_AUTH_TOKEN="$(cat "$TOKEN_FILE")"
else
    echo "ERROR: No auth token found."
    echo "  Either set CLAUDE_CODE_OAUTH_TOKEN or create $TOKEN_FILE"
    echo "  Generate one with: claude setup-token"
    exit 1
fi

# ---------------------------------------------------------------------------
# Job parallelism (portable — works on macOS bash 3.2)
# ---------------------------------------------------------------------------
declare -a JOB_PIDS=()
COMPLETED=0
FAILED=0

wait_for_slot() {
    while true; do
        if [ ${#JOB_PIDS[@]} -eq 0 ]; then
            return
        fi
        local new_pids=()
        for pid in "${JOB_PIDS[@]}"; do
            if kill -0 "$pid" 2>/dev/null; then
                new_pids+=("$pid")
            else
                wait "$pid" 2>/dev/null && COMPLETED=$((COMPLETED + 1)) || FAILED=$((FAILED + 1))
            fi
        done
        JOB_PIDS=("${new_pids[@]+"${new_pids[@]}"}")

        if [ ${#JOB_PIDS[@]} -lt "$MAX_JOBS" ]; then
            return
        fi
        sleep 0.5
    done
}

wait_all() {
    if [ ${#JOB_PIDS[@]} -eq 0 ]; then
        return
    fi
    for pid in "${JOB_PIDS[@]}"; do
        wait "$pid" 2>/dev/null && COMPLETED=$((COMPLETED + 1)) || FAILED=$((FAILED + 1))
    done
    JOB_PIDS=()
}

# ---------------------------------------------------------------------------
# Core: run one prompt through Docker
# ---------------------------------------------------------------------------
run_docker_claude() {
    local prompt_file=$1   # absolute path on host
    local output_file=$2   # absolute path on host
    local condition=$3     # "with-skill" or "without-skill"
    local claude_md=${4:-$CLAUDE_MD_GENERIC}

    mkdir -p "$(dirname "$output_file")"

    # Resume support: skip if output already exists and is non-empty
    if [ -f "$output_file" ] && [ -s "$output_file" ]; then
        echo "  SKIP (exists): $output_file"
        return 0
    fi

    if [ "$condition" = "with-skill" ]; then
        docker run --rm -i \
            -e CLAUDE_CODE_OAUTH_TOKEN="$DOCKER_AUTH_TOKEN" \
            -v "$PROJECT_DIR:/workspace:ro" \
            "$IMG_WITH" \
            sh -c 'cd /workspace && claude --print --dangerously-skip-permissions' \
            < "$prompt_file" > "$output_file" 2>/dev/null || true
    else
        # Create tmpdir with minimal CLAUDE.md (no skills, no laser source)
        local tmpdir
        tmpdir=$(mktemp -d)
        printf '%s\n' "$claude_md" > "$tmpdir/CLAUDE.md"

        docker run --rm -i \
            -e CLAUDE_CODE_OAUTH_TOKEN="$DOCKER_AUTH_TOKEN" \
            -v "$tmpdir:/eval-work:ro" \
            "$IMG_WITHOUT" \
            sh -c 'cd /eval-work && claude --print --dangerously-skip-permissions --disable-slash-commands' \
            < "$prompt_file" > "$output_file" 2>/dev/null || true

        rm -rf "$tmpdir"
    fi

    local bytes
    bytes=$(wc -c < "$output_file" 2>/dev/null | tr -d ' ')
    echo "  Done: $(basename "$output_file") ($bytes bytes) [$condition]"
}

# ---------------------------------------------------------------------------
# Build Docker images
# ---------------------------------------------------------------------------
build_images() {
    echo "============================================================"
    echo "  Building Docker images"
    echo "============================================================"
    local t0
    t0=$(date +%s)

    echo "  Building $IMG_WITHOUT (base, no laser)..."
    docker build --target without-laser -t "$IMG_WITHOUT" \
        -f "$SCRIPT_DIR/Dockerfile" "$SCRIPT_DIR" >/dev/null 2>&1

    echo "  Building $IMG_WITH (base + laser-generic)..."
    docker build --target with-laser -t "$IMG_WITH" \
        -f "$SCRIPT_DIR/Dockerfile" "$SCRIPT_DIR" >/dev/null 2>&1

    local elapsed=$(( $(date +%s) - t0 ))
    echo "  Images built in ${elapsed}s"
    echo ""
}

# ---------------------------------------------------------------------------
# Suite 1: Atomic Competency Map (20 prompts)
# ---------------------------------------------------------------------------
run_suite_atomic() {
    echo "============================================================"
    echo "  Suite 1: Atomic Competency Map (20 prompts x 2 conditions)"
    echo "============================================================"
    local prompt_dir="$EVAL_DIR/atomic"
    local out_base="$EVAL_DIR/atomic/outputs-containerized"

    for id in A01 A02 A03 A04 A05 A06 A07 A08 A09 A10 \
              A11 A12 A13 A14 A15 A16 A17 A18 A19 A20; do
        local prompt="$prompt_dir/prompt-${id}.txt"
        [ -f "$prompt" ] || { echo "  WARN: $prompt not found, skipping"; continue; }

        for condition in with-skill without-skill; do
            local outfile="$out_base/$condition/prompt-${id}.md"
            wait_for_slot
            run_docker_claude "$prompt" "$outfile" "$condition" "$CLAUDE_MD_GENERIC" &
            JOB_PIDS+=($!)
        done
    done

    wait_all
    echo "  Suite 1: Complete ($COMPLETED succeeded, $FAILED failed)"
    echo ""
}

# ---------------------------------------------------------------------------
# Suite 2: Statistical Re-Runs (2 diseases x 3 runs x 5 prompts)
# ---------------------------------------------------------------------------
run_suite_reruns() {
    echo "============================================================"
    echo "  Suite 2: Statistical Re-Runs (60 sessions)"
    echo "============================================================"
    local out_base="$EVAL_DIR/reruns/outputs-containerized"

    for disease in polio guinea-worm; do
        local prompt_dir claude_md
        if [ "$disease" = "polio" ]; then
            prompt_dir="$EVAL_DIR"
            claude_md="$CLAUDE_MD_POLIO"
        else
            prompt_dir="$EVAL_DIR/guinea-worm"
            claude_md="$CLAUDE_MD_GUINEA"
        fi

        for run_num in 1 2 3; do
            for prompt_num in 1 2 3 4 5; do
                local prompt="$prompt_dir/prompt-${prompt_num}.txt"
                [ -f "$prompt" ] || continue

                for condition in with-skill without-skill; do
                    local outfile="$out_base/$disease/run-${run_num}/$condition/prompt-${prompt_num}.md"
                    wait_for_slot
                    run_docker_claude "$prompt" "$outfile" "$condition" "$claude_md" &
                    JOB_PIDS+=($!)
                done
            done
        done
    done

    wait_all
    echo "  Suite 2: Complete ($COMPLETED succeeded, $FAILED failed)"
    echo ""
}

# ---------------------------------------------------------------------------
# Suite 3: Difficulty Curve (8 prompts)
# ---------------------------------------------------------------------------
run_suite_difficulty() {
    echo "============================================================"
    echo "  Suite 3: Difficulty Curve (8 prompts x 2 conditions)"
    echo "============================================================"
    local prompt_dir="$EVAL_DIR/difficulty-curve"
    local out_base="$EVAL_DIR/difficulty-curve/outputs-containerized"

    for id in D1 D2 D3 D4 D5 D6 D7 D8; do
        local prompt="$prompt_dir/prompt-${id}.txt"
        [ -f "$prompt" ] || { echo "  WARN: $prompt not found, skipping"; continue; }

        for condition in with-skill without-skill; do
            local outfile="$out_base/$condition/prompt-${id}.md"
            wait_for_slot
            run_docker_claude "$prompt" "$outfile" "$condition" "$CLAUDE_MD_GENERIC" &
            JOB_PIDS+=($!)
        done
    done

    wait_all
    echo "  Suite 3: Complete ($COMPLETED succeeded, $FAILED failed)"
    echo ""
}

# ---------------------------------------------------------------------------
# Suite 4: Multi-Turn (delegates to Python driver with --docker flag)
# ---------------------------------------------------------------------------
run_suite_multiturn() {
    local num_runs="${MULTITURN_RUNS:-3}"
    echo "============================================================"
    echo "  Suite 4: Multi-Turn (5 steps x 2 conditions x ${num_runs} runs)"
    echo "============================================================"
    local driver="$EVAL_DIR/multi-turn/run-multi-turn.py"

    if [ ! -f "$driver" ]; then
        echo "  ERROR: $driver not found"
        return 1
    fi

    # Run both conditions in parallel (steps are sequential within each).
    # Pass the auth token via env so the Python driver can forward it to containers.
    # The driver handles --runs N for statistical re-runs with resume support.
    CLAUDE_CODE_OAUTH_TOKEN="$DOCKER_AUTH_TOKEN" /opt/anaconda3/bin/python3 "$driver" with --docker --runs "$num_runs" &
    local pid_with=$!

    CLAUDE_CODE_OAUTH_TOKEN="$DOCKER_AUTH_TOKEN" /opt/anaconda3/bin/python3 "$driver" without --docker --runs "$num_runs" &
    local pid_without=$!

    wait "$pid_with"  && COMPLETED=$((COMPLETED + 1)) || FAILED=$((FAILED + 1))
    wait "$pid_without" && COMPLETED=$((COMPLETED + 1)) || FAILED=$((FAILED + 1))

    echo "  Suite 4: Complete"
    echo ""
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
echo "============================================================"
echo "  Containerized Eval Runner"
echo "============================================================"
echo "  Suite:       $SUITE"
echo "  Parallelism: $MAX_JOBS"
echo "  Project:     $PROJECT_DIR"
echo "  Auth:        $HOME/.claude/ (mounted read-only)"
echo "============================================================"
echo ""

START_TIME=$(date +%s)

# Build images (shared by all suites)
build_images

case "$SUITE" in
    atomic)
        run_suite_atomic
        ;;
    difficulty)
        run_suite_difficulty
        ;;
    reruns)
        run_suite_reruns
        ;;
    multi-turn)
        run_suite_multiturn
        ;;
    all)
        # Run all 4 suites — suites 1-3 share the job pool, suite 4 is separate
        run_suite_atomic
        COMPLETED=0; FAILED=0
        run_suite_difficulty
        COMPLETED=0; FAILED=0
        run_suite_reruns
        COMPLETED=0; FAILED=0
        run_suite_multiturn
        ;;
    *)
        echo "ERROR: Unknown suite '$SUITE'"
        echo "Valid suites: atomic | difficulty | reruns | multi-turn | all"
        exit 1
        ;;
esac

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo "============================================================"
echo "  ALL COMPLETE"
echo "============================================================"
echo "  Wall time: $(format_duration $ELAPSED)"
echo ""
echo "  Outputs:"
echo "    eval/atomic/outputs-containerized/"
echo "    eval/difficulty-curve/outputs-containerized/"
echo "    eval/reruns/outputs-containerized/"
echo "    eval/multi-turn/outputs-containerized/"
echo ""
echo "  Verification:"
echo "    1. Check all output files are non-empty"
echo "    2. Confirm no laser imports in WITHOUT outputs:"
echo "       grep -r 'laser\\.' eval/*/outputs-containerized/without-skill/ || echo 'Clean'"
echo "    3. Score using existing rubrics"
echo "============================================================"
