#!/usr/bin/env python3
"""Multi-turn iterative evaluation driver for LASER skill A/B testing.

Runs a 5-step iterative conversation where each step builds on the previous.
Executes both WITH and WITHOUT skill conditions.

Each step:
  1. Sends a prompt (with previous code as context) to Claude
  2. Extracts Python code from the response
  3. Attempts to execute the extracted code
  4. Uses the code + execution result to inform the next step

Step 4 is dynamic: its prompt depends on whether Step 3's code runs or errors.

Usage:
    python eval/multi-turn/run-multi-turn.py              # Run both conditions
    python eval/multi-turn/run-multi-turn.py with         # WITH skill only
    python eval/multi-turn/run-multi-turn.py without      # WITHOUT skill only
    python eval/multi-turn/run-multi-turn.py both         # Both (default)
    python eval/multi-turn/run-multi-turn.py with --docker    # WITH, containerized
    python eval/multi-turn/run-multi-turn.py both --docker    # Both, containerized
    python eval/multi-turn/run-multi-turn.py both --docker --runs 3  # 3 runs each
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile
import textwrap
import time
import uuid
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.resolve()
PROJECT_DIR = SCRIPT_DIR.parent.parent

# ---------------------------------------------------------------------------
# Global mode flags (set by CLI parsing)
# ---------------------------------------------------------------------------
USE_DOCKER = False
NUM_RUNS = 1

# Docker image names (must match eval/docker/run-all-containerized.sh)
IMG_WITH = "laser-eval-with"
IMG_WITHOUT = "laser-eval-without"

# Timeouts
CLAUDE_TIMEOUT_DIRECT = 300   # seconds for direct CLI
CLAUDE_TIMEOUT_DOCKER = 1200  # seconds for Docker (increased from 600)
CODE_EXEC_TIMEOUT = 120       # seconds for running extracted code

# ---------------------------------------------------------------------------
# Inline code instruction appended to every prompt
# ---------------------------------------------------------------------------
INLINE_CODE_INSTRUCTION = (
    "\n\nIMPORTANT: Output the complete Python code in a single ```python code "
    "block in your response. Do not write it to a file — include it directly "
    "in your answer."
)

# ---------------------------------------------------------------------------
# Step prompts
# ---------------------------------------------------------------------------

STEP_PROMPTS = {
    1: textwrap.dedent("""\
        Using the LASER framework (laser-generic package), build a basic 4-patch SEIR
        model for a respiratory disease. Patches have populations 100k, 200k, 150k,
        80k. Disease: R0~5, latent period 4 days, infectious period 10 days. Initialize
        with 90% susceptible, 1% infectious, 9% recovered. Run for 1 year.

        IMPORTANT: Build ONLY a basic SEIR model for this step. Do NOT add spatial
        coupling (no gravity model), seasonal forcing, or demographics yet — those
        will be added in later steps.

        Write complete Python code. Do not install any packages."""),

    2: textwrap.dedent("""\
        Add gravity-model spatial coupling to this model. The 4 patches are arranged
        in a line 75km apart. Use gravity parameters k=0.01, a=1, b=1, c=1.5.
        Row-normalize so no patch exports more than 15%. Update the code and show the
        complete modified script. Do not install any packages."""),

    3: textwrap.dedent("""\
        Add seasonal forcing to the transmission. Winter peak (days 0-90) at 1.3x
        baseline, summer trough (days 150-240) at 0.7x baseline. Use LASER's ValuesMap
        to create the seasonal profile. Update the code and show the complete modified
        script. Do not install any packages."""),

    # Step 4 is generated dynamically by make_step4_prompt()

    5: textwrap.dedent("""\
        Add births (CBR=30 per 1000/year) and deaths (CDR=10 per 1000/year) using
        LASER's BirthsByCBR and MortalityByCDR. Use calc_capacity to pre-allocate for
        10 years. Extend the simulation to 10 years. After running, print the total
        population at the start and end -- does the ~2% annual growth rate look
        correct? Show the complete script. Do not install any packages."""),
}

# ---------------------------------------------------------------------------
# Minimal CLAUDE.md for the WITHOUT-skill condition
# ---------------------------------------------------------------------------

WITHOUT_SKILL_CLAUDE_MD = textwrap.dedent("""\
    # Project Context
    This is a public health epidemiological modeling project using the
    open-source LASER framework (https://laser.idmod.org/) from the
    Institute for Disease Modeling. The project builds spatial disease
    transmission models for public health planning.
""")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def extract_python_code(response: str) -> str:
    """Extract the largest Python code block from Claude's response.

    Looks for ```python ... ``` blocks first, then bare ``` ... ``` blocks.
    Returns the longest match, or empty string if none found.
    """
    # Try ```python blocks first
    blocks = re.findall(r"```python\n(.*?)```", response, re.DOTALL)
    if not blocks:
        # Fall back to bare ``` blocks
        blocks = re.findall(r"```\n(.*?)```", response, re.DOTALL)
    if not blocks:
        return ""
    return max(blocks, key=len).strip()


def try_run_code(code: str, timeout: int = CODE_EXEC_TIMEOUT) -> tuple:
    """Attempt to run extracted Python code.

    Returns:
        (success: bool, output_or_error: str)
    """
    if not code.strip():
        return False, "NO_CODE: No Python code was extracted from the response"

    tmp = tempfile.NamedTemporaryFile(
        suffix=".py", mode="w", delete=False, prefix="laser_eval_"
    )
    try:
        tmp.write(code)
        tmp.flush()
        tmp.close()

        result = subprocess.run(
            ["/opt/anaconda3/bin/python3", tmp.name],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        if result.returncode == 0:
            output = result.stdout[:2000]
            if result.stderr:
                output += "\n--- stderr ---\n" + result.stderr[:1000]
            return True, output
        else:
            error = result.stderr[:2000]
            if result.stdout:
                error = "--- stdout ---\n" + result.stdout[:500] + "\n--- stderr ---\n" + error
            return False, error

    except subprocess.TimeoutExpired:
        return False, f"TIMEOUT: Script exceeded {timeout}-second limit"
    except Exception as exc:
        return False, f"EXECUTION_ERROR: {exc}"
    finally:
        try:
            os.unlink(tmp.name)
        except OSError:
            pass


def make_step4_prompt(step3_code: str, step3_ran: bool, step3_error: str) -> str:
    """Generate the Step 4 prompt based on Step 3 execution results.

    If Step 3 errored: ask to fix the actual error.
    If Step 3 ran: pose a debugging scenario about identical curves.
    """
    if not step3_ran:
        return (
            "The model crashes with the following error:\n"
            "\n"
            f"{step3_error}\n"
            "\n"
            "Fix the issue and show the complete corrected script. "
            "Do not install any packages."
        )
    else:
        return (
            "The model runs but all patches show identical infection curves, "
            "suggesting spatial coupling isn't working. Debug and fix the "
            "issue. Show the complete corrected script. "
            "Do not install any packages."
        )


def build_full_prompt(step_num: int, prompt: str, prev_code: str) -> str:
    """Build the full prompt with previous code context for steps 2+."""
    # Append inline code instruction to every prompt
    prompt_with_inline = prompt + INLINE_CODE_INSTRUCTION

    if step_num == 1 or not prev_code.strip():
        return prompt_with_inline

    return (
        "Here is the code from the previous step:\n"
        "\n"
        "```python\n"
        f"{prev_code}\n"
        "```\n"
        "\n"
        f"{prompt_with_inline}"
    )


# ---------------------------------------------------------------------------
# Claude invocation — direct (original) and Docker (containerized)
# ---------------------------------------------------------------------------


def call_claude(
    full_prompt: str, condition: str, tmpdir: str | None = None
) -> str:
    """Invoke Claude CLI and return the response text.

    Dispatches to Docker or direct invocation based on USE_DOCKER flag.
    Includes retry on timeout (one retry attempt).

    Args:
        full_prompt: The complete prompt to send.
        condition: "with-skill" or "without-skill".
        tmpdir: Working directory for without-skill (created by caller).

    Returns:
        Claude's response text.
    """
    for attempt in range(2):  # retry once on timeout
        try:
            if USE_DOCKER:
                return _call_claude_docker(full_prompt, condition)
            return _call_claude_direct(full_prompt, condition, tmpdir)
        except subprocess.TimeoutExpired:
            timeout = CLAUDE_TIMEOUT_DOCKER if USE_DOCKER else CLAUDE_TIMEOUT_DIRECT
            if attempt == 0:
                print(f"    TIMEOUT after {timeout}s, retrying...", flush=True)
                continue
            return f"ERROR: Claude call timed out after {timeout} seconds (2 attempts)."
    return "ERROR: Unexpected retry loop exit."


def _call_claude_direct(
    full_prompt: str, condition: str, tmpdir: str | None = None
) -> str:
    """Invoke Claude CLI directly (original behavior)."""
    env = os.environ.copy()
    env.pop("CLAUDECODE", None)

    if condition == "with-skill":
        cmd = ["claude", "--print", "--dangerously-skip-permissions"]
        cwd = str(PROJECT_DIR)
    else:
        cmd = [
            "claude",
            "--print",
            "--dangerously-skip-permissions",
            "--disable-slash-commands",
        ]
        cwd = tmpdir

    result = subprocess.run(
        cmd,
        input=full_prompt,
        capture_output=True,
        text=True,
        cwd=cwd,
        env=env,
        timeout=CLAUDE_TIMEOUT_DIRECT,
    )

    return result.stdout


def _get_docker_auth_token() -> str:
    """Read OAuth token from env or ~/.claude/oauth-token file."""
    token = os.environ.get("CLAUDE_CODE_OAUTH_TOKEN", "")
    if token:
        return token
    token_file = Path.home() / ".claude" / "oauth-token"
    if token_file.exists():
        return token_file.read_text(encoding="utf-8").strip()
    raise RuntimeError(
        "No auth token found. Set CLAUDE_CODE_OAUTH_TOKEN or create "
        f"{token_file} (generate with: claude setup-token)"
    )


def _call_claude_docker(full_prompt: str, condition: str) -> str:
    """Invoke Claude CLI inside a Docker container.

    WITH condition:  project mounted at /workspace, laser-generic installed.
    WITHOUT condition: clean /eval-work dir with minimal CLAUDE.md, no laser.

    Auth via CLAUDE_CODE_OAUTH_TOKEN env var passed to container.
    Uses named containers so we can explicitly kill them on timeout.
    """
    token = _get_docker_auth_token()
    container_name = f"laser-eval-{uuid.uuid4().hex[:8]}"
    docker_tmpdir = None

    if condition == "with-skill":
        cmd = [
            "docker", "run", "--rm", "-i",
            "--name", container_name,
            "-e", f"CLAUDE_CODE_OAUTH_TOKEN={token}",
            "-v", f"{PROJECT_DIR}:/workspace:ro",
            IMG_WITH,
            "sh", "-c",
            "cd /workspace && claude --print --dangerously-skip-permissions",
        ]
    else:
        docker_tmpdir = tempfile.mkdtemp(prefix="laser_eval_docker_")
        Path(docker_tmpdir, "CLAUDE.md").write_text(
            WITHOUT_SKILL_CLAUDE_MD, encoding="utf-8"
        )
        cmd = [
            "docker", "run", "--rm", "-i",
            "--name", container_name,
            "-e", f"CLAUDE_CODE_OAUTH_TOKEN={token}",
            "-v", f"{docker_tmpdir}:/eval-work:ro",
            IMG_WITHOUT,
            "sh", "-c",
            "cd /eval-work && claude --print --dangerously-skip-permissions --disable-slash-commands",
        ]

    try:
        result = subprocess.run(
            cmd,
            input=full_prompt,
            capture_output=True,
            text=True,
            timeout=CLAUDE_TIMEOUT_DOCKER,
        )
        return result.stdout
    except subprocess.TimeoutExpired:
        # Kill the Docker container explicitly — subprocess.run only kills
        # the docker client process, leaving the container running as a zombie.
        subprocess.run(
            ["docker", "kill", container_name],
            capture_output=True, timeout=10,
        )
        subprocess.run(
            ["docker", "rm", "-f", container_name],
            capture_output=True, timeout=10,
        )
        raise  # re-raise so call_claude() retry logic handles it
    finally:
        if docker_tmpdir:
            shutil.rmtree(docker_tmpdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------------------------


def run_step(
    step_num: int,
    prompt: str,
    condition: str,
    prev_code: str,
    output_dir: Path,
) -> tuple:
    """Run one step of the multi-turn conversation.

    Args:
        step_num: Step number (1-5).
        prompt: The step's prompt text (not yet wrapped with prev_code).
        condition: "with-skill" or "without-skill".
        prev_code: Python code from the previous step.
        output_dir: Directory to save outputs.

    Returns:
        (extracted_code: str, ran: bool, output_or_error: str)
    """
    # Build full prompt with context
    full_prompt = build_full_prompt(step_num, prompt, prev_code)

    # Save prompt
    prompt_file = output_dir / f"step-{step_num}-prompt.txt"
    prompt_file.write_text(full_prompt, encoding="utf-8")

    # Run Claude
    tmpdir = None
    try:
        if not USE_DOCKER and condition == "without-skill":
            tmpdir = tempfile.mkdtemp(prefix="laser_eval_ws_")
            Path(tmpdir, "CLAUDE.md").write_text(
                WITHOUT_SKILL_CLAUDE_MD, encoding="utf-8"
            )

        print(f"    Calling Claude for step {step_num}...", flush=True)
        t0 = time.time()

        response = call_claude(full_prompt, condition, tmpdir)

        elapsed = time.time() - t0
        print(f"    Claude responded in {elapsed:.0f}s ({len(response)} chars)")

    finally:
        if tmpdir is not None:
            shutil.rmtree(tmpdir, ignore_errors=True)

    # Save response
    response_file = output_dir / f"step-{step_num}-response.md"
    response_file.write_text(response, encoding="utf-8")

    # Extract code
    code = extract_python_code(response)
    code_file = output_dir / f"step-{step_num}-code.py"
    code_file.write_text(code if code else "# No code extracted\n", encoding="utf-8")

    # Try running the code
    ran, output = try_run_code(code)

    # Save execution log
    status = "SUCCESS" if ran else "FAILED"
    log_file = output_dir / f"step-{step_num}-execution.log"
    log_file.write_text(
        f"Status: {status}\n"
        f"Code length: {len(code)} chars\n"
        f"---\n"
        f"{output}\n",
        encoding="utf-8",
    )

    return code, ran, output


def run_condition(condition: str, run_label: str = "") -> dict:
    """Run all 5 steps for one condition (with-skill or without-skill).

    Args:
        condition: "with-skill" or "without-skill".
        run_label: Optional subdirectory label (e.g., "run-1") for multi-run mode.

    Returns a dict of per-step results for the summary.
    """
    # Use outputs-containerized/ when in Docker mode
    output_subdir = "outputs-containerized" if USE_DOCKER else "outputs"
    if run_label:
        output_dir = SCRIPT_DIR / output_subdir / run_label / condition
    else:
        output_dir = SCRIPT_DIR / output_subdir / condition
    output_dir.mkdir(parents=True, exist_ok=True)

    # Resume support: skip if all 5 response files exist and are non-trivial
    existing_responses = [
        output_dir / f"step-{s}-response.md" for s in range(1, 6)
    ]
    if all(f.exists() and f.stat().st_size > 100 for f in existing_responses):
        print(f"  SKIP (all 5 steps exist): {output_dir}", flush=True)
        # Reconstruct results from saved files
        results = {}
        for step in range(1, 6):
            code_file = output_dir / f"step-{step}-code.py"
            log_file = output_dir / f"step-{step}-execution.log"
            code_len = 0
            ran = False
            if code_file.exists():
                code_len = len(code_file.read_text(encoding="utf-8").strip())
            if log_file.exists():
                ran = log_file.read_text(encoding="utf-8").startswith("Status: SUCCESS")
            results[step] = {
                "ran": ran,
                "code_len": code_len,
                "output_preview": "(resumed from saved)",
            }
        return results

    results = {}
    prev_code = ""
    step3_ran = False
    step3_error = ""

    for step in [1, 2, 3, 4, 5]:
        print(f"  Step {step}/5:", flush=True)

        # Get the prompt
        if step == 4:
            prompt = make_step4_prompt(prev_code, step3_ran, step3_error)
        else:
            prompt = STEP_PROMPTS[step]

        # Run the step
        code, ran, output = run_step(
            step, prompt, condition, prev_code, output_dir
        )

        # Track Step 3 results for dynamic Step 4
        if step == 3:
            step3_ran = ran
            step3_error = output

        # Use extracted code as context for next step (even if it didn't run,
        # since the next step should try to fix or build on it)
        if code.strip():
            prev_code = code

        status_str = "RAN" if ran else "FAILED"
        print(f"    Result: {status_str} | Code: {len(code)} chars")

        results[step] = {
            "ran": ran,
            "code_len": len(code),
            "output_preview": output[:200],
        }

    return results


def print_summary(all_results: dict) -> None:
    """Print a summary table of all results."""
    mode_label = " (containerized)" if USE_DOCKER else ""
    output_subdir = "outputs-containerized" if USE_DOCKER else "outputs"

    print("\n" + "=" * 70)
    print(f"SUMMARY{mode_label}")
    print("=" * 70)

    header = f"{'Condition':<24} {'Step':>4} {'Runs?':>6} {'Code Len':>9}"
    print(header)
    print("-" * len(header))

    for condition, steps in all_results.items():
        for step_num in sorted(steps.keys()):
            info = steps[step_num]
            ran_str = "YES" if info["ran"] else "NO"
            print(
                f"{condition:<24} {step_num:>4} {ran_str:>6} {info['code_len']:>9}"
            )
        print()

    # Count runs per condition
    for condition, steps in all_results.items():
        ran_count = sum(1 for s in steps.values() if s["ran"])
        print(f"{condition}: {ran_count}/5 steps produced running code")

    print()
    print(f"Outputs saved to: {SCRIPT_DIR / output_subdir}")
    print(f"Score using the rubric: {SCRIPT_DIR / 'rubric.md'}")


def main():
    global USE_DOCKER, NUM_RUNS

    # Parse arguments: [with|without|both] [--docker] [--runs N]
    args = sys.argv[1:]
    target = "both"
    i = 0
    while i < len(args):
        arg = args[i]
        if arg == "--docker":
            USE_DOCKER = True
        elif arg == "--runs":
            i += 1
            if i < len(args):
                NUM_RUNS = int(args[i])
            else:
                print("ERROR: --runs requires a number")
                sys.exit(1)
        elif arg in ("with", "without", "both"):
            target = arg
        else:
            print("Usage: python run-multi-turn.py [with|without|both] [--docker] [--runs N]")
            print("  with         Run WITH-skill condition only")
            print("  without      Run WITHOUT-skill condition only")
            print("  both         Run both conditions (default)")
            print("  --docker     Use Docker containers for Claude calls")
            print("               (outputs go to outputs-containerized/)")
            print("  --runs N     Run N independent runs per condition (default 1)")
            print("               Outputs go to run-1/, run-2/, ... subdirectories")
            sys.exit(1)
        i += 1

    mode_label = " [Docker]" if USE_DOCKER else ""
    all_results = {}

    for run_num in range(1, NUM_RUNS + 1):
        run_label = f"run-{run_num}" if NUM_RUNS > 1 else ""
        run_prefix = f"[Run {run_num}/{NUM_RUNS}] " if NUM_RUNS > 1 else ""

        if target in ("with", "both"):
            label = f"{run_prefix}with-skill"
            print("=" * 70)
            print(f"CONDITION: {run_prefix}WITH SKILL{mode_label}")
            print("=" * 70)
            all_results[label] = run_condition("with-skill", run_label)

        if target in ("without", "both"):
            label = f"{run_prefix}without-skill"
            print("=" * 70)
            print(f"CONDITION: {run_prefix}WITHOUT SKILL{mode_label}")
            print("=" * 70)
            all_results[label] = run_condition("without-skill", run_label)

    print_summary(all_results)


if __name__ == "__main__":
    main()
