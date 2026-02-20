# LASER Skill A/B Evaluation: Polio Model for Pakistan

Evaluates whether a Claude Code skill for the LASER disease modeling framework
improves code generation quality when building a polio transmission model — a
task that's deliberately out-of-domain from the skill's measles/England&Wales
training scenario.

## Quick Start

```bash
cd laser-polio-pakistan-eval

# Run a single prompt (both conditions)
./eval/run-eval.sh 1

# Run all 5 prompts
./eval/run-eval.sh all

# Or run manually:
# WITH skill:
claude
# WITHOUT skill:
claude --deny "Skill(laser-spatial-disease-modeling)"
```

## Structure

```
.claude/skills/laser-spatial-disease-modeling/   # The skill under test
eval/
  prompts.md           # 5 test prompts (basic setup → full integration)
  rubric.md            # Scoring rubric (API, structure, domain, completeness)
  run-eval.sh          # Automated runner
  outputs/
    with-skill/        # Captured outputs from skill-enabled sessions
    without-skill/     # Captured outputs from skill-disabled sessions
reference/
  polio-model-requirements.md   # Ground truth for polio model features
```

## What This Tests

| Dimension | What the skill should help | Risk of negative transfer |
|-----------|---------------------------|--------------------------|
| API correctness | Import paths, function signatures, component ordering | Low |
| Structural correctness | GeoDataFrame setup, network init, model assembly | Low |
| Domain adaptation | Calibration patterns, spatial analysis | Skill may bias toward measles patterns |
| Completeness | Code templates, working examples | May produce rigid measles-shaped code |

## Scoring

See `eval/rubric.md` for the full rubric. Each prompt is scored 0-3 on four
dimensions (max 12 per prompt, max 60 total across all 5 prompts).
