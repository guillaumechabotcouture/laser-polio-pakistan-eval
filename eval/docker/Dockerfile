# eval/docker/Dockerfile
# Multi-stage build for containerized eval runner
#
# Targets:
#   without-laser — Base environment, NO laser packages (for WITHOUT condition)
#   with-laser    — Base + laser-generic installed (for WITH condition)
#
# Build:
#   docker build --target without-laser -t laser-eval-without -f eval/docker/Dockerfile eval/docker/
#   docker build --target with-laser    -t laser-eval-with    -f eval/docker/Dockerfile eval/docker/

FROM node:22-slim AS base

# Python 3 + system deps for scientific packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-venv \
    python3-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Scientific Python stack
RUN pip3 install --break-system-packages --no-cache-dir \
    numpy \
    geopandas \
    shapely \
    matplotlib

# Claude Code CLI
RUN npm install -g @anthropic-ai/claude-code

# Non-root user (--dangerously-skip-permissions refuses to run as root)
RUN useradd -m -s /bin/bash evaluser

# -------------------------------------------------------------------
# WITHOUT condition: no laser packages
# -------------------------------------------------------------------
FROM base AS without-laser
USER evaluser
WORKDIR /home/evaluser

# -------------------------------------------------------------------
# WITH condition: laser-core + laser-generic installed
# -------------------------------------------------------------------
FROM base AS with-laser

# laser-core 1.0.0 has no Linux wheels on PyPI (macOS ARM64 only).
# Install build tools, then build from GitHub source (still as root).
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install laser-core from GitHub (compiles C extension), then laser-generic from PyPI
RUN pip3 install --break-system-packages --no-cache-dir \
    "laser-core @ git+https://github.com/InstituteforDiseaseModeling/laser.git@v1.0.0" \
    laser-generic==1.0.0

# Switch to non-root user after all installs
USER evaluser
WORKDIR /home/evaluser
