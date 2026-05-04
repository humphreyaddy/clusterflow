# Phase 9.2 — multi-stage production image. See transmission_cluster_tool_plan.md §9.2.
# syntax=docker/dockerfile:1.6

# ----- builder stage -----
FROM python:3.11-slim AS builder

ENV PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# Build deps for igraph / leidenalg / numpy wheels
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        libxml2-dev \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY pyproject.toml README.md ./
COPY clusterflow ./clusterflow

# Install into a staging prefix; include the streaming + dashboard extras so
# the runtime image can serve. R-bridge is intentionally excluded.
RUN pip install --prefix=/install ".[dashboard]" || pip install --prefix=/install .


# ----- runtime stage -----
FROM python:3.11-slim AS runtime

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PATH="/usr/local/bin:${PATH}"

# Slim runtime libs
RUN apt-get update && apt-get install -y --no-install-recommends \
        libxml2 \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/* \
    && useradd --create-home --uid 1000 clusterflow

COPY --from=builder /install /usr/local

USER clusterflow
WORKDIR /home/clusterflow

EXPOSE 8000 8050

# Health check: validate the CLI loads
HEALTHCHECK --interval=30s --timeout=5s --start-period=10s \
    CMD clusterflow --help > /dev/null || exit 1

ENTRYPOINT ["clusterflow"]
CMD ["--help"]
