# syntax=docker/dockerfile:1
# Multi-stage build for pathway-activity backend
# - Base: python:3.11-slim (hard-pinned minor)
# - Preserves fail-closed release gate (gate defaults OFF; curated endpoints 403 without env)
# - Entrypoint is the real FastAPI app: backend.main:app via uvicorn

FROM python:3.11-slim AS builder

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /app

# System deps (if scipy/numpy need compilation; keep minimal - wheels usually available)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
  && rm -rf /var/lib/apt/lists/*

# Install Python deps first for layer caching
COPY requirements.txt .
RUN pip install --no-cache-dir --prefix=/install -r requirements.txt

# --- runtime stage ---
FROM python:3.11-slim AS runtime

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /app

# Copy installed deps from builder
COPY --from=builder /install /usr/local

# Copy application code
# backend/ contains FastAPI app, data_pipeline/ contains scoring logic
COPY backend/ ./backend/
COPY data_pipeline/ ./data_pipeline/
COPY requirements.txt .

# Create non-root user
RUN useradd -m -u 1000 appuser && chown -R appuser:appuser /app
USER appuser

EXPOSE 8000

# Healthcheck uses the real /health endpoint (honest gate status)
HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
  CMD python -c "import urllib.request,sys; sys.exit(0 if urllib.request.urlopen('http://localhost:8000/health', timeout=4).status==200 else 1)"

# NOTE: release gate is fail-closed by design. Do NOT set MODEL_RELEASE_APPROVED=true here.
# Curated artifact endpoints (/api/v1/curated/*) will return 403 until the deployer
# explicitly sets MODEL_RELEASE_APPROVED=true and APPROVED_ARTIFACT_REVISION.
# Scoring endpoints (/api/v1/score/*) remain available without the gate.
CMD ["uvicorn", "backend.main:app", "--host", "0.0.0.0", "--port", "8000"]
