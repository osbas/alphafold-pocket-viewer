# ==============================================================================
# Dockerfile — AlphaFold-Lite Pocket Viewer (Sistemas Programables)
# ==============================================================================
FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    DIRECTORIO_CACHE_PDB=/app/datos/estructuras_pdb \
    PORT=8001

WORKDIR /app

# 1. Instalar dependencias aprovechando el caché de capas de Docker
COPY requirements.txt .
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# 2. Crear usuario sin privilegios (appuser, UID 1000) y directorio de volumen
RUN groupadd --gid 1000 appgroup && \
    useradd --uid 1000 --gid appgroup --shell /bin/bash --create-home appuser && \
    mkdir -p /app/datos/estructuras_pdb && \
    chown -R appuser:appgroup /app

# 3. Copiar código fuente y pruebas
COPY --chown=appuser:appgroup app/ ./app/
COPY --chown=appuser:appgroup tests/ ./tests/

USER appuser

EXPOSE 8001

HEALTHCHECK --interval=15s --timeout=5s --start-period=5s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8001/api/salud', timeout=3)" || exit 1

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8001"]
