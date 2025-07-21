# ===================
# Stage 1: Build
# ===================
FROM python:3.11-slim as builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    python3.11-dev \
    python3.11-venv \
    libglib2.0-0 \
    libgomp1 \
    libboost-all-dev \
    libhdf5-dev \
    libnetcdf-dev \
    libxml2-dev \
    libxslt1-dev \
    libssl-dev \
    libffi-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . .

RUN pip install --no-cache-dir pip==25.1.1 && \
    pip install --no-cache-dir poetry==2.1.3 && \
    poetry build && \
    pip install --no-cache-dir dist/*.whl && \
    pip cache purge

# ===================
# Stage 2: Runtime
# ===================
FROM python:3.11-slim as final

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US:en
ENV LC_ALL=en_US.UTF-8

RUN apt-get update && apt-get install -y --no-install-recommends \
    libglib2.0-0 \
    libgomp1 \
    libhdf5-103-1 \
    libnetcdf19 \
    libxml2 \
    libxslt1.1 \
    libssl3 \
    libffi8 && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/local /usr/local

RUN useradd --create-home --shell /bin/bash app && \
    mkdir /data && chown -R app:app /data
USER app

WORKDIR /data

RUN python3.11 -c "import pyopenms; print('pyOpenMS imported successfully')"
