FROM ubuntu:22.04

# Some metadata
LABEL base_image="ubuntu:22.04"
LABEL version="1"
LABEL software="quantms-rescoring"
LABEL software.version="0.0.11"
LABEL about.summary="quantms-rescoring: Python scripts and helpers for the quantMS workflow"
LABEL about.home="https://github.com/bigbio/quantms-rescoring"
LABEL about.documentation="https://github.com/bigbio/quantms-rescoring"
LABEL about.license_file="https://github.com/bigbio/quantms-rescoring/blob/main/LICENSE"
LABEL about.tags="Proteomics,MS2PIP,DeepLC,OpenMS"
LABEL maintainer="Yasset Perez-Riverol <ypriverol@gmail.com>"

ENV DEBIAN_FRONTEND=noninteractive

# Update package lists and install necessary packages in a single RUN command
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3.11 \
    python3.11-dev \
    python3.11-venv \
    python3-pip \
    build-essential \
    curl \
    git \
    wget \
    locales \
    # pyOpenMS dependencies
    libglib2.0-0 \
    libgomp1 \
    # Additional dependencies for proteomics tools
    libboost-all-dev \
    libhdf5-dev \
    libnetcdf-dev \
    libxml2-dev \
    libxslt1-dev \
    libssl-dev \
    libffi-dev && \
    # Configure locale
    locale-gen en_US.UTF-8 && \
    update-locale LANG=en_US.UTF-8 && \
    # Clean up package cache
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set environment variables for locale
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US:en
ENV LC_ALL=en_US.UTF-8

# Set work directory
WORKDIR /app

# Copy application code first
COPY . .

# Install Poetry and dependencies using the same strategy as CI/CD
RUN python3.11 -m pip install --upgrade pip && \
    python3.11 -m pip install flake8==7.3.0 pytest==8.4.1 && \
    python3.11 -m pip install -r requirements.txt && \
    python3.11 -m pip install poetry==2.1.3 && \
    poetry build && \
    python3.11 -m pip install dist/*.whl && \
    python3.11 -m pip cache purge

# Test pyOpenMS import
RUN python3.11 -c "import pyopenms; print('pyOpenMS imported successfully')"

# Create a non-root user
RUN useradd --create-home --shell /bin/bash app && \
    chown -R app:app /app
USER app

WORKDIR /data/
