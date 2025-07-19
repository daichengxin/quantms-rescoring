FROM ubuntu:22.04

# Some metadata
LABEL base_image="ubuntu:22.04"
LABEL version="1"
LABEL software="quantms-rescoring"
LABEL software.version="0.0.10"
LABEL about.summary="quantms-rescoring: Python scripts and helpers for the quantMS workflow"
LABEL about.home="https://github.com/bigbio/quantms-rescoring"
LABEL about.documentation="https://github.com/bigbio/quantms-rescoring"
LABEL about.license_file="https://github.com/bigbio/quantms-rescoring/blob/main/LICENSE"
LABEL about.tags="Proteomics,MS2PIP,DeepLC,OpenMS"
LABEL maintainer="Yasset Perez-Riverol <ypriverol@gmail.com>"

ENV DEBIAN_FRONTEND=noninteractive

# Set environment variables
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    # pyOpenMS environment variables
    OPENMS_DATA_PATH=/usr/local/lib/python3.11/site-packages/pyopenms/share/OpenMS \
    # Disable CUDA warnings (since we're not using GPU)
    TF_CPP_MIN_LOG_LEVEL=2

# Update package lists and install necessary packages
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
    rm -rf /var/lib/apt/lists/*

# Set work directory
WORKDIR /app

# Configure locale to avoid runtime errors
RUN locale-gen en_US.UTF-8 && \
    update-locale LANG=en_US.UTF-8

# Set environment variables for locale
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US:en
ENV LC_ALL=en_US.UTF-8

# Set work directory
WORKDIR /app

# Copy requirements first for better caching
COPY requirements.txt .

# Install Python dependencies
RUN python3.11 -m pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY . .

# Install the package in development mode
RUN python3.11 -m pip install -e .

# Test pyOpenMS import
RUN python3.11 -c "import pyopenms; print('pyOpenMS imported successfully')"

# Remove unnecessary packages
RUN apt-get remove -y wget && apt-get autoremove -y && apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Create a non-root user
RUN useradd --create-home --shell /bin/bash app && \
    chown -R app:app /app
USER app

WORKDIR /data/

# NOTE: It is entirely the user's responsibility to ensure compliance with quantms-rescoring license terms.
# Please review the licensing terms for quantms-rescoring before using or distributing this Docker image.
