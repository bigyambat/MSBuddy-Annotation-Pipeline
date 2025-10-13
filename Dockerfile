# MS Annotation & QC Pipeline Docker Image
# Version: 2.0

FROM python:3.10-slim

LABEL maintainer="Bigy Ambat"
LABEL description="Docker image for MS Annotation & QC Pipeline with MSBuddy"
LABEL version="2.0"

# Set environment variables
ENV PYTHONUNBUFFERED=1
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    wget \
    git \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Create app directory
WORKDIR /app

# Install Python dependencies
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir \
    pandas>=1.5.0 \
    matplotlib>=3.6.0 \
    seaborn>=0.12.0 \
    pyteomics>=4.5.0 \
    numpy>=1.23.0 \
    scipy>=1.9.0

# Install MSBuddy
# Note: Adjust this based on actual MSBuddy installation method
# If MSBuddy is available via pip:
RUN pip install --no-cache-dir msbuddy || \
    echo "MSBuddy installation via pip failed. Please install manually or from source."

# Create MSBuddy data directory with write permissions
# This fixes the "Permission denied" error when MSBuddy tries to initialize its database
RUN mkdir -p /usr/local/lib/python3.10/site-packages/msbuddy/data && \
    chmod -R 777 /usr/local/lib/python3.10/site-packages/msbuddy/data

# Alternative: Install MSBuddy from GitHub if not available on PyPI
# Uncomment and adjust as needed:
# RUN git clone https://github.com/Philipbear/msbuddy.git && \
#     cd msbuddy && \
#     pip install -e . && \
#     cd .. && \
#     rm -rf msbuddy

# Copy pipeline scripts
COPY bin/generate_qc_report.py /usr/local/bin/
RUN chmod +x /usr/local/bin/generate_qc_report.py

# Add /usr/local/bin to PATH
ENV PATH="/usr/local/bin:${PATH}"

# Set working directory for pipeline execution
WORKDIR /work

# Default command
CMD ["/bin/bash"]
