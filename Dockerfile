# GNPS Reference Library Annotation Pipeline Docker Image
# Version: 3.0

FROM python:3.10-slim

LABEL maintainer="Bigy Ambat"
LABEL description="Docker image for GNPS Reference Library Annotation Pipeline"
LABEL version="3.0"

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
RUN pip install --no-cache-dir --upgrade pip

# Install all required Python packages
RUN pip install --no-cache-dir \
    pandas>=1.5.0 \
    matplotlib>=3.6.0 \
    seaborn>=0.12.0 \
    pyteomics>=4.5.0 \
    numpy>=1.23.0 \
    scipy>=1.9.0 \
    rdkit>=2022.9.1

# Verify installations
RUN python -c "import pandas; print('Pandas version:', pandas.__version__)" && \
    python -c "import matplotlib; print('Matplotlib version:', matplotlib.__version__)" && \
    python -c "import pyteomics.mgf; print('Pyteomics: OK')" && \
    python -c "import numpy; print('NumPy version:', numpy.__version__)" && \
    python -c "import scipy; print('SciPy version:', scipy.__version__)" && \
    python -c "import seaborn; print('Seaborn version:', seaborn.__version__)" && \
    python -c "from rdkit import Chem; print('RDKit: OK')"

# Copy GNPS pipeline scripts
COPY bin/parse_gnps_reference.py /usr/local/bin/
COPY bin/annotate_peaks_gnps.py /usr/local/bin/
COPY bin/calculate_cosine_similarity.py /usr/local/bin/
COPY bin/generate_qc_report_gnps.py /usr/local/bin/

# Make scripts executable
RUN chmod +x /usr/local/bin/parse_gnps_reference.py && \
    chmod +x /usr/local/bin/annotate_peaks_gnps.py && \
    chmod +x /usr/local/bin/calculate_cosine_similarity.py && \
    chmod +x /usr/local/bin/generate_qc_report_gnps.py

# Add /usr/local/bin to PATH
ENV PATH="/usr/local/bin:${PATH}"

# Set working directory for pipeline execution
WORKDIR /work

# Default command
CMD ["/bin/bash"]
