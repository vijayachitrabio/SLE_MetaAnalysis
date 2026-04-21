# SLE GWAS Meta-Analysis Pipeline - Docker Container
# Complete reproducible environment for R and Python pipeline scripts

FROM rocker/tidyverse:4.3.0

# Set working directory
WORKDIR /sle_analysis

# Install system dependencies and Python 3
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements file and install Python packages
# (Using --break-system-packages is necessary on recent Debian/Ubuntu versions for pip install)
COPY requirements.txt /sle_analysis/
RUN pip3 install --no-cache-dir --break-system-packages -r /sle_analysis/requirements.txt

# Install core R packages from CRAN
RUN Rscript -e "install.packages(c('vroom', 'data.table', 'cowplot', 'ggpubr', 'pheatmap', 'viridis', 'msigdbr', 'gprofiler2', 'rjson'), repos = 'http://cran.us.r-project.org')"

# Install GWAS-specific R packages from CRAN / GitHub / Bioconductor if needed
RUN Rscript -e "install.packages('httr')"
RUN Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install(c('coloc', 'biomaRt', 'fgsea'))"

# Copy analysis scripts
COPY scripts/ /sle_analysis/scripts/

# Copy data files (if available during build)
COPY SLE_Publication_Package /sle_analysis/SLE_Publication_Package 2>/dev/null || true

# Set environment variables
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1
ENV MPLBACKEND=Agg

# Create a startup script
RUN echo '#!/bin/bash\n\
echo "================================"\n\
echo "SLE GWAS Meta-Analysis Pipeline"\n\
echo "================================"\n\
echo ""\n\
echo "This container has both R (Tidyverse) and Python available."\n\
echo "Data location: /sle_analysis/data/"\n\
echo "Scripts location: /sle_analysis/scripts/"\n\
echo "Results location: /sle_analysis/results/"\n\
echo "Figures location: /sle_analysis/figures/"\n\
echo ""\n\
if [ "$#" -eq 0 ]; then\n\
  /bin/bash\n\
else\n\
  exec "$@"\n\
fi' > /entrypoint.sh && chmod +x /entrypoint.sh

# Set entrypoint
ENTRYPOINT ["/entrypoint.sh"]

# Default command
CMD ["bash"]

# Labels
LABEL maintainer="SLE GWAS Pipeline"
LABEL version="1.1"
LABEL description="Dual Python/R Docker container for SLE GWAS meta-analysis and therapeutic mapping"
