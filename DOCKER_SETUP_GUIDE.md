# SLE GWAS Pipeline - Docker Setup Guide

Complete Docker containerization for reproducible analysis across all systems.

**Status:** ✅ Production-Ready  
**Version:** 1.0  
**Date:** April 19, 2026

---

## 📋 Table of Contents

1. [Requirements](#requirements)
2. [Quick Start](#quick-start)
3. [Detailed Setup](#detailed-setup)
4. [Running Analysis Scripts](#running-analysis-scripts)
5. [Volume Mounting](#volume-mounting)
6. [Troubleshooting](#troubleshooting)
7. [Advanced Configuration](#advanced-configuration)

---

## Requirements

### System Requirements
- **Docker**: 20.10+ ([Install Docker](https://docs.docker.com/get-docker/))
- **Docker Compose**: 1.29+ ([Install Docker Compose](https://docs.docker.com/compose/install/))
- **Disk Space**: ~2GB for Docker image + data
- **RAM**: Minimum 4GB (8GB recommended)
- **OS**: Linux, macOS, Windows (with WSL2)

### Verify Installation
```bash
docker --version
docker-compose --version
```

---

## Quick Start

### 1. **Clone/Download the Repository**
```bash
cd /path/to/sle-analysis
```

### 2. **Build the Docker Image**
```bash
docker build -t sle-gwas-pipeline:latest .
```

Or using docker-compose:
```bash
docker-compose build
```

### 3. **Run the Container**

**Option A: Using Docker Compose (Recommended)**
```bash
docker-compose up -d
docker-compose exec sle-analysis bash
```

**Option B: Using Docker Directly**
```bash
docker run -it \
  -v $(pwd)/publication_package:/sle_analysis/data:ro \
  -v $(pwd)/output:/sle_analysis/output:rw \
  --name sle-analysis \
  sle-gwas-pipeline:latest
```

### 4. **Run Analysis Scripts Inside Container**
```bash
# Inside container
python /sle_analysis/scripts/generate_therapeutic_heatmap_v2.py
```

---

## Detailed Setup

### Directory Structure
```
project-root/
├── Dockerfile                    # Container definition
├── docker-compose.yml           # Docker Compose configuration
├── requirements.txt             # Python dependencies
├── .dockerignore                # Files to exclude from build
├── DOCKER_SETUP_GUIDE.md        # This file
├── scripts/                     # Analysis scripts
│   ├── generate_missing_figures.py
│   ├── generate_therapeutic_heatmap.py
│   └── generate_therapeutic_heatmap_v2.py
├── publication_package/        # Input data (TSV files, figures)
│   ├── Tables/
│   ├── Supplementary/
│   └── Figures/
└── output/                      # Generated output (mapped to /sle_analysis/output)
```

### Build Image with Tags
```bash
# Build with specific tag
docker build -t sle-gwas-pipeline:1.0 .

# Build with multiple tags
docker build -t sle-gwas-pipeline:latest -t sle-gwas-pipeline:1.0 .

# Build with metadata
docker build \
  --label "maintainer=your-email@example.com" \
  --label "version=1.0" \
  -t sle-gwas-pipeline:latest .
```

### Check Image
```bash
docker images | grep sle-gwas
```

---

## Running Analysis Scripts

### Inside Docker Container

**Method 1: Interactive Shell**
```bash
docker-compose up -d
docker-compose exec sle-analysis bash
cd /sle_analysis
python scripts/generate_therapeutic_heatmap_v2.py
```

**Method 2: Direct Command Execution**
```bash
docker run --rm \
  -v $(pwd)/publication_package:/sle_analysis/data:ro \
  -v $(pwd)/output:/sle_analysis/output:rw \
  sle-gwas-pipeline:latest \
  python /sle_analysis/scripts/generate_therapeutic_heatmap_v2.py
```

**Method 3: Using docker-compose**
```bash
docker-compose exec sle-analysis python scripts/generate_therapeutic_heatmap_v2.py
```

### Output Files
All generated files appear in the `output/` directory on your host machine:
```
output/
├── therapeutic_mapping_optimized.png
├── therapeutic_mapping_optimized.pdf
├── Fig_3_replication_scatter.png
├── Fig_5_CLIC1_colocalization.png
├── Fig_8_PPI_network.png
└── therapeutic_targets_detailed.tsv
```

---

## Volume Mounting

### Read-Only Input Data
```bash
# Mount publication_package as read-only
docker run -it \
  -v $(pwd)/publication_package:/sle_analysis/data:ro \
  sle-gwas-pipeline:latest
```

### Write-Enabled Output
```bash
# Mount output directory for writing results
docker run -it \
  -v $(pwd)/output:/sle_analysis/output:rw \
  sle-gwas-pipeline:latest
```

### Mount Scripts for Development
```bash
# Development: mount scripts directory for live editing
docker run -it \
  -v $(pwd)/scripts:/sle_analysis/scripts:rw \
  sle-gwas-pipeline:latest
```

### Combined Mount Setup
```bash
docker run -it \
  -v $(pwd)/publication_package:/sle_analysis/data:ro \
  -v $(pwd)/output:/sle_analysis/output:rw \
  -v $(pwd)/scripts:/sle_analysis/scripts:rw \
  --name sle-analysis \
  sle-gwas-pipeline:latest
```

---

## Troubleshooting

### Issue: "Cannot find Docker"
**Solution:** Install Docker from https://docs.docker.com/get-docker/

### Issue: Permission Denied Error
**Solution:** Add your user to docker group
```bash
sudo usermod -aG docker $USER
newgrp docker
```

### Issue: Container Exits Immediately
**Solution:** Check logs
```bash
docker-compose logs sle-analysis
# or
docker logs sle-analysis
```

### Issue: Output Files Not Appearing
**Solution:** Verify volume mount
```bash
# Check mounted volumes
docker-compose exec sle-analysis ls -la /sle_analysis/output

# Verify host directory permissions
ls -la $(pwd)/output
```

### Issue: Out of Memory
**Solution:** Increase Docker memory allocation
```bash
# In docker-compose.yml, increase memory limit:
# memory: 16G (max for your system)

# Or run with memory limit
docker run -it \
  -m 8g \
  sle-gwas-pipeline:latest
```

### Issue: Slow Performance on macOS/Windows
**Solution:** Use delegated mounts
```yaml
# In docker-compose.yml
volumes:
  - ./publication_package:/sle_analysis/data:ro,delegated
  - ./output:/sle_analysis/output:rw,delegated
```

---

## Advanced Configuration

### Multi-Stage Build (Optional)
Create `Dockerfile.multistage` for smaller image:
```dockerfile
FROM python:3.11-slim as base
RUN pip install --upgrade pip

FROM base as builder
COPY requirements.txt .
RUN pip wheel --no-cache-dir --no-deps --wheel-dir /wheels -r requirements.txt

FROM base
COPY --from=builder /wheels /wheels
COPY --from=builder requirements.txt .
RUN pip install --no-cache /wheels/*
COPY . /sle_analysis
WORKDIR /sle_analysis
CMD ["bash"]
```

Build with:
```bash
docker build -f Dockerfile.multistage -t sle-gwas-pipeline:slim .
```

### Environment Variables
Create `.env` file:
```bash
DOCKER_IMAGE=sle-gwas-pipeline:latest
DOCKER_CONTAINER=sle-analysis
PYTHON_VERSION=3.11
PYTHONUNBUFFERED=1
```

### Network Configuration
```yaml
# docker-compose.yml
networks:
  sle-network:
    driver: bridge
    ipam:
      config:
        - subnet: 172.20.0.0/16
```

---

## Docker Commands Reference

### Image Management
```bash
# List all images
docker images

# Remove image
docker rmi sle-gwas-pipeline:latest

# Tag image
docker tag sle-gwas-pipeline:latest sle-gwas-pipeline:v1.0

# Push to registry (if configured)
docker push your-registry/sle-gwas-pipeline:latest
```

### Container Management
```bash
# List running containers
docker ps

# List all containers
docker ps -a

# Stop container
docker stop sle-analysis

# Remove container
docker rm sle-analysis

# View logs
docker logs sle-analysis

# Execute command
docker exec sle-analysis python --version

# Access shell
docker exec -it sle-analysis bash
```

### Docker Compose
```bash
# Build services
docker-compose build

# Start services
docker-compose up

# Start in background
docker-compose up -d

# Stop services
docker-compose down

# View logs
docker-compose logs -f sle-analysis

# Execute command
docker-compose exec sle-analysis python --version
```

---

## Publishing to Docker Registry

### Local Registry
```bash
# Tag image
docker tag sle-gwas-pipeline:latest localhost:5000/sle-gwas-pipeline:latest

# Push to local registry
docker push localhost:5000/sle-gwas-pipeline:latest
```

### Docker Hub
```bash
# Login to Docker Hub
docker login

# Tag image
docker tag sle-gwas-pipeline:latest username/sle-gwas-pipeline:latest

# Push to Docker Hub
docker push username/sle-gwas-pipeline:latest
```

### Using Published Image
```bash
docker run -it \
  -v $(pwd)/publication_package:/sle_analysis/data:ro \
  -v $(pwd)/output:/sle_analysis/output:rw \
  username/sle-gwas-pipeline:latest
```

---

## Performance Optimization

### CPU Limits
```yaml
# docker-compose.yml
deploy:
  resources:
    limits:
      cpus: '4'
    reservations:
      cpus: '2'
```

### Memory Optimization
```bash
# Run with memory limit
docker run -it \
  -m 8g \
  --memory-swap 8g \
  sle-gwas-pipeline:latest
```

### Parallel Processing
```bash
# Use multiple CPU cores
docker run -it \
  --cpus="4.0" \
  sle-gwas-pipeline:latest
```

---

## Testing the Container

### Verify Installation
```bash
docker-compose exec sle-analysis python -c "import pandas, numpy, matplotlib; print('All packages installed!')"
```

### Run Test Script
```bash
docker-compose exec sle-analysis bash -c "
python scripts/generate_therapeutic_heatmap_v2.py
ls -la output/
"
```

### Check Data Access
```bash
docker-compose exec sle-analysis bash -c "
ls -la /sle_analysis/data/
echo 'Data accessible!'
"
```

---

## Cleanup

### Remove Stopped Containers
```bash
docker container prune

# or
docker-compose down
```

### Remove Unused Images
```bash
docker image prune
```

### Remove All Docker Resources
```bash
# WARNING: Remove all containers, images, networks
docker system prune -a
```

---

## Summary

✅ **Dockerfile**: Defines reproducible container environment  
✅ **docker-compose.yml**: Simplifies container orchestration  
✅ **requirements.txt**: Manages Python dependencies  
✅ **Volume Mounting**: Handles data I/O  
✅ **Network Configuration**: Isolated environment  

The SLE GWAS pipeline is now fully containerized and can run on any system with Docker installed!

---

## Support

For issues or questions:
1. Check Docker logs: `docker-compose logs sle-analysis`
2. Verify volume mounts: `docker-compose exec sle-analysis ls -la /sle_analysis/`
3. Test individual scripts: `docker-compose exec sle-analysis python scripts/...`
4. Review Dockerfile for base image details

---

**Generated:** April 19, 2026  
**Version:** 1.0  
**Status:** Production Ready ✅
