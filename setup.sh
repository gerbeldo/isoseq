#!/bin/bash

# Bioinformatics EC2 Setup Script
# This script installs micromamba and bioconda packages for PacBio analysis
# Designed for Ubuntu EC2 instances

set -e  # Exit on any error

echo "=========================================="
echo "Starting bioinformatics environment setup"
echo "=========================================="

# Update system packages
echo "Updating system packages..."
sudo apt-get update
sudo apt-get upgrade -y

# Install basic dependencies
echo "Installing system dependencies..."
sudo apt-get install -y \
    curl \
    bzip2 \
    ca-certificates \
    git \
    wget

# Install micromamba
echo "Installing micromamba..."
# The download URL returns a tar.bz2 archive
# We pipe it directly to tar which extracts the bin/micromamba file
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba

# Move to a permanent location
sudo mkdir -p /opt/micromamba
sudo mv bin/micromamba /opt/micromamba/
sudo chmod +x /opt/micromamba/micromamba
rmdir bin  # Clean up the temporary directory

# Initialize micromamba for current user
echo "Initializing micromamba..."
/opt/micromamba/micromamba shell init -s bash --root-prefix ~/micromamba

# Source the shell configuration to make micromamba available
export MAMBA_ROOT_PREFIX=~/micromamba
eval "$(/opt/micromamba/micromamba shell hook -s bash)"

# Configure conda channels
echo "Configuring bioconda channels..."
micromamba config append channels conda-forge
micromamba config append channels bioconda
micromamba config append channels defaults
micromamba config set channel_priority strict

# Create a base environment with bioinformatics tools
echo "Creating bioinformatics environment..."
micromamba create -n bioinfo -y python=3.11

# Activate the environment
micromamba activate bioinfo

# Install bioinformatics packages
echo "Installing bioinformatics tools from bioconda..."
micromamba install -y \
    isoseq3 \
    lima \
    pbmm2 \
    pbpigeon \
    samtools

# Verify installations
echo ""
echo "=========================================="
echo "Verifying installations..."
echo "=========================================="

micromamba activate bioinfo

echo "isoseq3 version:"
isoseq3 --version 2>&1 | head -n 1 || echo "isoseq3 check failed"

echo "lima version:"
lima --version 2>&1 | head -n 1 || echo "lima check failed"

echo "pbmm2 version:"
pbmm2 --version 2>&1 | head -n 1 || echo "pbmm2 check failed"

echo "pbpigeon version:"
pbpigeon --version 2>&1 | head -n 1 || echo "pbpigeon check failed"

echo "samtools version:"
samtools --version 2>&1 | head -n 1 || echo "samtools check failed"

echo ""
echo "=========================================="
echo "Setup complete!"
echo "=========================================="
echo ""
echo "To use the bioinformatics tools:"
echo "1. Start a new shell session or run: source ~/.bashrc"
echo "2. Activate the environment: micromamba activate bioinfo"
echo "3. Run your tools: isoseq3, lima, pbmm2, pbpigeon, samtools"
echo ""
echo "The environment will be automatically available after creating an AMI."
echo "=========================================="
