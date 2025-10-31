#!/bin/bash

# Bioinformatics EC2 Setup Script
# This script installs micromamba and bioconda packages globally
# Tools will be available system-wide without activating environments
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
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba

# Move to a permanent location
sudo mkdir -p /opt/micromamba/bin
sudo mv bin/micromamba /opt/micromamba/bin/
sudo chmod +x /opt/micromamba/bin/micromamba
rmdir bin  # Clean up the temporary directory

# Set up micromamba with a global installation directory
echo "Setting up micromamba..."
export MAMBA_ROOT_PREFIX=/opt/micromamba/envs

# Create the base environment directory
sudo mkdir -p $MAMBA_ROOT_PREFIX

# Configure conda channels (these commands don't require activation)
echo "Configuring bioconda channels..."
sudo /opt/micromamba/bin/micromamba config append channels conda-forge --root-prefix $MAMBA_ROOT_PREFIX
sudo /opt/micromamba/bin/micromamba config append channels bioconda --root-prefix $MAMBA_ROOT_PREFIX
sudo /opt/micromamba/bin/micromamba config append channels defaults --root-prefix $MAMBA_ROOT_PREFIX
sudo /opt/micromamba/bin/micromamba config set channel_priority strict --root-prefix $MAMBA_ROOT_PREFIX

# Create a base environment with bioinformatics tools
echo "Creating base environment with bioinformatics tools..."
echo "This may take several minutes..."
sudo /opt/micromamba/bin/micromamba create -n base -y \
    python=3.11 \
    isoseq3 \
    lima \
    pbmm2 \
    pbpigeon \
    samtools \
    --root-prefix $MAMBA_ROOT_PREFIX \
    -c conda-forge \
    -c bioconda \
    -c defaults

# Add tools to system PATH by creating symlinks
echo "Making tools available system-wide..."
sudo mkdir -p /usr/local/bin

# Create symlinks for all the tools
for tool in isoseq3 lima pbmm2 pbpigeon samtools; do
    if [ -f "$MAMBA_ROOT_PREFIX/base/bin/$tool" ]; then
        sudo ln -sf "$MAMBA_ROOT_PREFIX/base/bin/$tool" /usr/local/bin/$tool
        echo "Linked $tool to /usr/local/bin"
    else
        echo "Warning: $tool not found at $MAMBA_ROOT_PREFIX/base/bin/$tool"
    fi
done

# Also symlink commonly used dependencies
echo "Linking additional dependencies..."
for dep in python python3 pip; do
    if [ -f "$MAMBA_ROOT_PREFIX/base/bin/$dep" ]; then
        sudo ln -sf "$MAMBA_ROOT_PREFIX/base/bin/$dep" /usr/local/bin/$dep
    fi
done

# Verify installations
echo ""
echo "=========================================="
echo "Verifying installations..."
echo "=========================================="

export PATH="/usr/local/bin:$PATH"

echo "isoseq3 version:"
/usr/local/bin/isoseq3 --version 2>&1 | head -n 1 || echo "isoseq3 check failed"

echo ""
echo "lima version:"
/usr/local/bin/lima --version 2>&1 | head -n 1 || echo "lima check failed"

echo ""
echo "pbmm2 version:"
/usr/local/bin/pbmm2 --version 2>&1 | head -n 1 || echo "pbmm2 check failed"

echo ""
echo "pbpigeon version:"
/usr/local/bin/pbpigeon --version 2>&1 | head -n 1 || echo "pbpigeon check failed"

echo ""
echo "samtools version:"
/usr/local/bin/samtools --version 2>&1 | head -n 1 || echo "samtools check failed"

echo ""
echo "=========================================="
echo "Setup complete!"
echo "=========================================="
echo ""
echo "All bioinformatics tools are now available system-wide."
echo "You can run them directly from any terminal:"
echo "  - isoseq3"
echo "  - lima"
echo "  - pbmm2"
echo "  - pbpigeon"
echo "  - samtools"
echo ""
echo "No environment activation needed!"
echo "The tools will be automatically available after creating an AMI."
echo "=========================================="
