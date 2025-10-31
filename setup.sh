#!/bin/bash
set -euo pipefail

# Bioinformatics Tools Setup for Iso-Seq Analysis
# This script installs micromamba and required bioinformatics tools on Ubuntu EC2
# Tools are configured to be available in PATH automatically at startup

echo "Starting bioinformatics environment setup..."

# Update system packages
echo "Updating system packages..."
apt-get update
apt-get upgrade -y

# Install system dependencies required by bioinformatics tools
echo "Installing system dependencies..."
apt-get install -y \
    wget \
    curl \
    git \
    build-essential \
    ca-certificates \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libffi-dev

# Create a dedicated directory for micromamba installation
MICROMAMBA_DIR="/opt/micromamba"
MICROMAMBA_ENV_DIR="/opt/micromamba/envs/isoseq"
MICROMAMBA_BIN="${MICROMAMBA_DIR}/bin/micromamba"

echo "Installing micromamba to ${MICROMAMBA_DIR}..."

# Create directory structure
mkdir -p "${MICROMAMBA_DIR}"
mkdir -p "${MICROMAMBA_ENV_DIR}"

# Download and install micromamba
MICROMAMBA_VERSION=$(curl -s https://api.github.com/repos/mamba-org/micromamba-releases/releases/latest | grep tag_name | cut -d '"' -f 4)
MICROMAMBA_DOWNLOAD_URL="https://github.com/mamba-org/micromamba-releases/releases/download/${MICROMAMBA_VERSION}/micromamba-linux-64"

echo "Downloading micromamba version ${MICROMAMBA_VERSION}..."
curl -fsSL "${MICROMAMBA_DOWNLOAD_URL}" -o "${MICROMAMBA_BIN}"
chmod +x "${MICROMAMBA_BIN}"

# Create isoseq conda environment with required tools
echo "Creating isoseq conda environment with bioinformatics tools..."
export MAMBA_ROOT_PREFIX="${MICROMAMBA_DIR}"
"${MICROMAMBA_BIN}" create \
    -y \
    -p "${MICROMAMBA_ENV_DIR}" \
    -c bioconda \
    -c conda-forge \
    isoseq3 \
    lima \
    pbmm2 \
    pbpigeon \
    samtools

# Create activation wrapper script to make tools available in PATH by default
echo "Setting up automatic activation on shell startup..."

PROFILE_SCRIPT="/etc/profile.d/isoseq-env.sh"
cat > "${PROFILE_SCRIPT}" << EOF
# Automatically activate isoseq environment
export MAMBA_EXE="${MICROMAMBA_BIN}"
export MAMBA_ROOT_PREFIX="${MICROMAMBA_DIR}"
eval "\$("${MICROMAMBA_BIN}" shell hook --shell bash --root-prefix "${MICROMAMBA_DIR}" 2>/dev/null)" || true
micromamba activate "${MICROMAMBA_ENV_DIR}" 2>/dev/null || true
EOF

chmod 644 "${PROFILE_SCRIPT}"

# Create bashrc configuration for non-login shells and interactive sessions
echo "Configuring bashrc for interactive shell sessions..."

BASHRC_CONFIG="/opt/micromamba/bashrc-config.sh"
cat > "${BASHRC_CONFIG}" << EOF
# Micromamba and isoseq environment configuration
export MAMBA_EXE="${MICROMAMBA_BIN}"
export MAMBA_ROOT_PREFIX="${MICROMAMBA_DIR}"
eval "\$("${MICROMAMBA_BIN}" shell hook --shell bash --root-prefix "${MICROMAMBA_DIR}" 2>/dev/null)" || true
micromamba activate "${MICROMAMBA_ENV_DIR}" 2>/dev/null || true
EOF

chmod 644 "${BASHRC_CONFIG}"

# Update root's bashrc to source the isoseq environment
if ! grep -q "isoseq-env" /root/.bashrc; then
    echo "source ${BASHRC_CONFIG}" >> /root/.bashrc
fi

# Create a verification script to test installations
echo "Creating verification script..."

VERIFY_SCRIPT="/usr/local/bin/verify-isoseq-tools"
cat > "${VERIFY_SCRIPT}" << 'EOF'
#!/bin/bash

echo "Verifying bioinformatics tools installation..."
echo ""

tools=("isoseq3" "lima" "pbmm2" "pbpigeon" "samtools")
all_ok=true

for tool in "${tools[@]}"; do
    if command -v "${tool}" &> /dev/null; then
        version=$("${tool}" --version 2>&1 | head -n 1 || echo "version check not available")
        echo "âœ“ ${tool}: ${version}"
    else
        echo "âœ— ${tool}: NOT FOUND"
        all_ok=false
    fi
done

echo ""
if [ "$all_ok" = true ]; then
    echo "All tools are correctly installed!"
    exit 0
else
    echo "Some tools are missing. Please check the installation."
    exit 1
fi
EOF

chmod +x "${VERIFY_SCRIPT}"

# Final verification
echo ""
echo "Running final verification..."
export MAMBA_EXE="${MICROMAMBA_BIN}"
export MAMBA_ROOT_PREFIX="${MICROMAMBA_DIR}"
eval "$("${MICROMAMBA_BIN}" shell hook --shell bash --root-prefix "${MICROMAMBA_DIR}" 2>/dev/null)" || true
micromamba activate "${MICROMAMBA_ENV_DIR}" 2>/dev/null || true

"${VERIFY_SCRIPT}"

echo ""
echo "Setup complete!"
echo ""
echo "To use the isoseq environment in new shell sessions, tools will be automatically available."
echo "If needed, you can manually activate with:"
echo "  source /opt/micromamba/bashrc-config.sh"
echo ""
echo "To verify installation, run: verify-isoseq-tools"
