#!/usr/bin/env bash
set -euxo pipefail

# This script is intended to run as root (on an Ubuntu 24.04 builder EC2).
# If not root, re-run with:  sudo -H bash setup.sh

export DEBIAN_FRONTEND=noninteractive

### 0) Minimal deps
apt-get update
apt-get install -y --no-install-recommends curl ca-certificates bzip2 xz-utils tar
apt-get clean
rm -rf /var/lib/apt/lists/*

### 1) Micromamba (global install) + PATH for all users
MICRO=/usr/local/bin/micromamba

# Pick the right micromamba build based on arch
ARCH="$(uname -m)"
case "$ARCH" in
  x86_64)   MAMBA_URL="https://micro.mamba.pm/api/micromamba/linux-64/latest" ;;
  aarch64)  MAMBA_URL="https://micro.mamba.pm/api/micromamba/linux-aarch64/latest" ;;
  *) echo "Unsupported arch: $ARCH" >&2; exit 1 ;;
esac

# Download (bz2) and install
curl -L -o /tmp/micromamba.tar.bz2 "$MAMBA_URL"
tar -xvjf /tmp/micromamba.tar.bz2 -C /usr/local/bin --strip-components=1 bin/micromamba
chmod 0755 "$MICRO"
rm -f /tmp/micromamba.tar.bz2

# Global env root for micromamba
export MAMBA_ROOT_PREFIX=/opt/conda
mkdir -p "$MAMBA_ROOT_PREFIX"
chmod 755 "$MAMBA_ROOT_PREFIX"

# PATH for interactive shells
cat >/etc/profile.d/mamba.sh <<'EOF'
export MAMBA_ROOT_PREFIX=/opt/conda
export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
EOF
chmod 0644 /etc/profile.d/mamba.sh

# PATH for non-interactive env (no "export" here)
# (Ensure we don't duplicate entries)
grep -q '^MAMBA_ROOT_PREFIX=' /etc/environment || echo 'MAMBA_ROOT_PREFIX=/opt/conda' >> /etc/environment
if ! grep -q '^PATH=.*\/opt\/conda\/bin' /etc/environment; then
  sed -i 's#^PATH=.*#PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin#' /etc/environment || \
  echo 'PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin' >> /etc/environment
fi

### 2) Reproduce your Dockerfile env (base env with the tools)
# Matches: isoseq3 lima pbmm2 pbpigeon samtools
"$MICRO" create -y -n base -c conda-forge -c bioconda \
  isoseq3 lima pbmm2 pbpigeon samtools

# Clean caches
"$MICRO" clean -a -y

# (Optional) convenience symlinks
ln -sf /opt/conda/bin/* /usr/local/bin/ || true

### 3) Sanity checks (new login shells will pick PATH; we source it here)
source /etc/profile.d/mamba.sh
isoseq3 --version
lima --version
pbmm2 --version
pbpigeon --help >/dev/null
samtools --version

### 4) (Optional) export env for provenance
/opt/conda/bin/conda env export -n base --no-builds > /root/isoseq-base-env.yml || true

### 5) Tidy before baking AMI
apt-get clean
rm -rf /var/lib/apt/lists/*
cloud-init clean
rm -f /root/.bash_history
