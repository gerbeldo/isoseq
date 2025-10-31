#!/usr/bin/env bash
set -euxo pipefail

### 0) Basics + locale (fixes the LC_ALL warnings)
export DEBIAN_FRONTEND=noninteractive
apt-get update
apt-get install -y --no-install-recommends \
  curl ca-certificates bzip2 xz-utils tar locales

# Enable and set en_US.UTF-8
sed -i 's/^# *\(en_US.UTF-8 UTF-8\)/\1/' /etc/locale.gen
locale-gen en_US.UTF-8
update-locale LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8
echo -e 'LANG=en_US.UTF-8\nLC_ALL=en_US.UTF-8' >/etc/default/locale

### 1) Micromamba (global) + env path for all users
MICRO=/usr/local/bin/micromamba
curl -L https://micro.mamba.pm/api/micromamba/linux-64/latest \
  | tar -xJ -C /usr/local/bin --strip-components=1 bin/micromamba
chmod +x "$MICRO"

# Global install root
export MAMBA_ROOT_PREFIX=/opt/conda
mkdir -p "$MAMBA_ROOT_PREFIX"
chmod 755 "$MAMBA_ROOT_PREFIX"

# Make PATH persistent for everyone (login + non-login shells)
cat >/etc/profile.d/mamba.sh <<'EOF'
export MAMBA_ROOT_PREFIX=/opt/conda
export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
EOF
chmod 644 /etc/profile.d/mamba.sh
# Also for non-interactive shells invoked without sourcing /etc/profile
echo 'export MAMBA_ROOT_PREFIX=/opt/conda' >> /etc/environment
echo 'PATH="/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"' >> /etc/environment

### 2) Reproduce your Dockerfile env (base env with the tools)
# Channels explicit; keep it minimal like your Dockerfile
"$MICRO" create -y -n base -c conda-forge -c bioconda \
  isoseq3 lima pbmm2 pbpigeon samtools
"$MICRO" clean -a -y

# Convenience: symlink common tools into /usr/local/bin
ln -sf /opt/conda/bin/* /usr/local/bin/ || true

### 3) Sanity checks (will fail early if anything’s off)
source /etc/profile.d/mamba.sh
isoseq3 --version
lima --version
pbmm2 --version
pigeon --version
samtools --version
locale

### 4) (Optional) export env for provenance
/opt/conda/bin/conda env export -n base --no-builds > /root/isoseq-base-env.yml

### 5) Tidy up before imaging
apt-get clean
rm -rf /var/lib/apt/lists/*
cloud-init clean
rm -f /root/.bash_history
