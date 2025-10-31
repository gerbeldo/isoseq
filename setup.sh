#!/usr/bin/env bash
set -euxo pipefail

# Detect the real target user & home even if invoked via sudo
if [ -n "${SUDO_USER:-}" ] && [ "$SUDO_USER" != "root" ]; then
  TARGET_USER="$SUDO_USER"
else
  TARGET_USER="$(id -un)"
fi
TARGET_HOME="$(getent passwd "$TARGET_USER" | cut -d: -f6)"

# Convenience helpers to run as the target user
run_as_user() { sudo -u "$TARGET_USER" --preserve-env=PATH "$@"; }

# Optional: minimal deps (only if we have apt privileges)
if command -v apt-get >/dev/null 2>&1; then
  sudo apt-get update
  sudo apt-get install -y --no-install-recommends curl ca-certificates bzip2 xz-utils tar
  sudo apt-get clean
  sudo rm -rf /var/lib/apt/lists/*
fi

# Paths inside the user's home
BIN_DIR="$TARGET_HOME/.local/bin"
DL_DIR="$TARGET_HOME/.cache/mamba-dl"
export MAMBA_ROOT_PREFIX="$TARGET_HOME/mamba"

# Create dirs with correct ownership/permissions
sudo -u "$TARGET_USER" mkdir -p "$BIN_DIR" "$DL_DIR" "$MAMBA_ROOT_PREFIX"
sudo chown -R "$TARGET_USER":"$TARGET_USER" "$TARGET_HOME/.local" "$TARGET_HOME/.cache" "$MAMBA_ROOT_PREFIX"

# Pick micromamba URL by arch (Linux x86_64 vs aarch64)
ARCH="$(uname -m)"
case "$ARCH" in
  x86_64)  MAMBA_URL="https://micro.mamba.pm/api/micromamba/linux-64/latest" ;;
  aarch64) MAMBA_URL="https://micro.mamba.pm/api/micromamba/linux-aarch64/latest" ;;
  *) echo "Unsupported arch: $ARCH" >&2; exit 1 ;;
esac

# Download to the user's cache dir and install into ~/.local/bin
run_as_user curl -fL -o "$DL_DIR/micromamba.tar.bz2" "$MAMBA_URL"
tar -xvjf "$DL_DIR/micromamba.tar.bz2" -C "$BIN_DIR" --strip-components=1 bin/micromamba
sudo chown "$TARGET_USER":"$TARGET_USER" "$BIN_DIR/micromamba"
chmod 0755 "$BIN_DIR/micromamba"
rm -f "$DL_DIR/micromamba.tar.bz2"

# Ensure micromamba & env root are on PATH for the target user
BASHRC="$TARGET_HOME/.bashrc"
PROFILE="$TARGET_HOME/.profile"
append_export() {
  local file="$1"
  sudo -u "$TARGET_USER" bash -c "grep -q 'MAMBA_ROOT_PREFIX=' '$file' || echo 'export MAMBA_ROOT_PREFIX=$MAMBA_ROOT_PREFIX' >> '$file'"
  sudo -u "$TARGET_USER" bash -c "grep -q '\$HOME/.local/bin' '$file' || echo 'export PATH=\"\$HOME/.local/bin:\$PATH\"' >> '$file'"
  sudo -u "$TARGET_USER" bash -c "grep -q '\$MAMBA_ROOT_PREFIX/bin' '$file' || echo 'export PATH=\"\$MAMBA_ROOT_PREFIX/bin:\$PATH\"' >> '$file'"
}
append_export "$BASHRC"
append_export "$PROFILE"

# Create the base env with your tools (run as the target user)
run_as_user env MAMBA_ROOT_PREFIX="$MAMBA_ROOT_PREFIX" "$BIN_DIR/micromamba" create -y -n base -c conda-forge -c bioconda \
  isoseq3 lima pbmm2 pbpigeon samtools
run_as_user env MAMBA_ROOT_PREFIX="$MAMBA_ROOT_PREFIX" "$BIN_DIR/micromamba" clean -a -y

# Sanity check (new shells will source bashrc/profile; source here for the current session)
run_as_user bash -lc 'source "$HOME/.profile" 2>/dev/null || true; source "$HOME/.bashrc" 2>/dev/null || true; \
  command -v isoseq3 && isoseq3 --version && \
  command -v lima && lima --version && \
  command -v pbmm2 && pbmm2 --version && \
  command -v samtools && samtools --version && \
  pbpigeon --help >/dev/null'

echo "✅ Install complete. Open a new shell or run: source ~/.profile && source ~/.bashrc"
echo "Env root: $MAMBA_ROOT_PREFIX"
echo "Micromamba: $BIN_DIR/micromamba"
