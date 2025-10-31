#!/usr/bin/env bash
set -euxo pipefail

# This script targets a single login user. It is safe to re-run.
# It installs:
# - micromamba → ~/.local/bin/micromamba
# - root prefix (base env) → ~/mamba
# - auto shell init + auto-activate for bash/zsh

#### 0) Resolve the real target user & home, even if invoked with sudo
if [ -n "${SUDO_USER:-}" ] && [ "$SUDO_USER" != "root" ]; then
  TARGET_USER="$SUDO_USER"
else
  TARGET_USER="$(id -un)"
fi
TARGET_HOME="$(getent passwd "$TARGET_USER" | cut -d: -f6)"
BIN_DIR="$TARGET_HOME/.local/bin"
DL_DIR="$TARGET_HOME/.cache/mamba-dl"
export MAMBA_ROOT_PREFIX="$TARGET_HOME/mamba"   # root prefix == "base" env location

# Helper to run as the target user with their env
run_user() { sudo -u "$TARGET_USER" --preserve-env=PATH,HOME bash -lc "$*"; }

#### 1) Minimal deps (only uses sudo for apt)
if command -v apt-get >/dev/null 2>&1; then
  sudo apt-get update
  sudo apt-get install -y --no-install-recommends curl ca-certificates bzip2 xz-utils tar
  sudo apt-get clean
  sudo rm -rf /var/lib/apt/lists/*
fi

#### 2) Prepare user directories with correct ownership
sudo -u "$TARGET_USER" mkdir -p "$BIN_DIR" "$DL_DIR" "$MAMBA_ROOT_PREFIX"
sudo chown -R "$TARGET_USER":"$TARGET_USER" "$TARGET_HOME/.local" "$TARGET_HOME/.cache" "$MAMBA_ROOT_PREFIX"

#### 3) Install micromamba (no /tmp; correct archive; arch-aware)
ARCH="$(uname -m)"
case "$ARCH" in
  x86_64)  MAMBA_URL="https://micro.mamba.pm/api/micromamba/linux-64/latest" ;;
  aarch64) MAMBA_URL="https://micro.mamba.pm/api/micromamba/linux-aarch64/latest" ;;
  *) echo "Unsupported arch: $ARCH" >&2; exit 1 ;;
esac

run_user "curl -fL -o '$DL_DIR/micromamba.tar.bz2' '$MAMBA_URL'"
# Extract the single binary to ~/.local/bin
tar -xvjf "$DL_DIR/micromamba.tar.bz2" -C "$BIN_DIR" --strip-components=1 bin/micromamba
sudo chown "$TARGET_USER":"$TARGET_USER" "$BIN_DIR/micromamba"
chmod 0755 "$BIN_DIR/micromamba"
rm -f "$DL_DIR/micromamba.tar.bz2"

# Ensure ~/.local/bin is first on PATH for this user (bash + login shells)
for FILE in "$TARGET_HOME/.bashrc" "$TARGET_HOME/.profile"; do
  sudo -u "$TARGET_USER" touch "$FILE"
  if ! sudo -u "$TARGET_USER" grep -q '\$HOME/.local/bin' "$FILE"; then
    echo 'export PATH="$HOME/.local/bin:$PATH"' | sudo tee -a "$FILE" >/dev/null
  fi
done

#### 4) Create/upgrade the "base" env at the root prefix
# IMPORTANT: use -p $MAMBA_ROOT_PREFIX (not -n base) to target the root prefix
if [ -d "$MAMBA_ROOT_PREFIX/conda-meta" ]; then
  # already exists → ensure packages are present/updated
  run_user "MAMBA_ROOT_PREFIX='$MAMBA_ROOT_PREFIX' '$BIN_DIR/micromamba' install -y -p '$MAMBA_ROOT_PREFIX' -c conda-forge -c bioconda \
            isoseq3 lima pbmm2 pbpigeon samtools"
else
  run_user "MAMBA_ROOT_PREFIX='$MAMBA_ROOT_PREFIX' '$BIN_DIR/micromamba' create  -y -p '$MAMBA_ROOT_PREFIX' -c conda-forge -c bioconda \
            isoseq3 lima pbmm2 pbpigeon samtools"
fi
run_user "MAMBA_ROOT_PREFIX='$MAMBA_ROOT_PREFIX' '$BIN_DIR/micromamba' clean -a -y -p '$MAMBA_ROOT_PREFIX'"

#### 5) Shell init + auto-activation (bash & zsh if present)
# Add the micromamba hook and auto-activate base for interactive shells.
# (Keeps things simple for interactive use; comment out the 'micromamba activate' line if you prefer manual activation.)

# Bash
BASHRC="$TARGET_HOME/.bashrc"
if ! sudo -u "$TARGET_USER" grep -q 'micromamba shell hook -s bash' "$BASHRC"; then
  run_user "eval \"\$('$BIN_DIR/micromamba' shell hook -s bash)\" >/dev/null 2>&1 || true"
  {
    echo ''
    echo '# micromamba init (bash)'
    echo 'eval "$($HOME/.local/bin/micromamba shell hook -s bash)"'
    echo "export MAMBA_ROOT_PREFIX=\"$MAMBA_ROOT_PREFIX\""
    echo 'micromamba activate'   # auto-activate base
  } | sudo tee -a "$BASHRC" >/dev/null
fi

# Zsh (only if user has zsh config already)
if [ -f "$TARGET_HOME/.zshrc" ]; then
  ZSHRC="$TARGET_HOME/.zshrc"
  if ! sudo -u "$TARGET_USER" grep -q 'micromamba shell hook -s zsh' "$ZSHRC"; then
    {
      echo ''
      echo '# micromamba init (zsh)'
      echo 'eval "$($HOME/.local/bin/micromamba shell hook -s zsh)"'
      echo "export MAMBA_ROOT_PREFIX=\"$MAMBA_ROOT_PREFIX\""
      echo 'micromamba activate'
    } | sudo tee -a "$ZSHRC" >/dev/null
  fi
fi

#### 6) Sanity checks (new interactive shell to verify auto-activation + commands)
run_user 'bash -lc "isoseq3 --version && lima --version && pbmm2 --version && samtools --version && pbpigeon --help >/dev/null"'

echo "✅ Done. Open a NEW shell (or 'exec \$SHELL') and your base env will be active automatically."
echo "   micromamba: $BIN_DIR/micromamba"
echo "   root prefix (base env): $MAMBA_ROOT_PREFIX"
