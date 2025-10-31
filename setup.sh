#!/usr/bin/env bash
# install_pb_tools.sh
# Ubuntu (x86_64) micromamba + PacBio stack with global wrappers.

set -euo pipefail

MAMBA_ROOT_PREFIX="/opt/micromamba"
ENV_NAME="pb"
WRAP_DIR="/usr/local/bin"

# --- Assert amd64 Ubuntu ---
if [[ "$(dpkg --print-architecture)" != "amd64" ]]; then
  echo "[x] This script targets x86_64 (amd64). Detected: $(dpkg --print-architecture)"
  exit 1
fi

# --- Minimal deps ---
export DEBIAN_FRONTEND=noninteractive
sudo apt-get update -y
sudo apt-get install -y --no-install-recommends curl ca-certificates bzip2 tar
sudo update-ca-certificates

# --- Install micromamba binary on PATH ---
echo "[+] Installing micromamba to /usr/local/bin ..."
TMP_TAR="$(mktemp)"
curl -fsSL "https://micro.mamba.pm/api/micromamba/linux-64/latest" -o "${TMP_TAR}"
sudo tar -xj -f "${TMP_TAR}" -C /usr/local/bin --strip-components=1 bin/micromamba
rm -f "${TMP_TAR}"
/usr/local/bin/micromamba --version

# --- Prefix + profile hook (for interactive use) ---
echo "[+] Preparing prefix at ${MAMBA_ROOT_PREFIX} ..."
sudo mkdir -p "${MAMBA_ROOT_PREFIX}"
MAMBA_USER="${SUDO_USER:-$USER}"
sudo chown -R "${MAMBA_USER}:${MAMBA_USER}" "${MAMBA_ROOT_PREFIX}"

echo "[+] Writing /etc/profile.d/micromamba.sh ..."
sudo tee /etc/profile.d/micromamba.sh >/dev/null <<'EOF'
export MAMBA_ROOT_PREFIX=/opt/micromamba
if command -v micromamba >/dev/null 2>&1; then
  eval "$(micromamba shell hook -s bash)"
fi
EOF
sudo chmod 0644 /etc/profile.d/micromamba.sh

# --- Create env + install tools ---
echo "[+] Creating env '${ENV_NAME}' and installing tools ..."
sudo -u "${MAMBA_USER}" bash -lc '
  set -euo pipefail
  export MAMBA_ROOT_PREFIX='"${MAMBA_ROOT_PREFIX}"'
  eval "$(micromamba shell hook -s bash)"

  micromamba create -y -n '"${ENV_NAME}"' \
    -c conda-forge -c bioconda --channel-priority strict \
    isoseq3 lima pbmm2 pbpigeon samtools

  micromamba clean -a -y
'

# --- Create global wrapper scripts so no activation is needed ---
echo "[+] Creating global wrappers in ${WRAP_DIR} ..."
ENV_PREFIX="${MAMBA_ROOT_PREFIX}/envs/${ENV_NAME}"

make_wrapper () {
  local name="$1" real_cmd="$2"
  sudo tee "${WRAP_DIR}/${name}" >/dev/null <<EOF
#!/usr/bin/env bash
# Auto-run ${real_cmd} from micromamba env ${ENV_NAME}
exec micromamba run -p "${ENV_PREFIX}" ${real_cmd} "\$@"
EOF
  sudo chmod 0755 "${WRAP_DIR}/${name}"
}

# Detect which iso command exists inside the env and create sensible aliases.
HAVE_ISOSEQ3=$(sudo -u "${MAMBA_USER}" micromamba run -p "${ENV_PREFIX}" bash -lc 'command -v isoseq3 >/dev/null && echo yes || true')
HAVE_ISOSEQ=$(sudo -u "${MAMBA_USER}" micromamba run -p "${ENV_PREFIX}" bash -lc 'command -v isoseq >/dev/null && echo yes || true')

if [[ -n "${HAVE_ISOSEQ3}" ]]; then
  make_wrapper "isoseq3" "isoseq3"
  # also provide 'isoseq' alias if not present
  if [[ -z "${HAVE_ISOSEQ}" ]]; then
    make_wrapper "isoseq" "isoseq3"
  fi
fi

if [[ -n "${HAVE_ISOSEQ}" ]]; then
  make_wrapper "isoseq" "isoseq"
  # also provide 'isoseq3' alias if not present
  if [[ -z "${HAVE_ISOSEQ3}" ]]; then
    make_wrapper "isoseq3" "isoseq"
  fi
fi

# PacBio + samtools wrappers
for cmd in lima pbmm2 pigeon samtools; do
  make_wrapper "${cmd}" "${cmd}"
done

# --- Smoke check (non-fatal) ---
echo "[+] Verifying versions (non-fatal) ..."
set +e
isoseq3 --version 2>/dev/null || isoseq --version 2>/dev/null
lima --version 2>/dev/null
pbmm2 --version 2>/dev/null
pigeon --help 2>/dev/null | head -n1
samtools --version 2>/dev/null | head -n1
set -e

cat <<MSG

[✓] Done.

You can now run tools directly (no activation needed):
  isoseq3 --help
  lima --version
  pbmm2 --version
  pigeon --help
  samtools --version

(Interactive shells can still:  source /etc/profile && micromamba activate ${ENV_NAME})
MSG
