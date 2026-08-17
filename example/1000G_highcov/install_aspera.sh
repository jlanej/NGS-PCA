#!/usr/bin/env bash
# =============================================================================
# install_aspera.sh — Provision a working ascp inside WORK_DIR
# =============================================================================
#
# Two separate problems have to be solved to get Aspera working against ENA,
# and they are solved from two different places:
#
#   1. THE BINARY. fasp.sra.ebi.ac.uk now runs OpenSSH 8.7 and offers only
#      modern SSH transport algorithms (curve25519 / ecdh / DH-GEX-SHA256 for
#      key exchange, and encrypt-then-MAC MACs only). The libssh2 bundled with
#      ascp 3.9.x supports none of them, so its handshake dies at algorithm
#      negotiation with "Failure Event: -5 - Unable to exchange encryption
#      keys", which ascp then misreports as "failed to authenticate". Fixed by
#      installing a current Aspera Connect.
#
#   2. THE KEY. Connect releases after 4.1 no longer ship the anonymous
#      ENA/EBI key 'asperaweb_id_dsa.openssh', and it is not hosted standalone
#      anywhere authoritative. It only exists inside older installs. ENA does
#      still accept ssh-dss for user authentication, so the old key remains
#      valid — it just has to be found rather than downloaded.
#
# So: the NEW installer supplies the binary, an OLD install supplies the key.
# Nothing sensitive is committed to this repository.
#
# Usage:
#   bash install_aspera.sh          # idempotent; safe to re-run
#   module load aspera && bash install_aspera.sh   # helps step 2 find the key
#
# After a successful run, config.sh picks up both automatically on subsequent
# pipeline invocations.
# =============================================================================

set -euo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

ASPERA_HOME="${ASPERA_HOME:-${WORK_DIR}/aspera}"
ASPERA_BIN="${ASPERA_HOME}/.aspera/connect/bin/ascp"
ASPERA_KEY_DEST="${ASPERA_HOME}/etc/asperaweb_id_dsa.openssh"

# Pinned so a silent upstream change cannot alter what lands on the cluster.
# Override CONNECT_VERSION/CONNECT_URL to move to a newer build.
CONNECT_VERSION="${CONNECT_VERSION:-4.2.13.820}"
CONNECT_TARBALL="ibm-aspera-connect_${CONNECT_VERSION}_linux_x86_64.tar.gz"
CONNECT_URL="${CONNECT_URL:-https://d3gcli72yxqn2z.cloudfront.net/downloads/connect/latest/bin/${CONNECT_TARBALL}}"
# sha256 of 4.2.13.820; blank disables verification (required when overriding
# CONNECT_VERSION, since the checksum would no longer match).
CONNECT_SHA256="${CONNECT_SHA256:-870885149ccaf57a794b1e8f50885cac8d454768a615e4c11a0f1356872b989e}"

echo "============================================================"
echo " Aspera provisioning"
echo "   install root : ${ASPERA_HOME}"
echo "   Connect ver  : ${CONNECT_VERSION}"
echo "============================================================"

# ── Step 1: the binary ──────────────────────────────────────────────────────

install_connect() {
  local tmp
  tmp="$(mktemp -d "${WORK_DIR}/.aspera-install.XXXXXX")"
  # shellcheck disable=SC2064  # expand tmp now, not at trap time
  trap "rm -rf '${tmp}'" RETURN

  echo "[1/2] Downloading Aspera Connect ${CONNECT_VERSION} (~68 MB)..."
  if ! curl -fsSL --retry 5 --retry-delay 10 -o "${tmp}/${CONNECT_TARBALL}" "${CONNECT_URL}"; then
    echo "  ERROR: download failed: ${CONNECT_URL}"
    echo "  IBM rotates these paths; check https://www.ibm.com/products/aspera/downloads"
    echo "  and re-run with CONNECT_URL=... CONNECT_SHA256= set."
    return 1
  fi

  if [[ -n "${CONNECT_SHA256}" ]]; then
    local actual
    actual="$(sha256sum "${tmp}/${CONNECT_TARBALL}" 2>/dev/null | awk '{print $1}')"
    [[ -z "${actual}" ]] && actual="$(shasum -a 256 "${tmp}/${CONNECT_TARBALL}" | awk '{print $1}')"
    if [[ "${actual}" != "${CONNECT_SHA256}" ]]; then
      echo "  ERROR: checksum mismatch for ${CONNECT_TARBALL}"
      echo "    expected ${CONNECT_SHA256}"
      echo "    got      ${actual}"
      echo "  If you intentionally changed CONNECT_VERSION, set CONNECT_SHA256= to skip this."
      return 1
    fi
    echo "  Checksum verified."
  fi

  tar xzf "${tmp}/${CONNECT_TARBALL}" -C "${tmp}"
  local installer=("${tmp}"/ibm-aspera-connect*.sh)
  if (( ${#installer[@]} == 0 )); then
    echo "  ERROR: no installer script inside ${CONNECT_TARBALL}"
    return 1
  fi

  # The installer always targets $HOME/.aspera/connect, so point HOME at the
  # work-dir root to keep everything self-contained and shared across nodes.
  mkdir -p "${ASPERA_HOME}"
  echo "  Running installer (non-interactive, no root required)..."
  if ! HOME="${ASPERA_HOME}" bash "${installer[0]}" >"${ASPERA_HOME}/install.log" 2>&1; then
    echo "  ERROR: installer failed; see ${ASPERA_HOME}/install.log"
    return 1
  fi
  # Desktop-integration warnings on a headless node are expected and harmless;
  # the binary landing is what actually matters.
  if [[ ! -x "${ASPERA_BIN}" ]]; then
    echo "  ERROR: ascp not found at ${ASPERA_BIN} after install"
    echo "  See ${ASPERA_HOME}/install.log"
    return 1
  fi
}

if [[ -x "${ASPERA_BIN}" ]]; then
  echo "[1/2] ascp already installed: ${ASPERA_BIN}"
else
  mkdir -p "${WORK_DIR}"
  install_connect || exit 1
fi
echo "      version: $("${ASPERA_BIN}" -A 2>&1 | head -1)"

# ── Step 2: the anonymous ENA key ───────────────────────────────────────────

# Emit every plausible location for asperaweb_id_dsa.openssh, cheapest first.
key_candidates() {
  [[ -n "${ASPERA_SSH_KEY:-}" ]] && echo "${ASPERA_SSH_KEY}"
  echo "${ASPERA_KEY_DEST}"
  echo "${HOME}/.aspera/connect/etc/asperaweb_id_dsa.openssh"
  [[ -n "${CONDA_PREFIX:-}" ]] && echo "${CONDA_PREFIX}/etc/asperaweb_id_dsa.openssh"

  # Relative to any ascp already on PATH — this is what catches an HPC
  # 'module load aspera' pointing at a pre-4.2 install.
  local on_path resolved
  on_path="$(command -v ascp 2>/dev/null || true)"
  if [[ -n "${on_path}" ]]; then
    resolved="$(readlink -f "${on_path}" 2>/dev/null || echo "${on_path}")"
    echo "$(dirname "$(dirname "${resolved}")")/etc/asperaweb_id_dsa.openssh"
  fi

  # Common site-module and conda roots.
  local root
  for root in /common/software/install/migrated/aspera /opt/aspera /usr/local/aspera \
              /usr/share/aspera "${HOME}/miniconda3/envs" "${HOME}/miniforge3/envs"; do
    local p
    for p in "${root}"/*/etc/asperaweb_id_dsa.openssh "${root}"/etc/asperaweb_id_dsa.openssh; do
      echo "${p}"
    done
  done
}

if [[ -s "${ASPERA_KEY_DEST}" ]]; then
  echo "[2/2] Key already cached: ${ASPERA_KEY_DEST}"
else
  FOUND_KEY=""
  while IFS= read -r cand; do
    [[ -n "${cand}" && -s "${cand}" ]] || continue
    if grep -q "PRIVATE KEY" "${cand}" 2>/dev/null; then
      FOUND_KEY="${cand}"
      break
    fi
  done < <(key_candidates)

  if [[ -n "${FOUND_KEY}" ]]; then
    mkdir -p "$(dirname "${ASPERA_KEY_DEST}")"
    cp "${FOUND_KEY}" "${ASPERA_KEY_DEST}"
    chmod 600 "${ASPERA_KEY_DEST}"
    echo "[2/2] Key found: ${FOUND_KEY}"
    echo "      cached to: ${ASPERA_KEY_DEST}"
  else
    cat <<EOF
[2/2] Could not locate 'asperaweb_id_dsa.openssh'.

  Connect ${CONNECT_VERSION} no longer ships it, and it is not hosted
  standalone anywhere authoritative — it exists only inside pre-4.2 installs.
  ENA still accepts it, so any of these will do:

    a) An Aspera module already on this cluster (most likely):
         module load aspera
         ls \$(dirname \$(dirname \$(readlink -f \$(command -v ascp))))/etc/
       then re-run this script, or pass the path directly:
         ASPERA_SSH_KEY=/path/to/asperaweb_id_dsa.openssh bash install_aspera.sh

    b) Conda, which bundles it:
         conda create -n asperakey -c hcc aspera-cli
         cp \$(conda run -n asperakey printenv CONDA_PREFIX)/etc/asperaweb_id_dsa.openssh \\
            ${ASPERA_KEY_DEST}

  Until then the pipeline uses its parallel-HTTPS paths, which are already
  several times faster than single-stream wget. Set USE_ASPERA=0 to skip the
  ascp attempt entirely and keep the logs clean.
EOF
    exit 2
  fi
fi

# ── Verify end to end ───────────────────────────────────────────────────────

echo "Verifying against ENA with a real transfer..."
VERIFY_DIR="$(mktemp -d "${WORK_DIR}/.aspera-verify.XXXXXX")"
trap 'rm -rf "${VERIFY_DIR}"' EXIT
# A small public CRAI, so success/failure is quick and unambiguous.
if "${ASPERA_BIN}" -i "${ASPERA_KEY_DEST}" -T -Q -l "${ASPERA_BANDWIDTH}" -P"${ASPERA_PORT}" \
     "${ENA_ASPERA_USER}:vol1/run/ERR324/ERR3242211/HG02121.final.cram.crai" \
     "${VERIFY_DIR}/" >"${VERIFY_DIR}/verify.log" 2>&1; then
  echo "  SUCCESS — Aspera transfer works."
  echo
  echo "Add these to your environment (or let config.sh auto-detect them):"
  echo "  export ASPERA_BIN=\"${ASPERA_BIN}\""
  echo "  export ASPERA_SSH_KEY=\"${ASPERA_KEY_DEST}\""
else
  echo "  FAILED — transfer did not complete. Log:"
  sed 's/^/    /' "${VERIFY_DIR}/verify.log" | tail -20
  echo
  echo "  The pipeline will fall back to its parallel-HTTPS paths."
  exit 1
fi
