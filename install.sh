#!/usr/bin/env bash
#
# install.sh — one-shot dependency setup for the rrna-operon pipeline.
#
# Installs, idempotently (anything already present is detected and skipped):
#   1. a conda/mamba package manager. mamba is preferred when available
#      (its solver is much faster); if neither mamba nor conda is found,
#      Miniforge is installed to ~/miniforge3 (it bundles mamba + conda-forge).
#   2. barrnap (+ aragorn/infernal/diamond) in a dedicated 'barrnap' conda env,
#      used for the optional rRNA/tRNA operon annotation step (run with --fast).
#   3. the DuckDB v1.5.3 CLI + the 'miint' extension.
#
# Usage:
#   ./install.sh                          # install everything
#   ./install.sh --no-barrnap             # DuckDB + miint only (skip conda/barrnap)
#   DUCKDB_DEST=~/.local/bin ./install.sh # change where the duckdb CLI lands
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- pinned versions (keep in lockstep with .github/workflows/test.yml) ------
DUCKDB_VERSION="v1.5.3"
MIINT_REPO_URL="https://ftp.microbio.me/pub/miint"
BARRNAP_ENV="barrnap"
DUCKDB_DEST="${DUCKDB_DEST:-${SCRIPT_DIR}/bin}"   # where the duckdb CLI is placed

c_grn=$'\033[1;32m'; c_yel=$'\033[1;33m'; c_red=$'\033[1;31m'; c_off=$'\033[0m'
log()  { printf '%s[install]%s %s\n' "$c_grn" "$c_off" "$*"; }
warn() { printf '%s[install]%s %s\n' "$c_yel" "$c_off" "$*" >&2; }
die()  { printf '%s[install] ERROR:%s %s\n' "$c_red" "$c_off" "$*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1; }

INSTALL_BARRNAP=1
while [ $# -gt 0 ]; do
  case "$1" in
    --no-barrnap) INSTALL_BARRNAP=0; shift ;;
    -h|--help)    sed -n '2,20p' "$0"; exit 0 ;;
    *)            die "unknown option: $1 (supported: --no-barrnap)" ;;
  esac
done
BN_VER=""

need curl || die "required tool 'curl' not found on PATH"
OS="$(uname -s)"; ARCH="$(uname -m)"

# one scratch dir for downloads (project-local, never /tmp), cleaned on exit
WORK="$(mktemp -d "${SCRIPT_DIR}/.install.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

if [ "$INSTALL_BARRNAP" = 1 ]; then
  #########################################################################
  # 1. conda / mamba
  #########################################################################
  source_conda_from() {                 # $1 = base dir; source conda+mamba if present
    [ -f "$1/etc/profile.d/conda.sh" ] || return 1
    # shellcheck disable=SC1091
    . "$1/etc/profile.d/conda.sh"
    [ -f "$1/etc/profile.d/mamba.sh" ] && . "$1/etc/profile.d/mamba.sh"
    return 0
  }
  resolve_pkgmgr() {                    # echo 'mamba' or 'conda' (preferring mamba), or ''
    need mamba && { echo mamba; return; }
    need conda && { echo conda; return; }
    local base
    for base in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3" /opt/conda; do
      if source_conda_from "$base"; then
        need mamba && { echo mamba; return; }
        need conda && { echo conda; return; }
      fi
    done
    echo ""
  }

  PKG="$(resolve_pkgmgr)"
  if [ -z "$PKG" ]; then
    log "No mamba/conda found — installing Miniforge (bundles mamba + conda-forge)"
    case "$OS" in Linux|Darwin) : ;; *) die "auto-install unsupported on OS '$OS'; install Miniforge manually" ;; esac
    mf="Miniforge3-${OS}-${ARCH}.sh"
    log "Downloading $mf …"
    curl -fsSL "https://github.com/conda-forge/miniforge/releases/latest/download/${mf}" -o "$WORK/mf.sh" \
      || die "Miniforge download failed for ${OS}-${ARCH}"
    bash "$WORK/mf.sh" -b -p "$HOME/miniforge3"
    source_conda_from "$HOME/miniforge3" || die "Miniforge installed but conda.sh not found"
    PKG="$(need mamba && echo mamba || echo conda)"
  fi
  log "Package manager: $PKG ($("$PKG" --version | head -1))"

  #########################################################################
  # 2. barrnap conda env
  #########################################################################
  if "$PKG" env list | awk '{print $1}' | grep -qx "$BARRNAP_ENV"; then
    log "conda env '$BARRNAP_ENV' already exists — skipping create"
  else
    log "Creating '$BARRNAP_ENV' env (barrnap from bioconda; pulls aragorn/infernal/diamond)…"
    "$PKG" create -y -n "$BARRNAP_ENV" -c conda-forge -c bioconda barrnap \
      || die "barrnap env creation failed"
  fi
  "$PKG" run -n "$BARRNAP_ENV" barrnap --version >/dev/null 2>&1 \
    || die "barrnap installed but won't run"
  BN_VER="$("$PKG" run -n "$BARRNAP_ENV" barrnap --version 2>&1 | head -1)"
  log "barrnap OK: $BN_VER"
  "$PKG" run -n "$BARRNAP_ENV" barrnap --help 2>&1 | grep -q -- '--trna' \
    || warn "this barrnap lacks --trna (need >=1.10 for tRNA annotation)"
else
  log "Skipping conda/barrnap (--no-barrnap); installing DuckDB + miint only"
fi

#############################################################################
# 3. DuckDB v1.5.3 CLI + miint extension
#############################################################################
duckdb_ok() { [ -x "$1" ] && [ "$("$1" --version 2>/dev/null | awk '{print $1}')" = "$DUCKDB_VERSION" ]; }

DUCKDB_BIN=""
for cand in "$DUCKDB_DEST/duckdb" "$(command -v duckdb 2>/dev/null || true)"; do
  if [ -n "$cand" ] && duckdb_ok "$cand"; then DUCKDB_BIN="$cand"; break; fi
done

if [ -n "$DUCKDB_BIN" ]; then
  log "DuckDB $DUCKDB_VERSION already present: $DUCKDB_BIN"
else
  case "${OS}-${ARCH}" in
    Linux-x86_64)              asset="duckdb_cli-linux-amd64.zip" ;;
    Linux-aarch64|Linux-arm64) asset="duckdb_cli-linux-arm64.zip" ;;
    Darwin-*)                  asset="duckdb_cli-osx-universal.zip" ;;
    *) die "no prebuilt DuckDB asset for ${OS}-${ARCH}; install ${DUCKDB_VERSION} manually and re-run" ;;
  esac
  need unzip || die "'unzip' is required to install the DuckDB CLI"
  mkdir -p "$DUCKDB_DEST"
  url="https://github.com/duckdb/duckdb/releases/download/${DUCKDB_VERSION}/${asset}"
  log "Downloading DuckDB ${DUCKDB_VERSION} ($asset)…"
  curl -fsSL "$url" -o "$WORK/duckdb.zip" || die "DuckDB download failed: $url"
  unzip -qo "$WORK/duckdb.zip" -d "$WORK"
  install -m 0755 "$WORK/duckdb" "$DUCKDB_DEST/duckdb"
  DUCKDB_BIN="$DUCKDB_DEST/duckdb"
  log "Installed DuckDB CLI → $DUCKDB_BIN"
fi

log "Installing the 'miint' extension (unsigned, from $MIINT_REPO_URL)…"
"$DUCKDB_BIN" -unsigned \
    -c "INSTALL miint FROM '${MIINT_REPO_URL}'; LOAD miint; SELECT 'miint loaded' AS status;" \
  || die "miint install failed — confirm ${MIINT_REPO_URL}/${DUCKDB_VERSION}/ publishes a build for ${OS}-${ARCH} and the host is reachable"
log "miint extension installed and loads cleanly"

#############################################################################
# done
#############################################################################
printf '\n%s[install] All dependencies ready.%s\n\n' "$c_grn" "$c_off"
printf '  DuckDB CLI : %s  (%s, miint loaded)\n' "$DUCKDB_BIN" "$DUCKDB_VERSION"
[ "$INSTALL_BARRNAP" = 1 ] && printf "  barrnap    : conda env '%s' (%s)\n" "$BARRNAP_ENV" "$BN_VER"
printf '\nRun the pipeline (point it at the duckdb CLI just installed):\n\n'
printf '  DUCKDB="%s" ./run.sh --output out/ reads.fastq.gz\n' "$DUCKDB_BIN"
printf '  # or put it on PATH:  export PATH="%s:$PATH"\n' "$DUCKDB_DEST"
if [ "$INSTALL_BARRNAP" = 1 ]; then
  printf '\nAnnotate operon consensus (rRNA + ITS tRNA, fast mode):\n\n'
  printf '  conda run -n %s barrnap --kingdom bac --trna --fast out/consensus.fa > out/operons.gff\n' "$BARRNAP_ENV"
fi
printf '\n'
