#!/usr/bin/env bash
set -euo pipefail

# --- Defaults ---
BIN_DEFAULT="/app/bin/z"            # path inside container
CONF_DEFAULT="/app/conf/z.conf"     # default config path
SUB_DEFAULT="z"                     # subcommand: z | z_os | valley_scanner | refine
T_DEFAULT="4400123.227930"          # t (just a test)

usage() {
  cat <<EOF
Usage: ./run.sh [options] [SUBCOMMAND] [t]
  SUBCOMMAND        default: ${SUB_DEFAULT}   (e.g. z, z_os)
  t                 default: ${T_DEFAULT}

Options:
  -b, --bin PATH        Path to binary (default ${BIN_DEFAULT})
  -c, --conf PATH       Path to conf (default ${CONF_DEFAULT})
  -j, --threads N       OMP threads (default from OMP_NUM_THREADS or 2)
  --                    End of options; remaining args passed as SUBCOMMAND and t

Environment:
  OMP_NUM_THREADS       Overrides threads if set
  OMP_PROC_BIND         Default: close
  OMP_PLACES            Default: cores
  Z_CONF                If set, passed as conf path
Examples:
  ./run.sh z_os 4.400123411112e11
  ./run.sh -j 2 z 440012341111.2279308671282884
EOF
}

BIN="${BIN_DEFAULT}"
CONF="${Z_CONF:-$CONF_DEFAULT}"
SUB="${SUB_DEFAULT}"
TVAL="${T_DEFAULT}"
THREADS="${OMP_NUM_THREADS:-2}"

# --- Parse args ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    -b|--bin) BIN="$2"; shift 2 ;;
    -c|--conf) CONF="$2"; shift 2 ;;
    -j|--threads) THREADS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    --) shift; break ;;
    -*)
      echo "Unknown option: $1" >&2; usage; exit 1 ;;
    *)
      # positional: SUBCOMMAND and t
      SUB="$1"; shift
      if [[ $# -gt 0 ]]; then TVAL="$1"; shift; fi
      ;;
  esac
done

export OMP_NUM_THREADS="${THREADS}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-close}"
export OMP_PLACES="${OMP_PLACES:-cores}"

echo "---------------------------------------------"
echo "BIN         : $BIN"
echo "SUBCOMMAND  : $SUB"
echo "t           : $TVAL"
echo "CONF        : $CONF"
echo "Threads     : $OMP_NUM_THREADS"
echo "Bind/Places : $OMP_PROC_BIND / $OMP_PLACES"
echo "---------------------------------------------"

if [[ ! -x "$BIN" ]]; then
  echo "ERROR: binary not found or not executable: $BIN" >&2
  exit 2
fi
if [[ ! -f "$CONF" ]]; then
  echo "WARN : conf file not found at $CONF"
fi

# Run + wall-time
START=$(date +%s)
"$BIN" "$SUB" "$TVAL"
RC=$?
END=$(date +%s)
echo "Exit code: $RC"
echo "Wall time: $((END-START)) s"
exit $RC
