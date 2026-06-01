#!/usr/bin/env bash
set -euo pipefail

# --- Defaults ---
BIN_DEFAULT="/app/bin/z"            # path inside container
CONF_DEFAULT="/app/conf/z.conf"     # default config path

echo "---------------------------------------------"
echo "BIN         : $BIN_DEFAULT"
echo "CONF        : $CONF_DEFAULT"
echo "ARGS        : $@"
echo "---------------------------------------------"

if [[ ! -x "$BIN_DEFAULT" ]]; then
  echo "ERROR: binary not found or not executable: $BIN" >&2
  exit 2
fi
if [[ ! -f "$CONF_DEFAULT" ]]; then
  echo "WARN : conf file not found at $CONF"
fi

# Run + wall-time
START=$(date +%s)
"$BIN_DEFAULT" "$@"
RC=$?
END=$(date +%s)
echo "Exit code: $RC"
echo "Wall time: $((END-START)) s"
exit $RC
