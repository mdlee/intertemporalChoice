#!/usr/bin/env bash
# Append-only MCMC status watcher → models/logs/mcmcRuns.status.log
# Follow with:  tail -f models/logs/mcmcRuns.status.log
set -u

ROOT="$(cd "$(dirname "$0")" && pwd)"
LOG_DIR="${ROOT}/logs"
OUT="${LOG_DIR}/mcmcRuns.status.log"
PIDFILE="${LOG_DIR}/watchMcmcStatus.pid"
INTERVAL="${MCMC_STATUS_INTERVAL:-30}"

mkdir -p "$LOG_DIR"
echo $$ > "$PIDFILE"

stamp() { date '+%Y-%m-%d %H:%M:%S'; }
append() { printf '[%s] %s\n' "$(stamp)" "$*" >> "$OUT"; }

read_status() {
  local f="$1"
  [[ -f "$f" ]] || return 1
  local msg
  msg="$(sed -n '2p' "$f" 2>/dev/null | tr -d '\r')"
  [[ -n "$msg" ]] || return 1
  printf '%s' "$msg"
}

matlab_line() {
  local line pid etime
  # Prefer the real MATLAB binary (ignore Cursor sandbox shells that echo the same text)
  line="$(ps -u mdlee -o pid=,etime=,cmd= 2>/dev/null | grep -F 'glnxa64/MATLAB' | grep -F 'runLatentMixtureRobustnessSequential' | grep -v grep | head -1 || true)"
  if [[ -z "$line" ]]; then
    line="$(ps -u mdlee -o pid=,etime=,cmd= 2>/dev/null | grep -F 'glnxa64/MATLAB' | grep -F 'runHierarchicalExecutionSequential' | grep -v grep | head -1 || true)"
  fi
  if [[ -z "$line" ]]; then
    line="$(ps -u mdlee -o pid=,etime=,cmd= 2>/dev/null | grep -F 'glnxa64/MATLAB' | grep -F 'runLatentMixtureSequential' | grep -v Robustness | grep -v grep | head -1 || true)"
  fi
  if [[ -n "$line" ]]; then
    pid="$(awk '{print $1}' <<<"$line")"
    etime="$(awk '{print $2}' <<<"$line")"
    echo "matlab pid=$pid etime=$etime"
  else
    echo "matlab: none matching sequential runners"
  fi
}

last_batch_line() {
  local diary="$1"
  [[ -f "$diary" ]] || return 1
  grep -E '^\-\-\- .* \| batch=' "$diary" 2>/dev/null | tail -1 || true
}

last_key=""
declare -A LAST_BATCH_KEY=()
append "watcher started (interval=${INTERVAL}s pid=$$ user=$(id -un)) → $OUT"

while true; do
  job="unknown"
  msg=""
  if msg="$(read_status "${LOG_DIR}/runLatentMixtureRobustnessSequential.status")"; then
    # Prefer robustness while its snapshot does not say finished
    if [[ "$msg" == *finished* ]]; then
      if msg2="$(read_status "${LOG_DIR}/runHierarchicalExecutionSequential.status")"; then
        job="hierarchical"
        msg="$msg2"
      else
        job="robustness"
      fi
    else
      job="robustness"
    fi
  elif msg="$(read_status "${LOG_DIR}/runHierarchicalExecutionSequential.status")"; then
    job="hierarchical"
  elif msg="$(read_status "${LOG_DIR}/runLatentMixtureSequential.status")"; then
    job="mixture"
  else
    msg="(no .status snapshot)"
  fi

  ml="$(matlab_line)"
  key="${job}|${msg}|${ml}"
  if [[ "$key" != "$last_key" ]]; then
    append "[$job] $msg | $ml"
    last_key="$key"
  fi

  for diary in \
    "${LOG_DIR}/runLatentMixtureRobustnessSequential.log" \
    "${LOG_DIR}/runHierarchicalExecutionSequential.log" \
    "${LOG_DIR}/runLatentMixtureSequential.log"
  do
    batch_line="$(last_batch_line "$diary")"
    if [[ -n "${batch_line:-}" ]]; then
      prev="${LAST_BATCH_KEY[$diary]:-}"
      if [[ "$prev" != "$batch_line" ]]; then
        append "diary: $batch_line"
        LAST_BATCH_KEY[$diary]="$batch_line"
      fi
    fi
  done

  sleep "$INTERVAL" || break
done
