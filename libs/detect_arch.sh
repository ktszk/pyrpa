#!/bin/sh
# detect_arch.sh -- list the -march targets needed by the machines this library
# will actually RUN on (not the one it is built on).
#
#   sh detect_arch.sh            one target per line, sorted, deduplicated
#   sh detect_arch.sh --table    node / CPU model / target, for the build log
#   sh detect_arch.sh --local    probe only this machine
#   sh detect_arch.sh --refresh  ignore the cache and re-probe
#
# WHY this exists: `-march=native` targets the BUILD host.  Building on a Zen 4
# node (avx512) and running on a Zen 3 one (avx2) gives SIGILL, so the target
# set has to come from the cluster's node inventory.  `make auto` builds one
# libfmod_<target>.so per line printed here and libs/flibs/_loader.py picks the
# best one the running CPU supports at import time.
#
# Detection is delegated to gcc itself (`gcc -march=native -Q --help=target`),
# which knows the CPU tables; the /proc/cpuinfo flag fallback below is only for
# nodes without gcc.  Hand-written flag tables are brittle -- this kernel spells
# Zen 4's flags avx512_vbmi2/avx512_vnni, other kernels use avx512vbmi2.

set -u
CACHE="$(dirname "$0")/.arch_cache"
MODE=list
REFRESH=0
for a in "$@"; do
  case "$a" in
    --table)   MODE=table ;;
    --local)   MODE=${MODE}; LOCAL_ONLY=1 ;;
    --refresh) REFRESH=1 ;;
    -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
  esac
done
LOCAL_ONLY=${LOCAL_ONLY:-0}

# ---- the probe, run ON each node (stdin of sh -s) --------------------------
PROBE='
  m=$(awk -F: "/model name/{sub(/^ /,\"\",\$2); print \$2; exit}" /proc/cpuinfo)
  a=$(gcc -march=native -Q --help=target 2>/dev/null |
      awk "/^ *-march=/{gsub(/[ \t]/,\"\",\$0); sub(/^-march=/,\"\"); print; exit}")
  if [ -z "$a" ] || [ "$a" = "native" ]; then     # no gcc: fall back to the flags
    f=$(awk "/^flags/{print; exit}" /proc/cpuinfo)
    case " $f " in
      *" avx512f "*) case " $f " in *znver*|*" avx512_bf16 "*|*" avx512bf16 "*) a=znver4 ;;
                                    *) a=skylake-avx512 ;; esac ;;
      *" avx2 "*)    a=znver3 ;;
      *)             a=x86-64-v2 ;;
    esac
  fi
  echo "${m:-unknown}|$a"
'

probe_local() { printf "%s" "$PROBE" | sh -s; }

probe_node() {   # $1 = node name; try the scheduler first, then ssh
  _part=$(sinfo -h -N -o '%n %P' 2>/dev/null | sed 's/\*//' |
          awk -v n="$1" '$1==n{print $2; exit}')
  if [ -n "${_part:-}" ] &&
     out=$(printf "%s" "$PROBE" |
           timeout 60 srun -p "$_part" -w "$1" -n1 --immediate=30 -Q sh -s 2>/dev/null) &&
     [ -n "$out" ]; then
    echo "$out"; return 0
  fi
  out=$(printf "%s" "$PROBE" |
        timeout 30 ssh -o BatchMode=yes -o ConnectTimeout=5 "$1" sh -s 2>/dev/null)
  [ -n "$out" ] && { echo "$out"; return 0; }
  return 1
}

list_nodes() {
  if command -v sinfo >/dev/null 2>&1; then
    sinfo -h -N -o '%n' 2>/dev/null | sort -u
  elif command -v pbsnodes >/dev/null 2>&1; then
    pbsnodes -l all 2>/dev/null | awk '{print $1}' | sort -u
  fi
}

# ---- collect ---------------------------------------------------------------
if [ "$LOCAL_ONLY" = 1 ]; then
  NODES=$(hostname -s)
else
  NODES=$(list_nodes)
  [ -n "${NODES:-}" ] || NODES=$(hostname -s)
fi
KEY=$(echo "$NODES" | tr '\n' ' ')

if [ "$REFRESH" = 0 ] && [ -f "$CACHE" ] &&
   [ "$(head -1 "$CACHE")" = "#nodes: $KEY" ]; then
  TABLE=$(tail -n +2 "$CACHE")
else
  TABLE=""
  for n in $NODES; do
    if [ "$n" = "$(hostname -s)" ] && [ "$LOCAL_ONLY" = 1 ]; then
      r=$(probe_local)
    else
      r=$(probe_node "$n") || r=$(: )
    fi
    if [ -z "${r:-}" ]; then
      # unreachable node: do NOT guess, and do not silently drop it either
      echo "detect_arch: WARNING: could not probe $n -- keeping the portable default" >&2
      r="unreachable|znver3"
    fi
    TABLE="${TABLE}${n}|${r}
"
  done
  TABLE=$(echo "$TABLE" | sed '/^$/d')
  { echo "#nodes: $KEY"; echo "$TABLE"; } > "$CACHE" 2>/dev/null || true
fi

# ---- report ----------------------------------------------------------------
if [ "$MODE" = table ]; then
  echo "$TABLE" | awk -F'|' 'BEGIN{printf "%-12s %-40s %s\n","NODE","CPU","-march"}
                             {printf "%-12s %-40s %s\n",$1,$2,$3}'
else
  echo "$TABLE" | awk -F'|' '{print $3}' | sort -u
fi
