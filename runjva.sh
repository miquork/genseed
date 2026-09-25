#!/bin/bash
# genseed - Copyright 2026 Mikko Voutilainen - Apache License 2.0
# runjva.sh - the JVA study (JVA.h, jvassoc.h) on N cores, as run.sh does genseed.
#
#   ./runjva.sh v1                   # all events, 8 jobs, text/jec_v1.txt -> rootfiles/JVA_v1.root
#   ./runjva.sh v1 4                 # 4 jobs
#   LIST=data/other.txt NEV=2000 ./runjva.sh test 2      # another list, a cap on the events
#   JEC= ./runjva.sh raw 8           # no correction (the dijet part is then uncorrected)
#   TAUTRK=2 TAUALL=8 ./runjva.sh tau2 8                 # other prices
#   CANDMODE=2 RECOILMODE=1 MHTMIN=10 ./runjva.sh c2r1 8  # v2: in-cluster charged candidates, MHT recoil
#
# Knobs (defaults = v1): TAUTRK 4, TAUALL 12, SIGMAK 1, TRKMINPT 0.5, TUPLE 1,
# CANDMODE 0 (0 in-cone tracks, 1 in-cluster charged, 2 either), CHMINPT 0.5,
# RECOILMODE 0 (0 particle MET, 1 MHT of the vertex-resolved jets), MHTMIN 10,
# SIGMAJ 0.35 (see jvassoc.h).
#
# The correction is required for the dijet part and defaults to the genseed
# pass-1 table text/jec_v1.txt (made by drawGenSeed.C, see run.sh); a table
# that is named but missing stops the script before anything runs.  The
# library is built ONCE with a two-event run before the jobs start (into
# $TMPDIR/aclic_genseed/<backend>, see runJVA.C), so that N jobs do not all
# decide it is stale and compile over one another; no source or header is to
# be edited until every job has printed "wrote".  The jobs read disjoint,
# even-aligned entry ranges (an empty range writes a valid empty file), the
# merge globs only the _job*of<n> files of this tag, the job files of an
# earlier run with the same tag are removed before the launch, and nothing is
# merged unless every job printed "wrote".  LIST and JEC are relative to the
# package directory, wherever runjva.sh is started from.
set -u
TAG="${1:-v1}"
NJ="${2:-8}"
LIST="${LIST:-files.txt}"
NEV="${NEV:-0}"
JEC="${JEC-text/jec_v1.txt}"
TAUTRK="${TAUTRK:-4}"
TAUALL="${TAUALL:-12}"
SIGMAK="${SIGMAK:-1}"
TRKMINPT="${TRKMINPT:-0.5}"
TUPLE="${TUPLE:-1}"
CANDMODE="${CANDMODE:-0}"
CHMINPT="${CHMINPT:-0.5}"
RECOILMODE="${RECOILMODE:-0}"
MHTMIN="${MHTMIN:-10}"
SIGMAJ="${SIGMAJ:-0.35}"
cd "$(dirname "$0")" || exit 1
TAG="_${TAG#_}"
[ -e "$LIST" ] || { echo "no input list $LIST (one NanoAOD path per line, or a .root file; relative to $(pwd))"; exit 1; }
if [ -n "$JEC" ]; then
  [ -s "$JEC" ] || { echo "no $JEC in $(pwd) - the dijet part needs the genseed correction (run.sh p1 + drawGenSeed), or set JEC= to run without"; exit 1; }
else
  echo "WARNING: no correction, the dijet part is uncorrected"
fi
case "$NJ" in ''|*[!0-9]*) echo "usage: $0 tag [njobs]"; exit 1 ;; esac
[ "$NJ" -ge 1 ] || { echo "usage: $0 tag [njobs >= 1]"; exit 1; }
mkdir -p rootfiles
ARGS="$TAUTRK,$TAUALL,$SIGMAK,$TRKMINPT,$TUPLE,$CANDMODE,$CHMINPT,$RECOILMODE,$MHTMIN,$SIGMAJ"
echo "== build once"
root -l -b -q "runJVA.C(\"$LIST\",\"${TAG}_build\",0,1,2,\"$JEC\",$ARGS)" > rootfiles/buildjva${TAG}.log 2>&1 \
  || { echo "build failed, see rootfiles/buildjva${TAG}.log"; tail -20 rootfiles/buildjva${TAG}.log; exit 1; }
grep -q "wrote rootfiles" rootfiles/buildjva${TAG}.log || { echo "build run did not finish, see rootfiles/buildjva${TAG}.log"; tail -20 rootfiles/buildjva${TAG}.log; exit 1; }
grep "clustering with" rootfiles/buildjva${TAG}.log
rm -f rootfiles/JVA${TAG}_build.root
echo "== launch $NJ jobs, tag $TAG, list $LIST, correction ${JEC:-none}, tauTrk $TAUTRK tauAll $TAUALL sigmaK $SIGMAK trkMinPt $TRKMINPT, candMode $CANDMODE chMinPt $CHMINPT, recoilMode $RECOILMODE mhtMin $MHTMIN sigmaJ $SIGMAJ"
rm -f rootfiles/JVA${TAG}_job*of${NJ}.root
if [ "$NJ" -eq 1 ]; then
  root -l -b -q "runJVA.C(\"$LIST\",\"$TAG\",0,1,$NEV,\"$JEC\",$ARGS)" > rootfiles/filljva${TAG}_job0of1.log 2>&1
else
  for ((j=0; j<NJ; ++j)); do
    nohup root -l -b -q "runJVA.C(\"$LIST\",\"$TAG\",$j,$NJ,$NEV,\"$JEC\",$ARGS)" \
      > rootfiles/filljva${TAG}_job${j}of${NJ}.log 2>&1 &
  done
  wait
fi
echo "== jobs finished"
NFAIL=0
for ((j=0; j<NJ; ++j)); do
  grep -q "wrote rootfiles" rootfiles/filljva${TAG}_job${j}of${NJ}.log \
    || { echo "  job $j did not finish cleanly, see rootfiles/filljva${TAG}_job${j}of${NJ}.log"; NFAIL=$((NFAIL+1)); }
  grep -q "VIOLATED" rootfiles/filljva${TAG}_job${j}of${NJ}.log && echo "  job $j: an identity is VIOLATED, see its log"
done
[ "$NFAIL" -eq 0 ] || { echo "$NFAIL of $NJ jobs failed - nothing merged"; exit 1; }
if [ "$NJ" -gt 1 ]; then
  hadd -f rootfiles/JVA${TAG}.root rootfiles/JVA${TAG}_job*of${NJ}.root > rootfiles/haddjva${TAG}.log 2>&1 \
    && echo "  merged rootfiles/JVA${TAG}.root" || { echo "hadd failed"; tail rootfiles/haddjva${TAG}.log; exit 1; }
fi
# counters, and hopts divided by the jobs: every knob must come out as set
root -l -b -q -e "TFile f(\"rootfiles/JVA${TAG}.root\"); TH1 *h=(TH1*)f.Get(\"hist/hcount\"); for(int i=1;i<=h->GetNbinsX();++i) printf(\"  %-10s %12.0f\n\", h->GetXaxis()->GetBinLabel(i), h->GetBinContent(i)); TTree *t=(TTree*)f.Get(\"jva\"); printf(\"  %-10s %12lld\n\", \"tuple\", t ? t->GetEntries() : 0LL); TH1 *o=(TH1*)f.Get(\"hist/hopts\"); const double nj=h->GetBinContent(6); printf(\"  knobs (hopts / jobs):\"); for(int i=1;i<=o->GetNbinsX();++i) printf(\" %s=%g\", o->GetXaxis()->GetBinLabel(i), nj>0 ? o->GetBinContent(i)/nj : 0.); printf(\"\\n\");" 2>/dev/null
grep -h "ms/crossing" rootfiles/filljva${TAG}_job*of${NJ}.log | grep done | tail -$NJ
