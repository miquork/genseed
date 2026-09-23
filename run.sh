#!/bin/bash
# genseed - Copyright 2026 Mikko Voutilainen - Apache License 2.0
# run.sh - the two passes of genseed on N cores.
#
#   ./run.sh p1              # pass 1: all events, no correction -> rootfiles/GenSeed_p1.root
#   root -l -b -q 'drawGenSeed.C+("rootfiles/GenSeed_p1.root","_v1")'   # -> text/jec_v1.txt
#   ./run.sh p2              # pass 2: all events with text/jec_v1.txt -> rootfiles/GenSeed_v1.root
#   ./run.sh p1 4            # with 4 jobs instead of 8
#   LIST=data/other.txt NEV=2000 ./run.sh p1 2     # another list, a cap on the events
#
# The library is built ONCE with a two-event run before the jobs start, so
# that N jobs do not all decide it is stale and compile over one another; no
# source or header is to be edited until every job has printed "wrote".  The
# jobs read disjoint, even-aligned entry ranges of the input (see
# runGenSeed.C), so the a/b halves and their counters land in every job and
# hadd sums them consistently; hadd merges the TTree "genseed" too, so the
# merged file carries the tuple of every job.  The merge globs only the
# _job*of<n> files, so a stale merged file cannot be counted twice; the job
# files of an earlier run are removed before the launch, and nothing is
# merged unless every job printed "wrote", so a job that dies early cannot
# bring back its predecessor's file either.  LIST and the pass-2 table are
# relative to the package directory, wherever run.sh is started from.
set -u
PASS="${1:-p1}"
NJ="${2:-8}"
LIST="${LIST:-files.txt}"
NEV="${NEV:-0}"
cd "$(dirname "$0")" || exit 1
case "$PASS" in
  p1) TAG=_p1; JEC="" ;;
  p2) TAG=_v1; JEC=text/jec_v1.txt; [ -s "$JEC" ] || { echo "no $JEC in $(pwd) - run pass 1 and drawGenSeed first"; exit 1; } ;;
  *)  echo "usage: $0 p1|p2 [njobs]"; exit 1 ;;
esac
[ -e "$LIST" ] || { echo "no input list $LIST (one NanoAOD path per line, or a .root file; relative to $(pwd))"; exit 1; }
mkdir -p rootfiles text
echo "== build once"
root -l -b -q "runGenSeed.C(\"$LIST\",\"_build\",0,1,2)" > rootfiles/build.log 2>&1 \
  || { echo "build failed, see rootfiles/build.log"; tail -20 rootfiles/build.log; exit 1; }
grep -q "wrote rootfiles" rootfiles/build.log || { echo "build run did not finish, see rootfiles/build.log"; tail -20 rootfiles/build.log; exit 1; }
grep "clustering with" rootfiles/build.log
rm -f rootfiles/GenSeed_build.root
echo "== launch $NJ jobs, pass $PASS, tag $TAG, list $LIST"
rm -f rootfiles/GenSeed${TAG}_job*of${NJ}.root
for ((j=0; j<NJ; ++j)); do
  nohup root -l -b -q "runGenSeed.C(\"$LIST\",\"$TAG\",$j,$NJ,$NEV,\"$JEC\")" \
    > rootfiles/fill${TAG}_job${j}of${NJ}.log 2>&1 &
done
wait
echo "== jobs finished"
NFAIL=0
for ((j=0; j<NJ; ++j)); do
  grep -q "wrote rootfiles" rootfiles/fill${TAG}_job${j}of${NJ}.log \
    || { echo "  job $j did not finish cleanly, see rootfiles/fill${TAG}_job${j}of${NJ}.log"; NFAIL=$((NFAIL+1)); }
done
[ "$NFAIL" -eq 0 ] || { echo "$NFAIL of $NJ jobs failed - nothing merged"; exit 1; }
if [ "$NJ" -gt 1 ]; then
  hadd -f rootfiles/GenSeed${TAG}.root rootfiles/GenSeed${TAG}_job*of${NJ}.root > rootfiles/hadd${TAG}.log 2>&1 \
    && echo "  merged rootfiles/GenSeed${TAG}.root" || { echo "hadd failed"; tail rootfiles/hadd${TAG}.log; exit 1; }
fi
root -l -b -q -e "TFile f(\"rootfiles/GenSeed${TAG}.root\"); TH1 *h=(TH1*)f.Get(\"hist/hcount\"); for(int i=1;i<=h->GetNbinsX();++i) printf(\"  %-10s %12.0f\n\", h->GetXaxis()->GetBinLabel(i), h->GetBinContent(i)); TTree *t=(TTree*)f.Get(\"genseed\"); printf(\"  %-10s %12lld\n\", \"tuple\", t ? t->GetEntries() : 0LL);" 2>/dev/null
grep -h "ms/event" rootfiles/fill${TAG}_job*of${NJ}.log | grep done | tail -$NJ
