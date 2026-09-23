// genseed - Copyright 2026 Mikko Voutilainen - Apache License 2.0
// runGenSeed.C - the ROOT driver of GenSeed: knobs, job split by entry range,
// clustering backend, build directory.  Every numeric knob is set here and
// every one of them ends up in hist/hopts, so a histogram file always says how
// it was made.
//
// TWO PASSES, because the correction is derived from this sample:
//   root -l -b -q 'runGenSeed.C("files.txt","_p1")'
//   root -l -b -q 'drawGenSeed.C+("rootfiles/GenSeed_p1.root","_v1")'          // -> text/jec_v1.txt
//   root -l -b -q 'runGenSeed.C("files.txt","_v1",0,1,0,"text/jec_v1.txt")'
// `list` is a text file with one NanoAOD path per line (# comments), or a
// single .root file given directly.
//
// JOBS SPLIT BY ENTRY RANGE, NOT BY FILE.  The input is typically one large
// file, so a split by file would give job 0 everything, and a split by entry
// modulo would read every basket in every job.  Job k of n takes the
// contiguous range [first, last) of the requested entries, aligned to EVEN
// entries so that each job holds equal numbers of half a (even entries) and
// half b (odd entries) and hadd sums them consistently; the half is taken
// from the GLOBAL entry number.
//   root -l -b -q 'runGenSeed.C("files.txt","_p1",3,8)'                        // job 3 of 8, pass 1
//   root -l -b -q 'runGenSeed.C("files.txt","_v1",3,8,0,"text/jec_v1.txt")'    // job 3 of 8, pass 2
//   hadd -f rootfiles/GenSeed_p1.root rootfiles/GenSeed_p1_job*of8.root
// hcount 'jobs' counts the merged jobs, hist/hopts holds value x jobs, and
// the TTree "genseed" is merged by hadd like everything else.  With fewer
// than two events per job some ranges are empty; such a job writes a valid
// file with no events (it counts in 'jobs'), and the driver warns.
//
// BUILD THE LIBRARY FIRST IF YOU LAUNCH JOBS IN PARALLEL (run.sh does): N
// jobs starting at once all decide the library is out of date and compile
// on top of each other.  The library goes OUTSIDE the source tree, into
// $TMPDIR/aclic_genseed/<backend>: a shared library is memory-mapped, and a
// sync client rewriting it under a running process turns the next page fault
// into a bus error; one directory per backend because ACLiC decides on
// timestamps, not on compile flags, and would happily reuse a FastJet build
// for a tiledjet request.
//
// THE CLUSTERING BACKEND.  tiledjet.h (clean-room, tiled anti-kT) is the
// default.  FastJet is used, with -DGENSEED_USE_FASTJET, when an installation
// is found ($FASTJET, $CONDA_PREFIX, ~/miniforge3/envs/eejet, /usr/local,
// /opt/homebrew) AND tiledjet.h says it is not at least as fast
// (tiledjet::kAtLeastAsFastAsFastJet = false, the verdict of
// bench_tiledjet.C), or when tiledjet.h is not there at all.  Without either
// there is nothing to cluster with and the driver says so.  The choice is
// printed and recorded in hopts (useFastJet).
#include <TROOT.h>
#include <TChain.h>
#include <TSystem.h>
#include <TString.h>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

namespace {
  // The FastJet prefix, or "" if there is none to be had.
  TString FindFastJet()
  {
    std::vector<TString> tries;
    if (gSystem->Getenv("FASTJET")) tries.push_back(gSystem->Getenv("FASTJET"));
    if (gSystem->Getenv("CONDA_PREFIX")) tries.push_back(gSystem->Getenv("CONDA_PREFIX"));
    tries.push_back(TString(gSystem->Getenv("HOME")) + "/miniforge3/envs/eejet");
    tries.push_back("/usr/local");
    tries.push_back("/opt/homebrew");
    for (size_t i = 0; i < tries.size(); ++i) {
      if (tries[i].IsNull()) continue;
      if (!gSystem->AccessPathName(tries[i] + "/include/fastjet/ClusterSequence.hh")) return tries[i];
    }
    return "";
  }
  // The benchmark verdict compiled into tiledjet.h: 1 if tiledjet is at
  // least as fast as FastJet, 0 if not, -1 if the header is missing or does
  // not say.  Read as text: the header need not be interpretable by cling.
  // The verdict is set by hand, so the parse is strict: a trailing comment
  // is dropped first ("= false;  // true on the M2" is false), the token
  // right after the '=' (or '{') following the name must be exactly true or
  // false, and a line that only mentions the name (kPreferFastJet =
  // !kAtLeastAsFastAsFastJet) is not a definition.
  int TiledjetVerdict(const char *header)
  {
    std::ifstream in(header);
    if (!in) return -1;
    const std::string name = "kAtLeastAsFastAsFastJet";
    std::string line;
    while (std::getline(in, line)) {
      line = line.substr(0, line.find("//"));                  // the code part only
      const size_t k = line.find(name);
      if (k == std::string::npos) continue;
      size_t p = line.find_first_of("={", k + name.size());
      if (p == std::string::npos) continue;                    // a mention, not the definition
      p = line.find_first_not_of(" \t={", p);
      if (p == std::string::npos) continue;
      const size_t e = line.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_", p);
      const std::string tok = line.substr(p, e == std::string::npos ? std::string::npos : e - p);
      if (tok == "true")  return 1;
      if (tok == "false") return 0;
    }
    return -1;
  }
}

void runGenSeed(const char *list = "files.txt",
                const char *tag  = "_p1",
                int job = 0, int njobs = 1,
                double maxEvents = 0,
                const char *jecFile = "",
                double ptMinCluster = 1.0,
                int maxVertices = -1,
                double kernelScale = 1.0,
                int useDz = 1,
                int parallaxRef = 1,
                int bendSign = +1,
                double dzRecoSeed = 0.2)
{
   gSystem->mkdir("rootfiles", kTRUE);
   gSystem->mkdir("text", kTRUE);

   // ---- input -------------------------------------------------------------
   std::vector<TString> all;
   if (TString(list).EndsWith(".root")) all.push_back(list);
   else {
      std::ifstream in(list); TString line;
      while (line.ReadLine(in)) { line = line.Strip(TString::kBoth);
         if (line.Length() && !line.BeginsWith("#")) all.push_back(line); }
   }
   if (all.empty()) { printf("runGenSeed: no input files in %s\n", list); return; }
   if (njobs < 1 || job < 0 || job >= njobs) { printf("runGenSeed: job %d of %d makes no sense\n", job, njobs); return; }
   TChain *c = new TChain("Events");
   for (size_t i = 0; i < all.size(); ++i) c->Add(all[i]);
   const Long64_t nall = c->GetEntries();
   printf("runGenSeed: %zu files, %lld entries, root %s\n", all.size(), nall,
          gSystem->Getenv("ROOTSYS") ? gSystem->Getenv("ROOTSYS") : "(ROOTSYS unset)");
   if (!nall) return;

   // ---- the job's entry range, even-aligned --------------------------------
   const Long64_t nmax  = maxEvents > 0 ? std::min(nall, (Long64_t)maxEvents) : nall;
   const Long64_t first = 2*(Long64_t)floor(double(job)*nmax/(2.*njobs));
   const Long64_t last  = (job == njobs-1) ? nmax : 2*(Long64_t)floor(double(job+1)*nmax/(2.*njobs));
   printf("runGenSeed: job %d of %d, entries [%lld, %lld) of the %lld requested\n", job, njobs, first, last, nmax);
   // Ranges are whole event pairs, so with fewer than two events per job some
   // are empty: those jobs write an empty file (GenSeed::Loop stops at
   // last <= first) instead of reading the chain to its end.
   if (njobs > 1 && nmax < 2*(Long64_t)njobs)
      printf("runGenSeed: WARNING %lld events for %d jobs, fewer than two per job: some ranges are empty%s\n",
             nmax, njobs, last <= first ? ", this one included" : "");

   // ---- the clustering backend and the build directory ---------------------
   const TString src = gSystem->DirName(__FILE__);      // the package directory, wherever it is run from
   const TString fj = FindFastJet();
   const int verdict = TiledjetVerdict(src + "/tiledjet.h");
   bool useFJ = false;
   if (verdict < 0 && !gSystem->AccessPathName(src + "/tiledjet.h"))
      printf("runGenSeed: tiledjet.h does not state kAtLeastAsFastAsFastJet; treating it as false\n");
   if (verdict < 0 && gSystem->AccessPathName(src + "/tiledjet.h"))
      printf("runGenSeed: WARNING tiledjet.h is missing from %s\n", src.Data());
   if (!fj.IsNull() && verdict != 1) useFJ = true;
   if (!useFJ && gSystem->AccessPathName(src + "/tiledjet.h")) {
      printf("runGenSeed: neither tiledjet.h nor FastJet is available - nothing to cluster with.  Set $FASTJET.\n");
      return;
   }
   {
      const char *tmp = gSystem->Getenv("TMPDIR");
      const TString bdir = Form("%s/aclic_genseed/%s", tmp && tmp[0] ? tmp : "/tmp", useFJ ? "fastjet" : "tiledjet");
      gSystem->mkdir(bdir, kTRUE);
      gSystem->SetBuildDir(bdir, kTRUE);
      gSystem->AddIncludePath(Form("-I%s", src.Data()));
   }
   if (useFJ) {
      printf("runGenSeed: clustering with FastJet at %s (tiledjet verdict %d)\n", fj.Data(), verdict);
      gSystem->AddIncludePath(Form("-DGENSEED_USE_FASTJET -I%s/include", fj.Data()));
      gSystem->AddLinkedLibs(Form("-L%s/lib -Wl,-rpath,%s/lib -lfastjet", fj.Data(), fj.Data()));
      gSystem->AddDynamicPath(Form("%s/lib", fj.Data()));
   } else
      printf("runGenSeed: clustering with tiledjet.h (verdict %d, FastJet %s)\n", verdict,
             fj.IsNull() ? "not found" : "found but not used");

   const TString out = Form("rootfiles/GenSeed%s%s.root", tag, njobs > 1 ? Form("_job%dof%d", job, njobs) : "");

   if (gROOT->LoadMacro(src + "/GenSeed.C+") < 0) { printf("runGenSeed: build failed\n"); return; }
   gROOT->ProcessLine(Form("GenSeed a((TTree*)%p,\"%s\");", (void*)c, out.Data()));

   // ---- knobs ---------------------------------------------------------------
   gROOT->ProcessLine(Form("a.SetOpt(\"firstEntry\",%lld);", first));
   gROOT->ProcessLine(Form("a.SetOpt(\"lastEntry\",%lld);", last));
   // The correction table from pass one.  Empty means pass one itself.
   if (jecFile && jecFile[0]) gROOT->ProcessLine(Form("a.SetJEC(\"%s\");", jecFile));
   // One radius for the pure gen jets, the reco jets and the RecoJetSeeds.
   gROOT->ProcessLine("a.SetOpt(\"R\",0.4);");
   // Reco jets are clustered from this raw (PUPPI-weighted) pT, far below the
   // 5 GeV the measurement is reported from, so that the response is not
   // truncated where the spectrum starts and the 2-5 GeV bins have content.
   gROOT->ProcessLine(Form("a.SetOpt(\"ptMinCluster\",%g);", ptMinCluster));
   // Pure gen jets from 1 GeV, so that the 2-3 GeV bin is populated too.
   gROOT->ProcessLine("a.SetOpt(\"ptMinGen\",1.0);");
   // The tuple keeps jets above these; a partner below them is index -1.
   gROOT->ProcessLine("a.SetOpt(\"ptStoreReco\",3.0);");
   gROOT->ProcessLine("a.SetOpt(\"ptStoreGen\",3.0);");
   // Generated particles kept (neutrinos are always dropped).
   gROOT->ProcessLine("a.SetOpt(\"etaMaxPart\",5.5);");
   // A pure gen jet feeds a reco jet when it gives it more than this fraction
   // of the seed, and a reco jet takes a pure jet when it receives more than
   // this fraction of it.  1/3: two can share a jet, three cannot all be dominant.
   gROOT->ProcessLine("a.SetOpt(\"shareFrac\",0.333333333);");
   // A reco vertex OWNS the interaction nearest to it if that is within this in
   // z; the z recovery is good to 25-49 um, so 0.05 cm is loose on purpose.
   gROOT->ProcessLine("a.SetOpt(\"dzVertexGen\",0.05);");
   // An interaction is looked at on the gen side (topology, dR partner,
   // RecoJetSeed) at its nearest usable vertex if that is within this.
   gROOT->ProcessLine(Form("a.SetOpt(\"dzRecoSeed\",%g);", dzRecoSeed));
   // Vertices with fewer tracks than this have no usable z and are skipped.
   gROOT->ProcessLine("a.SetOpt(\"minTrkPerVertex\",3);");
   // Copies of the same jet from different vertex hypotheses: closer than this
   // in dR and in relative pT is one jet, kept once.
   gROOT->ProcessLine("a.SetOpt(\"dupDR\",0.15);");
   gROOT->ProcessLine("a.SetOpt(\"dupPtFrac\",0.10);");
   // The linker (genlink.h).  kernelScale scales every calo sigma and cone;
   // useDz turns the pass-1 dz cut off for the blind dz test; parallaxRef 1
   // references the calo directions to PV_z (link/hparallax must then have
   // slope +1); bendSign fixes the sign of the charged bend (link/hbend).
   gROOT->ProcessLine(Form("a.SetOpt(\"kernelScale\",%g);", kernelScale));
   gROOT->ProcessLine(Form("a.SetOpt(\"useDz\",%d);", useDz));
   gROOT->ProcessLine("a.SetOpt(\"useCharge\",1);");
   gROOT->ProcessLine(Form("a.SetOpt(\"parallaxRef\",%d);", parallaxRef));
   gROOT->ProcessLine(Form("a.SetOpt(\"bendSign\",%d);", bendSign));
   // Stop after this many vertices, for timing tests only.
   if (maxVertices > 0) gROOT->ProcessLine(Form("a.SetOpt(\"maxVertices\",%d);", maxVertices));
   gROOT->ProcessLine("a.SetOpt(\"progressEvery\",200);");

   gROOT->ProcessLine("a.Loop();");
}
