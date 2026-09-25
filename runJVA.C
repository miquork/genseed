// genseed - Copyright 2026 Mikko Voutilainen - Apache License 2.0
// runJVA.C - the ROOT driver of JVA (JVA.h, jvassoc.h): knobs, job split by
// entry range, clustering backend, build directory, exactly as runGenSeed.C
// does them.  Every numeric knob ends up in hist/hopts.
//
//   root -l -b -q 'runJVA.C("files.txt","_test",0,1,300)'             // 300 crossings, one job
//   root -l -b -q 'runJVA.C("files.txt","_v1",3,8)'                   // job 3 of 8
//   hadd -f rootfiles/JVA_v1.root rootfiles/JVA_v1_job*of8.root
//   ./runjva.sh v1 8                                                  // all of it
//   root -l -b -q 'runJVA.C("files.txt","_c2r1",0,1,300,"text/jec_v1.txt",4,12,1,0.5,1,2,0.5,1,10,0.35)'
//                                  // v2: candMode 2, recoilMode 1 (MHT above 10 GeV), sigmaJ 0.35
// `list` is a text file with one NanoAOD path per line (# comments), or a
// single .root file.  The output is rootfiles/JVA_<tag>.root (a leading '_'
// is added to the tag if it has none), with _job<k>of<n> for a split job.
//
// THE CORRECTION IS REQUIRED for the dijet part: the tag-probe balance and the
// JER are in pT^corr, and pT^raw at 5 GeV is a factor 1.5-2 below the link
// scale.  The default is the genseed pass-1 table text/jec_v1.txt, looked up
// from the working directory and then from the package directory; a table
// that is named but not found stops the job.  An empty jecFile runs without
// one (pT^corr = pT^raw), and says so: the MET, spike and assignment parts do
// not need it.
//
// JOBS, BUILD DIRECTORY, BACKEND: see runGenSeed.C.  Job k of n takes the
// even-aligned range [first, last) of the requested entries, so that halves
// a and b land in every job; an empty range writes a valid empty file; the
// library is built into $TMPDIR/aclic_genseed/<backend>, next to GenSeed's
// and never in the source tree (a sync client rewriting a mapped library is a
// bus error), and runjva.sh builds it once before launching jobs in parallel.
#include <TROOT.h>
#include <TChain.h>
#include <TSystem.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TString.h>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

namespace {
  // runGenSeed.C's FindFastJet and TiledjetVerdict, copied so that the two
  // drivers choose the same backend by the same rule.
  TString JvaFindFastJet()
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
  int JvaTiledjetVerdict(const char *header)
  {
    std::ifstream in(header);
    if (!in) return -1;
    const std::string name = "kAtLeastAsFastAsFastJet";
    std::string line;
    while (std::getline(in, line)) {
      line = line.substr(0, line.find("//"));
      const size_t k = line.find(name);
      if (k == std::string::npos) continue;
      size_t p = line.find_first_of("={", k + name.size());
      if (p == std::string::npos) continue;
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

void runJVA(const char *list = "files.txt",
            const char *tag  = "_test",
            int job = 0, int njobs = 1,
            double maxEvents = 0,
            const char *jecFile = "text/jec_v1.txt",
            double tauTrk = 4.0,
            double tauAll = 12.0,
            double sigmaK = 1.0,
            double trkMinPt = 0.5,
            int storeTuple = 1,
            int candMode = 0,
            double chMinPt = 0.5,
            int recoilMode = 0,
            double mhtMin = 10.,
            double sigmaJ = 0.35)
{
   gSystem->mkdir("rootfiles", kTRUE);
   const TString src = gSystem->DirName(__FILE__);      // the package directory, wherever it is run from

   // ---- input -------------------------------------------------------------
   std::vector<TString> all;
   if (TString(list).EndsWith(".root")) all.push_back(list);
   else {
      std::ifstream in(list); TString line;
      while (line.ReadLine(in)) { line = line.Strip(TString::kBoth);
         if (line.Length() && !line.BeginsWith("#")) all.push_back(line); }
   }
   if (all.empty()) { printf("runJVA: no input files in %s\n", list); return; }
   if (njobs < 1 || job < 0 || job >= njobs) { printf("runJVA: job %d of %d makes no sense\n", job, njobs); return; }
   TChain *c = new TChain("Events");
   for (size_t i = 0; i < all.size(); ++i) c->Add(all[i]);
   const Long64_t nall = c->GetEntries();
   printf("runJVA: %zu files, %lld entries, root %s\n", all.size(), nall,
          gSystem->Getenv("ROOTSYS") ? gSystem->Getenv("ROOTSYS") : "(ROOTSYS unset)");
   if (!nall) return;

   // ---- the correction ------------------------------------------------------
   TString jec = jecFile ? jecFile : "";
   if (!jec.IsNull() && gSystem->AccessPathName(jec)) {
      const TString alt = src + "/" + jec;
      if (!jec.BeginsWith("/") && !gSystem->AccessPathName(alt)) jec = alt;
      else { printf("runJVA: correction %s not found (also not in %s) - the dijet part needs it; stopping\n", jecFile, src.Data()); return; }
   }
   if (jec.IsNull()) printf("runJVA: WARNING no correction: pT^corr = pT^raw, the dijet results are uncorrected\n");

   // ---- the job's entry range, even-aligned --------------------------------
   const Long64_t nmax  = maxEvents > 0 ? std::min(nall, (Long64_t)maxEvents) : nall;
   const Long64_t first = 2*(Long64_t)floor(double(job)*nmax/(2.*njobs));
   const Long64_t last  = (job == njobs-1) ? nmax : 2*(Long64_t)floor(double(job+1)*nmax/(2.*njobs));
   printf("runJVA: job %d of %d, entries [%lld, %lld) of the %lld requested\n", job, njobs, first, last, nmax);
   if (njobs > 1 && nmax < 2*(Long64_t)njobs)
      printf("runJVA: WARNING %lld events for %d jobs, fewer than two per job: some ranges are empty%s\n",
             nmax, njobs, last <= first ? ", this one included" : "");

   // ---- the clustering backend and the build directory ---------------------
   const TString fj = JvaFindFastJet();
   const int verdict = JvaTiledjetVerdict(src + "/tiledjet.h");
   bool useFJ = false;
   if (verdict < 0 && !gSystem->AccessPathName(src + "/tiledjet.h"))
      printf("runJVA: tiledjet.h does not state kAtLeastAsFastAsFastJet; treating it as false\n");
   if (gSystem->AccessPathName(src + "/tiledjet.h")) {
      // jvassoc.h clusters the forward energy with tiledjet whatever the backend
      printf("runJVA: tiledjet.h is missing from %s and jvassoc.h needs it\n", src.Data()); return; }
   if (!fj.IsNull() && verdict != 1) useFJ = true;
   {
      const char *tmp = gSystem->Getenv("TMPDIR");
      const TString bdir = Form("%s/aclic_genseed/%s", tmp && tmp[0] ? tmp : "/tmp", useFJ ? "fastjet" : "tiledjet");
      gSystem->mkdir(bdir, kTRUE);
      gSystem->SetBuildDir(bdir, kTRUE);
      gSystem->AddIncludePath(Form("-I%s", src.Data()));
   }
   if (useFJ) {
      printf("runJVA: clustering with FastJet at %s (tiledjet verdict %d)\n", fj.Data(), verdict);
      gSystem->AddIncludePath(Form("-DGENSEED_USE_FASTJET -I%s/include", fj.Data()));
      gSystem->AddLinkedLibs(Form("-L%s/lib -Wl,-rpath,%s/lib -lfastjet", fj.Data(), fj.Data()));
      gSystem->AddDynamicPath(Form("%s/lib", fj.Data()));
   } else
      printf("runJVA: clustering with tiledjet.h (verdict %d, FastJet %s)\n", verdict,
             fj.IsNull() ? "not found" : "found but not used");

   TString t = tag ? tag : "";
   if (!t.BeginsWith("_")) t = "_" + t;
   const TString out = Form("rootfiles/JVA%s%s.root", t.Data(), njobs > 1 ? Form("_job%dof%d", job, njobs) : "");

   if (gROOT->LoadMacro(src + "/JVA.C+") < 0) { printf("runJVA: build failed\n"); return; }
   gROOT->ProcessLine(Form("JVA a((TTree*)%p,\"%s\");", (void*)c, out.Data()));

   // ---- knobs ---------------------------------------------------------------
   gROOT->ProcessLine(Form("a.SetOpt(\"firstEntry\",%lld);", first));
   gROOT->ProcessLine(Form("a.SetOpt(\"lastEntry\",%lld);", last));
   if (!jec.IsNull()) gROOT->ProcessLine(Form("a.SetJEC(\"%s\");", jec.Data()));
   // One radius for the forward clusters, their track cone and every jet.
   gROOT->ProcessLine("a.SetOpt(\"R\",0.4);");
   // The per-vertex jets are clustered from 1 GeV raw, as in genseed.
   gROOT->ProcessLine("a.SetOpt(\"ptMinCluster\",1.0);");
   // Past this |eta| PUPPI is vertex-blind (puppi::Config::etaTracker).
   gROOT->ProcessLine("a.SetOpt(\"etaFwd\",2.5);");
   // The candidates of a forward cluster (jvassoc.h, THE CANDIDATES): 0 the
   // vertices with trkMinPt of tracker tracks in its cone (v1), 1 those with
   // chMinPt of its own charged constituents (their vertexRef, which PUPPI
   // ignores past the tracker), 2 either.
   gROOT->ProcessLine(Form("a.SetOpt(\"candMode\",%d);", candMode));
   gROOT->ProcessLine(Form("a.SetOpt(\"trkMinPt\",%g);", trkMinPt));
   gROOT->ProcessLine(Form("a.SetOpt(\"chMinPt\",%g);", chMinPt));
   // The recoil a cluster balances (jvassoc.h, THE RECOIL): 0 the particle MET
   // of the vertex-resolved part, sigma^2 = sigma0^2 + sigmaK^2 S_v, measured
   // 1.0 and 1.0 (sqrt(1 + S) GeV per component; v1); 1 the MHT of the
   // vertex-resolved jets above mhtMin (raw), sigma^2 = sigma0^2 + sigmaJ^2
   // sum pT_j^2 + sigmaK^2 (S_v - HT_v).
   gROOT->ProcessLine(Form("a.SetOpt(\"recoilMode\",%d);", recoilMode));
   gROOT->ProcessLine(Form("a.SetOpt(\"mhtMin\",%g);", mhtMin));
   gROOT->ProcessLine("a.SetOpt(\"sigma0\",1.0);");
   gROOT->ProcessLine(Form("a.SetOpt(\"sigmaK\",%g);", sigmaK));
   gROOT->ProcessLine(Form("a.SetOpt(\"sigmaJ\",%g);", sigmaJ));
   // The price of assigning a cluster: to one of its track candidates, and to
   // any vertex when it has none.
   gROOT->ProcessLine(Form("a.SetOpt(\"tauTrk\",%g);", tauTrk));
   gROOT->ProcessLine(Form("a.SetOpt(\"tauAll\",%g);", tauAll));
   // Exhaustive search per component up to this many combinations, else descent.
   gROOT->ProcessLine("a.SetOpt(\"maxCombos\",200000);");
   gROOT->ProcessLine("a.SetOpt(\"maxSweeps\",50);");
   // Dijets: jets above this pT^corr; type-1 MET from jets above t1Min; back to back.
   gROOT->ProcessLine("a.SetOpt(\"ptMinJet\",3.0);");
   gROOT->ProcessLine("a.SetOpt(\"t1Min\",3.0);");
   gROOT->ProcessLine("a.SetOpt(\"dphiMin\",2.7);");
   // genseed's ownership (0.05 cm) and the oracle's vertex window (0.2 cm).
   gROOT->ProcessLine("a.SetOpt(\"dzVertexGen\",0.05);");
   gROOT->ProcessLine("a.SetOpt(\"dzOracle\",0.2);");
   gROOT->ProcessLine("a.SetOpt(\"minTrkPerVertex\",3);");
   gROOT->ProcessLine(Form("a.SetOpt(\"storeTuple\",%d);", storeTuple));
   gROOT->ProcessLine("a.SetOpt(\"progressEvery\",200);");
   // Any other knob, or an override of one above: JVAOPT="key=value,key=value"
   // in the environment (runjva.sh passes it on to every job), e.g.
   //   JVAOPT="etaFwd=3,pupEtaTracker=3,pupEtaVtxAssoc=3" ./runjva.sh vtx30 8
   if (const char *jo = gSystem->Getenv("JVAOPT")) {
     TObjArray *kv = TString(jo).Tokenize(",");
     for (int i = 0; i < kv->GetEntries(); ++i) {
       TString item = ((TObjString*)kv->At(i))->GetString();
       const int eq = item.Index("=");
       if (eq <= 0) { printf("runJVA: JVAOPT item '%s' is not key=value, ignored\n", item.Data()); continue; }
       const TString key = item(0, eq), val = item(eq + 1, item.Length());
       printf("runJVA: JVAOPT %s = %s\n", key.Data(), val.Data());
       gROOT->ProcessLine(Form("a.SetOpt(\"%s\",%s);", key.Data(), val.Data()));
     }
     delete kv;
   }

   gROOT->ProcessLine("a.Loop();");
}
