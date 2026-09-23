// genseed - Copyright 2026 Mikko Voutilainen - Apache License 2.0
//
// bench_tiledjet.C - is tiledjet.h right, and is it fast enough to be the
// default?  Exactness against a plain N^2 reference and against FastJet, and
// wall time per event, on real events of the pileup sample.
//
// WHY.  A jet finder that is wrong in one merge out of a million is wrong in
// one jet per event, and that jet is exactly the odd one the analysis will
// notice years later in a tail.  So the tiled, lazy, heap-driven clusterer is
// not trusted on its own: every jet above 1 GeV of 200 events - number of jets,
// pT to 1e-9 relative, and the SET of constituent indices - has to agree with
// tiledjet::cluster_n2, which is the same algorithm written the slow, obvious
// way with no tiles and no heap, and with FastJet (anti-kT, R = 0.4, E-scheme,
// Best strategy) when it is installed.  Three inputs are used because the
// analysis runs these regimes: the PF candidates with the stored PUPPI weight
// applied as pT*w and m*w and w <= 0 dropped (N ~ 130 in this neutrino-gun
// sample, where every interaction is pileup and PUPPI keeps little: what the
// analyzer clusters per vertex), all PF candidates unweighted (N ~ 2100, the
// size of a whole event), and all PF plus the generated particles of ALL
// interactions as ghosts at pT x 1e-9 (N ~ 10000; the ghost trick of v3/v4 that
// assigns every generated particle to the reco jet whose catchment area it
// falls in).  Ghosts stress the finder in the way it will be used: a
// heap where 1e18 sits next to 1e-2, tiles where 4/5 of the occupants are
// invisible, and merges that must not move a hard jet by more than rounding.
//
// WHAT 'IDENTICAL TO FASTJET' CAN MEAN.  The merge sequence is fixed by the
// ORDER of the d_ij, and two correct implementations disagree only where two
// d_ij tie.  Real PF input has ties: calorimeter candidates sit on the tower
// and crystal grids, so HF towers or ECAL photons share the same float eta or
// phi exactly, and a third particle mirrored between two of them is
// equidistant to the bit (an EXACT tie, which tiledjet and cluster_n2 both
// break towards the lower input index - so those two agree always) or to
// 1e-13 relative, the level at which two rapidity formulas of the same
// algebra round differently (a NEAR tie, decided by the last bit, and a
// different last bit is not a different algorithm).  The 'PF all, unweighted'
// input shows this in a few percent of the events (a handful of jets whose
// constituents are regrouped, occasionally one jet more or less above 1 GeV);
// the same candidates with the PUPPI weights, with the ghosts, or with pT,
// eta and phi jittered by a relative 1e-7 (the 'jittered' input, the same
// jitter to every backend) agree with FastJet jet by jet, which is the
// evidence that the differences are ties and not logic.  The tiledjet vs
// cluster_n2 pair tests the tiles, the lazy updates and the heap; the
// FastJet pair tests the algorithm.
//
// SPEED is measured the way the analysis pays for it: wall time (TStopwatch)
// per event including the constituent lists, for tiledjet, for cluster_n2 and
// for FastJet, as mean and rms over the events, per input size, and the ratio
// tiledjet/FastJet.  The verdict - every input within 1.0x of FastJet or not -
// is what the constant tiledjet::kAtLeastAsFastAsFastJet in tiledjet.h must
// say (it is set by hand from this output, with the numbers in kBenchmark).
// runGenSeed.C parses that constant from the text of tiledjet.h to choose the
// default backend: FastJet only if it is found and the verdict is false.
// (tiledjet::kPreferFastJet is derived from it and read by nothing.)
//
// USAGE (from the genseed directory, or with the path):
//   root -l -b -q 'bench_tiledjet.C(200, "../data/MC24NanoV15_PU_IT_OOT/NANOAODSIM_1.root")'
// Arguments: number of events, file, events for the N^2 reference (default
// all; it costs ~0.35 s per event at N ~ 10000), repetitions per event for the
// timing (the minimum over repetitions is used when > 1; default 1 = plain
// mean and rms).  The macro compiles ITSELF with ACLiC (-DBENCH_COMPILED,
// plus -DTILEDJET_USE_FASTJET and the include/link paths when FindFastJet
// locates an installation: $FASTJET, $CONDA_PREFIX, ~/miniforge3/envs/eejet,
// /usr/local, /opt/homebrew) into $TMPDIR/aclic_genseed, never into the
// source tree.
//
// The interpreted part is at the bottom; everything above it is compiled.
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TStopwatch.h>
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

#ifdef BENCH_COMPILED
#include "tiledjet.h"
#ifdef TILEDJET_USE_FASTJET
#include "fastjet/ClusterSequence.hh"
#endif

namespace {

  const int kMaxPF   = 16384;   // 2103 on average, 3878 at the top of NANOAODSIM_1.root
  const int kMaxGen  = 4096;    // largest single interaction there: 1121
  const int kMaxSlot = 200;     // the production wrote 120 slots
  const double kR = 0.4, kPtMin = 1.0, kGhost = 1e-9, kEtaMaxGhost = 10.0;

  // one event, the four inputs (see kInputName; the third is the second
  // with pT, eta and phi jittered by a relative 1e-7 - see the header)
  const int kNIn = 4;
  const char *kInputName[kNIn] = {"PF x PUPPI (w > 0)", "PF all, unweighted", "PF all, jittered 1e-7", "PF all + gen ghosts"};
  const double kJitter = 1e-7;
  struct Event {
    std::vector<tiledjet::PseudoJet> in[kNIn];
  };

  // a jet reduced to what is compared: pT and the sorted constituent set
  struct JetOut {
    double pt;
    std::vector<int> cons;
  };

  struct Stat {
    double sum = 0, sum2 = 0; int n = 0;
    void add(double x) { sum += x; sum2 += x*x; ++n; }
    double mean() const { return n ? sum/n : 0; }
    double rms()  const { if (n < 2) return 0; const double m = mean(), v = sum2/n - m*m; return v > 0 ? std::sqrt(v) : 0; }
  };

  std::vector<JetOut> Reduce(const std::vector<tiledjet::Jet> &js) {
    std::vector<JetOut> out(js.size());
    for (size_t k = 0; k < js.size(); ++k) {
      out[k].pt = js[k].p.pt();
      out[k].cons = js[k].constituents;
      std::sort(out[k].cons.begin(), out[k].cons.end());
    }
    return out;
  }

  // Returns the number of mismatching jets (or 1 if the jet counts differ)
  // and prints the first few.
  int Compare(const std::vector<JetOut> &a, const std::vector<JetOut> &b,
              const char *na, const char *nb, const char *input, int evt, int &printed) {
    if (a.size() != b.size()) {
      if (printed++ < 8) printf("  MISMATCH %s evt %d: %s has %zu jets, %s has %zu\n", input, evt, na, a.size(), nb, b.size());
      return 1;
    }
    int bad = 0;
    for (size_t k = 0; k < a.size(); ++k) {
      const double d = std::fabs(a[k].pt - b[k].pt), tol = 1e-9*std::max(a[k].pt, b[k].pt);
      const bool ptbad = d > tol, cbad = a[k].cons != b[k].cons;
      if (ptbad || cbad) {
        ++bad;
        if (printed++ < 8)
          printf("  MISMATCH %s evt %d jet %zu: pT %.9g vs %.9g (%s), constituents %zu vs %zu (%s)\n",
                 input, evt, k, a[k].pt, b[k].pt, ptbad ? "DIFFER" : "ok",
                 a[k].cons.size(), b[k].cons.size(), cbad ? "DIFFER" : "same");
      }
    }
    return bad;
  }

#ifdef TILEDJET_USE_FASTJET
  std::vector<JetOut> RunFastJet(const std::vector<fastjet::PseudoJet> &in, std::string *strategy) {
    static const fastjet::JetDefinition jd(fastjet::antikt_algorithm, kR, fastjet::E_scheme, fastjet::Best);
    fastjet::ClusterSequence cs(in, jd);
    std::vector<fastjet::PseudoJet> js = fastjet::sorted_by_pt(cs.inclusive_jets(kPtMin));
    std::vector<JetOut> out(js.size());
    for (size_t k = 0; k < js.size(); ++k) {
      out[k].pt = js[k].pt();
      const std::vector<fastjet::PseudoJet> cv = cs.constituents(js[k]);
      out[k].cons.reserve(cv.size());
      for (size_t m = 0; m < cv.size(); ++m) out[k].cons.push_back(cv[m].user_index());
      std::sort(out[k].cons.begin(), out[k].cons.end());
    }
    if (strategy) *strategy = cs.strategy_string();
    return out;
  }
#endif

  // Read the inputs of nevt events into memory: the timing must not
  // include the I/O of a 3.6 GB file.
  bool ReadEvents(const char *file, int nevt, std::vector<Event> &evs) {
    TFile *f = TFile::Open(file);
    if (!f || f->IsZombie()) { printf("bench: cannot open %s\n", file); return false; }
    TTree *t = (TTree*)f->Get("Events");
    if (!t) { printf("bench: no Events tree in %s\n", file); return false; }
    t->SetBranchStatus("*", 0);
    Int_t nPF = 0;
    static Float_t pf_pt[kMaxPF], pf_eta[kMaxPF], pf_phi[kMaxPF], pf_m[kMaxPF], pf_w[kMaxPF];
    static Int_t   ng[kMaxSlot];
    static Float_t g_pt[kMaxSlot][kMaxGen], g_eta[kMaxSlot][kMaxGen], g_phi[kMaxSlot][kMaxGen], g_m[kMaxSlot][kMaxGen];
    auto on = [&](const char *b, void *addr) {
      if (!t->GetBranch(b)) { printf("bench: missing branch %s\n", b); return false; }
      t->SetBranchStatus(b, 1); t->SetBranchAddress(b, addr); return true;
    };
    // The buffers are static and fixed, and GetEntry fills them before any
    // test on nPF or ng[s] can run: a test after it comes after the damage
    // (an oversized entry runs into the next slot's row or past the end).
    // So the file is checked up front, on the counter leaves' file-wide
    // maxima.  ROOT never reads more than that maximum into an array (the
    // leaf's ReadBasket clips to it), so this test is exact, not a guess.
    auto fits = [&](const char *counter, int limit) {
      const TLeaf *l = t->GetLeaf(counter);
      const int mx = l ? l->GetMaximum() : 0;
      if (mx <= limit) return true;
      printf("bench: %s reaches %d in %s, more than the %d the buffers hold - raise the limit in bench_tiledjet.C\n",
             counter, mx, file, limit);
      return false;
    };
    bool ok = on("nPFCandV2", &nPF) && fits("nPFCandV2", kMaxPF) && on("PFCandV2_pt", pf_pt) && on("PFCandV2_eta", pf_eta) &&
              on("PFCandV2_phi", pf_phi) && on("PFCandV2_mass", pf_m) && on("PFCandV2_puppiWeight", pf_w);
    if (!ok) return false;
    int nslot = 0;
    for (; nslot < kMaxSlot; ++nslot) {
      if (!t->GetBranch(Form("nGenPartCandBX0PUEvent%d", nslot))) break;
      if (!fits(Form("nGenPartCandBX0PUEvent%d", nslot), kMaxGen)) return false;
      on(Form("nGenPartCandBX0PUEvent%d", nslot), &ng[nslot]);
      on(Form("GenPartCandBX0PUEvent%d_pt", nslot),   g_pt[nslot]);
      on(Form("GenPartCandBX0PUEvent%d_eta", nslot),  g_eta[nslot]);
      on(Form("GenPartCandBX0PUEvent%d_phi", nslot),  g_phi[nslot]);
      on(Form("GenPartCandBX0PUEvent%d_mass", nslot), g_m[nslot]);
    }
    printf("bench: %s, %lld entries, %d gen slots; reading %d events\n", file, t->GetEntries(), nslot, nevt);
    const Long64_t n = std::min<Long64_t>(nevt, t->GetEntries());
    evs.clear(); evs.reserve(n);
    long nGhostEta = 0, nGhostPt0 = 0;
    unsigned long long rnd = 20260921ULL;
    for (Long64_t e = 0; e < n; ++e) {
      t->GetEntry(e);
      // (cannot fire after fits() above; kept so that a file whose counter
      // disagrees with its own leaf maximum never indexes past the buffers)
      if (nPF > kMaxPF) { printf("bench: entry %lld has %d PF candidates > %d, skipped\n", e, nPF, kMaxPF); continue; }
      Event ev;
      std::vector<tiledjet::PseudoJet> &pfw = ev.in[0], &pfall = ev.in[1], &pfjit = ev.in[2], &all = ev.in[3];
      pfw.reserve(nPF); pfall.reserve(nPF); pfjit.reserve(nPF);
      for (int i = 0; i < nPF; ++i) {
        const double w = pf_w[i];
        if (w > 0) pfw.push_back(tiledjet::PseudoJet::PtEtaPhiM(pf_pt[i]*w, pf_eta[i], pf_phi[i], pf_m[i]*w, (int)pfw.size()));
        if (pf_pt[i] > 0) {
          pfall.push_back(tiledjet::PseudoJet::PtEtaPhiM(pf_pt[i], pf_eta[i], pf_phi[i], pf_m[i], (int)pfall.size()));
          // deterministic jitter (a fixed-seed LCG), the same for every backend
          rnd = rnd*6364136223846793005ULL + 1442695040888963407ULL; const double j1 = (double)(rnd >> 11)/9007199254740992.0 - 0.5;
          rnd = rnd*6364136223846793005ULL + 1442695040888963407ULL; const double j2 = (double)(rnd >> 11)/9007199254740992.0 - 0.5;
          rnd = rnd*6364136223846793005ULL + 1442695040888963407ULL; const double j3 = (double)(rnd >> 11)/9007199254740992.0 - 0.5;
          pfjit.push_back(tiledjet::PseudoJet::PtEtaPhiM(pf_pt[i]*(1 + kJitter*j1), pf_eta[i]*(1 + kJitter*j2), pf_phi[i]*(1 + kJitter*j3), pf_m[i], (int)pfjit.size()));
        }
      }
      all = pfall;
      for (int s = 0; s < nslot; ++s) {
        // (the same kind of belt-and-braces test as for nPF)
        if (ng[s] > kMaxGen) { printf("bench: entry %lld slot %d has %d particles > %d, truncated\n", e, s, ng[s], kMaxGen); ng[s] = kMaxGen; }
        for (int g = 0; g < ng[s]; ++g) {
          if (!(g_pt[s][g] > 0)) { ++nGhostPt0; continue; }
          if (std::fabs(g_eta[s][g]) > kEtaMaxGhost) { ++nGhostEta; continue; }
          all.push_back(tiledjet::PseudoJet::PtEtaPhiM(g_pt[s][g]*kGhost, g_eta[s][g], g_phi[s][g], g_m[s][g]*kGhost, (int)all.size()));
        }
      }
      evs.push_back(std::move(ev));
    }
    printf("bench: %zu events in memory; ghosts dropped: %ld with pT <= 0, %ld with |eta| > %g\n",
           evs.size(), nGhostPt0, nGhostEta, kEtaMaxGhost);
    f->Close();
    return !evs.empty();
  }

} // namespace

void bench_tiledjet_run(int nevt, const char *file, int nN2, int nrep)
{
  std::vector<Event> evs;
  if (!ReadEvents(file, nevt, evs)) return;
  if (nN2 < 0 || nN2 > (int)evs.size()) nN2 = (int)evs.size();
  if (nrep < 1) nrep = 1;

  tiledjet::Config cfg; cfg.R = kR; cfg.p = -1; cfg.ptmin = kPtMin;
  const char **inputName = kInputName;
  Stat tN[kNIn], tTJ[kNIn], tN2[kNIn], tFJ[kNIn], nJ[kNIn];
  int badTJN2[kNIn] = {0}, badTJFJ[kNIn] = {0}, badN2FJ[kNIn] = {0}, evtTJFJ[kNIn] = {0}, printed = 0;
  const int kUnweighted = 1;   // the input where exact grid ties are expected
  std::string strategy[kNIn];
  const bool haveFJ =
#ifdef TILEDJET_USE_FASTJET
    true;
#else
    false;
#endif
  printf("bench: R = %g, p = %g (anti-kT), ptmin = %g GeV, %d events, N^2 reference on %d, %d repetition(s)%s\n",
         kR, cfg.p, kPtMin, (int)evs.size(), nN2, nrep, haveFJ ? ", FastJet available" : ", FastJet NOT available");

  // warm-up (page faults, library loading) on the first event, untimed
  { std::vector<tiledjet::Jet> j = tiledjet::cluster(evs[0].in[kNIn - 1], cfg); (void)j; }

  TStopwatch sw;
  for (size_t e = 0; e < evs.size(); ++e) {
    for (int in = 0; in < kNIn; ++in) {
      const std::vector<tiledjet::PseudoJet> &input = evs[e].in[in];
      tN[in].add((double)input.size());

      // tiledjet
      std::vector<JetOut> jTJ;
      double best = 1e30;
      for (int r = 0; r < nrep; ++r) {
        sw.Start(kTRUE);
        std::vector<tiledjet::Jet> js = tiledjet::cluster(input, cfg);
        sw.Stop();
        best = std::min(best, 1000.*sw.RealTime());
        if (r == 0) jTJ = Reduce(js);
      }
      tTJ[in].add(best); nJ[in].add((double)jTJ.size());

      // FastJet
#ifdef TILEDJET_USE_FASTJET
      std::vector<JetOut> jFJ;
      {
        std::vector<fastjet::PseudoJet> fin; fin.reserve(input.size());
        for (size_t i = 0; i < input.size(); ++i) {
          fin.push_back(fastjet::PseudoJet(input[i].px, input[i].py, input[i].pz, input[i].E));
          fin.back().set_user_index(input[i].user_index);
        }
        best = 1e30;
        for (int r = 0; r < nrep; ++r) {
          sw.Start(kTRUE);
          std::vector<JetOut> j = RunFastJet(fin, r == 0 ? &strategy[in] : 0);
          sw.Stop();
          best = std::min(best, 1000.*sw.RealTime());
          if (r == 0) jFJ = std::move(j);
        }
        tFJ[in].add(best);
        const int b = Compare(jTJ, jFJ, "tiledjet", "FastJet", inputName[in], (int)e, printed);
        badTJFJ[in] += b; if (b) ++evtTJFJ[in];
      }
#endif

      // N^2 reference
      if ((int)e < nN2) {
        sw.Start(kTRUE);
        std::vector<tiledjet::Jet> js = tiledjet::cluster_n2(input, cfg);
        sw.Stop();
        tN2[in].add(1000.*sw.RealTime());
        std::vector<JetOut> jN2 = Reduce(js);
        badTJN2[in] += Compare(jTJ, jN2, "tiledjet", "cluster_n2", inputName[in], (int)e, printed);
#ifdef TILEDJET_USE_FASTJET
        badN2FJ[in] += Compare(jN2, jFJ, "cluster_n2", "FastJet", inputName[in], (int)e, printed);
#endif
      }
    }
    if ((e + 1) % 50 == 0) printf("bench: %zu events done\n", e + 1);
  }

  // ---- report -------------------------------------------------------------
  printf("\nExactness (jets with pT >= %g GeV: count, pT to 1e-9 relative, constituent sets), %zu events:\n",
         kPtMin, evs.size());
  printf("  %-22s %14s %14s %14s %14s\n", "input", "tj vs n2", "tj vs FastJet", "n2 vs FastJet", "events tj!=FJ");
  for (int in = 0; in < kNIn; ++in)
    printf("  %-22s %14d %14s %14s %14s\n", inputName[in], badTJN2[in],
           haveFJ ? Form("%d", badTJFJ[in]) : "n/a", haveFJ ? Form("%d", badN2FJ[in]) : "n/a",
           haveFJ ? Form("%d", evtTJFJ[in]) : "n/a");
  printf("  (mismatching jets; a jet count that differs counts as one.  Differences against FastJet on\n"
         "   the unweighted input alone are the grid ties discussed in the header - the jittered copy\n"
         "   of the same candidates must show none.)\n");

  printf("\nTiming, ms/event (mean +- rms over events%s):\n", nrep > 1 ? Form(", min of %d repetitions per event", nrep) : "");
  printf("  %-22s %7s %7s %18s %18s %18s %8s  %s\n", "input", "<N>", "<jets>", "tiledjet", "cluster_n2", "FastJet", "tj/FJ", "FastJet strategy");
  bool within = haveFJ;
  for (int in = 0; in < kNIn; ++in) {
    const double ratio = tFJ[in].n && tFJ[in].mean() > 0 ? tTJ[in].mean()/tFJ[in].mean() : 0;
    if (haveFJ && !(ratio <= 1.0)) within = false;
    printf("  %-22s %7.0f %7.1f %8.3f +- %-7.3f %8.2f +- %-7.2f %8.3f +- %-7.3f %8s  %s\n",
           inputName[in], tN[in].mean(), nJ[in].mean(),
           tTJ[in].mean(), tTJ[in].rms(), tN2[in].mean(), tN2[in].rms(),
           tFJ[in].mean(), tFJ[in].rms(), haveFJ ? Form("%.2f", ratio) : "n/a",
           haveFJ ? strategy[in].c_str() : "-");
  }
  int badLogic = 0;
  for (int in = 0; in < kNIn; ++in) badLogic += badTJN2[in] + (in == kUnweighted ? 0 : badTJFJ[in]);
  printf("\nVerdict: %s; %d mismatching jets outside the grid-tie input%s; %d tie-resolution differences "
         "against FastJet in %d of %zu unweighted events.\n",
         !haveFJ ? "no FastJet to compare with" :
         within ? "tiledjet within 1.0x of FastJet on every input -> kAtLeastAsFastAsFastJet = true"
                : "tiledjet slower than FastJet on at least one input -> kAtLeastAsFastAsFastJet = false",
         badLogic, badLogic ? " - DO NOT USE until fixed" : "", badTJFJ[kUnweighted], evtTJFJ[kUnweighted], evs.size());
  printf("Recorded in tiledjet.h: kAtLeastAsFastAsFastJet = %s\n  kBenchmark = \"%s\"\n",
         tiledjet::kAtLeastAsFastAsFastJet ? "true" : "false", tiledjet::kBenchmark);
  printf("Machine: %s, %s\n", gSystem->GetFromPipe("uname -mrs").Data(), gSystem->GetBuildCompilerVersion());
}

#else  // ---- interpreted entry: compile this file and run ------------------

namespace {
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
      if (!gSystem->AccessPathName(tries[i] + "/include/fastjet/ClusterSequence.hh"))
        return tries[i];
    }
    return "";
  }
}

void bench_tiledjet(int nevt = 200,
                    const char *file = "../data/MC24NanoV15_PU_IT_OOT/NANOAODSIM_1.root",
                    int nN2 = -1, int nrep = 1)
{
  const TString dir = gSystem->DirName(__FILE__);
  const char *tmp = gSystem->Getenv("TMPDIR");
  const TString bdir = TString(tmp && tmp[0] ? tmp : "/tmp") + "/aclic_genseed";
  gSystem->mkdir(bdir, kTRUE);
  gSystem->SetBuildDir(bdir, kTRUE);

  TString opt = Form("-DBENCH_COMPILED -I%s", dir.Data());
  const TString fj = FindFastJet();
  if (fj.IsNull()) {
    printf("bench: no FastJet found ($FASTJET unset) - exactness against cluster_n2 only\n");
  } else {
    printf("bench: FastJet at %s\n", fj.Data());
    opt += Form(" -DTILEDJET_USE_FASTJET -I%s/include", fj.Data());
    gSystem->AddLinkedLibs(Form("-L%s/lib -Wl,-rpath,%s/lib -lfastjet", fj.Data(), fj.Data()));
    gSystem->AddDynamicPath(Form("%s/lib", fj.Data()));
  }
  gSystem->AddIncludePath(opt);
  if (!gSystem->CompileMacro(Form("%s/bench_tiledjet.C", dir.Data()), "k")) {
    printf("bench: compilation failed\n");
    return;
  }
  gROOT->ProcessLine(Form("bench_tiledjet_run(%d, \"%s\", %d, %d);", nevt, file, nN2, nrep));
}
#endif
