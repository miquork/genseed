// genseed - Copyright 2026 Mikko Voutilainen - Apache License 2.0
// scanJVA.C - the forward-cluster association of JVA re-run offline from the
// tuple that JVA.C writes (tree "jva" in rootfiles/JVA_<tag>.root), with the
// analyzer's own optimiser: jva::Assign of jvassoc.h on jva::Vertex and
// jva::Cluster filled the way JVA.C fills them.  Nothing here re-implements
// the objective or the search.
//
//   root -l -b -q -e 'gSystem->SetBuildDir("'$TMPDIR'/aclic_scanjva",true);' \
//        'scanJVA.C+("rootfiles/JVA_v2.root")' | tee rootfiles/scanjva_v2.log
//   optional arguments: threads (8), maxEntries (0 = all)
//
// Build it outside the source tree (the SetBuildDir above): ACLiC otherwise
// drops its .so/.d/.pcm next to the macro, in a synced directory.
//
// WHAT IT PRINTS.
// (0) The check that the offline rebuild is the analyzer: the file's own knobs
//     (hist/hopts / jobs) applied to the tuple give back its candidates
//     (cl_nopt), its jva hosts (cl_host_jva) and its trk choices
//     (cl_host_trk), mismatches counted.  The tuple is float, the analyzer
//     double; a near-tie can flip, and the count says how often it did.
// (a) The recoil.  jva moves cluster k to vertex v when |R_v - p_k|^2 beats
//     |R_v|^2 by the price, so R_v must predict what v's forward clusters
//     carry.  Sign convention of the tuple: M_v = vtx_mx/my = -(sum of the
//     vertex-resolved weighted p), MHT_T = vtx_mht<T> = -(sum of the
//     vertex-resolved jets above T in pT^raw), both "what is missing";
//     genfwd = +(generated pT of the owner at etaFwd < |eta| < 5), what its
//     forward clusters carry, so a perfect recoil is R_v = genfwd + genmet
//     with genmet = MET_true(owner) = -(gencen + genfwd), and R_v ~ genfwd
//     up to the owner's true imbalance.  For every owned vertex and both
//     components (pooled; the azimuthal symmetry makes the means 0, so the
//     fits go through the origin): the slope b of R on the truth, the
//     correlation r, the residual rms of R - b*truth, and the rms of truth -
//     (sum R.truth/sum R^2) R, the resolution on the truth that R gives;
//     overall and in bins of |genfwd| and of S_v.  Then the chain that says
//     where the prediction is lost: the owner's central generated recoil
//     -gencen and the truth jets -genmht_T against genfwd (what an ideal
//     detector would have), MHT_T against -genmht_T' (the jet response and
//     the threshold it effectively applies), and F_v = sum of p_k over the
//     clusters whose true vertex is v (what the ojet oracle gives v, on the
//     reco scale: what jva actually has to find).
// (b) The candidates.  Single-interaction clusters (cl_true >= 0), in bins of
//     their pT: is the true vertex among the candidates, and does it carry the
//     most evidence (the analyzer's trk rule: T_ku, C_ku, T_ku + C_ku, the
//     lowest index on a tie), for candMode 0, 1, 2, and the argmax of the
//     charged evidence over every u with C_ku > 0.
// (c) The grid: recoil {MET, MHT3, MHT5, MHT7, MHT10, MHT15} x candMode
//     {0, 1, 2} x (tauTrk, tauAll) {(4,12), (2,8), (0,6), (0,1e9)} x sigma
//     scale {1, 1.5} (sigma0, sigmaK and sigmaJ all scaled, sigma_v^2 by
//     2.25), and the references none (every cluster null), ojet (the
//     oracle), the stored jva, trk and lv, and tracks alone without the MET:
//     the trk rule of candMode 0, 1, 2 with null instead of the PV for a
//     cluster without candidates.  For clusters above 5 GeV the analyzer's
//     outcome categories; net right = right - wrong, with its binomial error
//     sqrt(R + W - (R - W)^2/N_single) and, as a check, the error from the
//     spread across crossings; right - kept_combo, the change of "correct"
//     (right + nulled_combo, the implementer's measure) against none; the
//     MET: MET_v = M_v - sum of the p_k given to v (the particle MET,
//     whatever the recoil), its rms per component over every vertex and at
//     the PV (the circular figure: jva minimises it), and its residual
//     against MET_true(owner) over the owned vertices, with the paired
//     difference to none, d<res^2> per component, and its error from the
//     spread of the per-crossing sums; and d|F - A|^2, the sum over the
//     vertices of |F_v - A_v|^2 - |F_v|^2 per crossing, A_v the forward pT
//     the method gave v and F_v what the oracle gives it: an energy-weighted
//     net right that does not depend on the MET at all (negative when what is
//     handed out is more the right vertex's energy than noise; ojet reaches
//     -sum |F_v|^2).
//     The residual against MET_true is NOT the non-circular measure it looks
//     like: MET_true is 3.2 GeV rms per component against ~7 GeV of MET noise,
//     so whatever lowers |MET_v| lowers the residual about as much: on the v2
//     production d<res^2> is d<MET^2> to ~10% at every grid point and for the
//     tracks-only rules, and only the oracle lowers the residual twice as much
//     as the MET.  Almost every grid point passes "no worse than none"
//     whatever it does to the clusters; net right, right - kept_combo and
//     d|F - A|^2 are the ones that discriminate.
// (d) The configuration with the largest net right among those whose residual
//     does not exceed none's, and what it gains over none; the largest net
//     right per recoil and the ten largest overall.
//
// THE TUPLE AS INPUT.  Candidates in candMode 1 and 2 are rebuilt from the
// in-cluster charged evidence ch_* (every u with C_ku > 0, so any chMinPt
// works) and the in-cone track evidence cand_* (only u with T_ku >= the
// production's trkMinPt, so the grid takes trkMinPt = chMinPt = 0.5 and needs
// a production at trkMinPt <= 0.5).  The one thing the tuple cannot give back
// is a T_ku below trkMinPt of a candidate that C_ku made: it enters the
// candMode-2 trk ranking (T + C) in the analyzer and is 0 here, which can move
// a trk choice (counted in (0) when the production is candMode 2); jva's
// assignment does not use the evidence, only the candidate list.
#include "jvassoc.h"
#include <TFile.h>
#include <TH1.h>
#include <TROOT.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <atomic>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <memory>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>
#include <algorithm>

namespace scanjva {

  const int kNT = 5;
  const int kT[kNT] = {3, 5, 7, 10, 15};             // the MHT thresholds of the tuple
  const int kNRec = 1 + kNT;                          // recoil: 0 particle MET, 1 + t MHT above kT[t]
  const int kNMode = 3;                               // candMode 0 1 2
  struct PriceSet { double tt, ta; const char *name; };
  const PriceSet kPrice[] = {{4,12,"4/12"}, {2,8,"2/8"}, {0,6,"0/6"}, {0,1e9,"0/1e9"}};
  const int kNPrice = 4;
  const double kSigScale[] = {1., 1.5};
  const int kNSig = 2;
  const int kNCfg = kNRec*kNMode*kNPrice*kNSig;       // 144
  inline int Cfg(int r, int m, int p, int s) { return ((r*kNMode + m)*kNPrice + p)*kNSig + s; }
  inline std::string RecName(int r) { return r == 0 ? std::string("MET") : "MHT" + std::to_string(kT[r-1]); }
  const double kTrkMin = 0.5, kChMin = 0.5;           // the grid's candidate thresholds (the analyzer's defaults)

  // references: the stored hosts, and the trk rule of candMode 0 1 2 with null
  // instead of the PV for a cluster without candidates (tracks alone, no MET)
  enum { kRefNone = 0, kRefOjet, kRefJva, kRefTrk, kRefLv, kRefT0, kRefT1, kRefT2, kNRef };
  const char *kRefName[kNRef] = {"none", "ojet", "jva (stored)", "trk (stored)", "lv (stored)",
                                 "trk c0, else null", "trk c1, else null", "trk c2, else null"};
  enum { kRight = 0, kWrong, kNulledSingle, kNulledCombo, kKeptCombo, kNOut };
  const double kClMin = 5.;                           // GeV, the clusters of the outcome counts (the analyzer's summary)

  // (a): predictors 0 M_v, 1..5 MHT_T, 6..10 -genmht_T, 11 -gencen, 12 F_v;
  // targets 0 genfwd, 1 genmet, 2 -gencen, 3..7 -genmht_T, 8 F_v;
  // categories 0 all owned, 1..4 |genfwd| bins, 5..9 S_v bins
  const int kNP = 13, kNG = 9;
  const double kGf[] = {0, 5, 10, 20, 1e30};  const int kNGf = 4;
  const double kSv[] = {0, 10, 20, 40, 80, 1e30}; const int kNSv = 5;
  const int kNCatA = 1 + kNGf + kNSv;
  inline int Bin(double x, const double *e, int n) { for (int i = 0; i < n; ++i) if (x < e[i+1]) return i; return n - 1; }
  const char *PName(int p) {
    static const char *s[kNP] = {"M_v (particle MET)","MHT3","MHT5","MHT7","MHT10","MHT15",
      "-genmht3","-genmht5","-genmht7","-genmht10","-genmht15","-gencen","F_v (oracle clusters)"};
    return s[p]; }
  const char *GName(int g) {
    static const char *s[kNG] = {"genfwd","genmet","-gencen","-genmht3","-genmht5","-genmht7","-genmht10","-genmht15","F_v"};
    return s[g]; }

  // (b): cluster pT bins, rules 0..2 candMode, 3 the charged argmax (C_ku > 0)
  const double kPtB[] = {2, 5, 10, 20, 1e30}; const int kNPtB = 4;
  const int kNRule = 4;

  // Accumulators: plain doubles only, so that the threads' copies add up as arrays.
  struct Reg  { double pg[kNCatA][kNP][kNG], pp[kNCatA][kNP], gg[kNCatA][kNG], n[kNCatA], nz[kNT], nown; };
  struct Cand { double n[kNPtB], has[kNRule][kNPtB], among[kNRule][kNPtB], leads[kNRule][kNPtB]; };
  struct Grid {
    double out[2][kNOut];                             // above kClMin: [0] with candidates under the rule, [1] without
    double outPt[3][kNOut];                           // 5-10, 10-20, > 20 GeV
    double res2, resN, res2pv, resNpv, met2, metN, met2pv, metNpv;
    double dSum, dSum2, nrSum, nrSum2, nCross;        // per crossing: d(res^2) to none, right - wrong
    double qSum, qSum2;                               // per crossing: sum_v |F_v - A_v|^2 - |F_v|^2 (A_v = what v got)
  };
  template <class T> void AddTo(T &a, const T &b) {
    static_assert(std::is_trivially_copyable<T>::value && sizeof(T) % sizeof(double) == 0, "doubles only");
    double *x = reinterpret_cast<double*>(&a); const double *y = reinterpret_cast<const double*>(&b);
    for (size_t i = 0; i < sizeof(T)/sizeof(double); ++i) x[i] += y[i];
  }
  template <class T> void Zero(T &a) { memset(&a, 0, sizeof(T)); }

  // The production's knobs, from hist/hopts divided by the jobs.
  struct Knobs { bool ok = false; jva::Config c; int njobs = 0; };
  inline Knobs FileKnobs(TFile *f) {
    Knobs k;
    TH1 *o = (TH1*)f->Get("hist/hopts"), *h = (TH1*)f->Get("hist/hcount");
    if (!o || !h || !(h->GetBinContent(6) > 0)) return k;
    const double nj = h->GetBinContent(6);
    auto get = [&](const char *lab, double def) {
      for (int i = 1; i <= o->GetNbinsX(); ++i)
        if (!strcmp(o->GetXaxis()->GetBinLabel(i), lab)) return std::round(o->GetBinContent(i)/nj*1e6)/1e6;
      return def; };
    jva::Config &c = k.c;
    c.R = get("R", 0.4); c.etaFwd = get("etaFwd", 2.5);
    c.candMode = int(get("candMode", 0)); c.trkMinPt = get("trkMinPt", 0.5); c.chMinPt = get("chMinPt", 0.5);
    c.sigma0 = get("sigma0", 1.); c.sigmaK = get("sigmaK", 1.); c.sigmaJ = get("sigmaJ", 0.35);
    c.recoilMode = int(get("recoilMode", 0)); c.mhtMin = get("mhtMin", 10.);
    c.tauTrk = get("tauTrk", 4.); c.tauAll = get("tauAll", 12.);
    c.maxCombos = get("maxCombos", 200000); c.maxSweeps = int(get("maxSweeps", 50));
    k.njobs = int(nj); k.ok = true;
    return k;
  }

  // One crossing of the tuple.
  struct Raw {
    int nv = 0, pv = -1, nk = 0;
    std::vector<int> owner;
    std::vector<double> mx, my, s, gmx, gmy, gfx, gfy, gcx, gcy;
    std::vector<double> ht[kNT], sj2[kNT], hx[kNT], hy[kNT], ghx[kNT], ghy[kNT];
    std::vector<double> px, py, pt, eta, phi;
    std::vector<int> nopt, tru, host[kNRef];
    std::vector< std::vector<int> > tv, cv; std::vector< std::vector<double> > tp, cp;
  };

  // Cluster k as jva::BuildClusters makes it, from the evidence in the tuple:
  // u ascending, T_ku and C_ku of every candidate, the candMode rule.
  inline void MakeClusters(const Raw &x, int mode, double trkMin, double chMin, std::vector<jva::Cluster> &out) {
    out.assign(x.nk, jva::Cluster());
    std::vector<double> t(x.nv), h(x.nv);
    for (int k = 0; k < x.nk; ++k) {
      jva::Cluster &q = out[k];
      q.px = x.px[k]; q.py = x.py[k]; q.pt = x.pt[k]; q.eta = x.eta[k]; q.phi = x.phi[k];
      std::fill(t.begin(), t.end(), 0.); std::fill(h.begin(), h.end(), 0.);
      for (size_t j = 0; j < x.tv[k].size(); ++j) if (x.tv[k][j] >= 0 && x.tv[k][j] < x.nv) t[x.tv[k][j]] = x.tp[k][j];
      for (size_t j = 0; j < x.cv[k].size(); ++j) if (x.cv[k][j] >= 0 && x.cv[k][j] < x.nv) h[x.cv[k][j]] = x.cp[k][j];
      for (int u = 0; u < x.nv; ++u) {
        const bool bt = t[u] >= trkMin && t[u] > 0, bh = h[u] >= chMin && h[u] > 0;
        if (bt) { q.trkVtx.push_back(u); q.trkPt.push_back(t[u]); }
        if (h[u] > 0) { q.chVtx.push_back(u); q.chPt.push_back(h[u]); }
        if (mode == 1 ? bh : (mode == 2 ? (bt || bh) : bt)) {
          q.cand.push_back(u); q.candTrkPt.push_back(t[u]); q.candChPt.push_back(h[u]); }
      }
    }
  }
  // the analyzer's trk rule: the candidate with the most evidence (strict >,
  // ascending, so the lowest index on a tie), else the PV
  inline int TrkChoice(const jva::Cluster &q, int mode, int pv) {
    jva::Config c; c.candMode = mode;
    int t = pv; double best = -1;
    for (size_t j = 0; j < q.cand.size(); ++j) { const double e = jva::Evidence(q, j, c); if (e > best) { best = e; t = q.cand[j]; } }
    return t;
  }
  inline int ArgMax(const std::vector<int> &v, const std::vector<double> &e) {
    int a = -1; double b = -1; for (size_t j = 0; j < v.size(); ++j) if (e[j] > b) { b = e[j]; a = v[j]; } return a; }

  // One assignment's outcomes and MET into g; rn2 = |M_v - MET_true|^2 of none per vertex,
  // (fx, fy) = F_v, the clusters whose true vertex is v (what ojet gives v).
  inline void Tally(Grid &g, const Raw &x, const std::vector<int> &h, const std::vector<char> &hasCand,
                    const std::vector<double> &rn2, const std::vector<double> &fx, const std::vector<double> &fy,
                    std::vector<double> &mx, std::vector<double> &my) {
    double nr = 0;
    for (int k = 0; k < x.nk; ++k) {
      if (!(x.pt[k] >= kClMin)) continue;
      const int tv = x.tru[k], hk = h[k];
      const int o = tv >= 0 ? (hk == tv ? kRight : (hk < 0 ? kNulledSingle : kWrong)) : (hk < 0 ? kNulledCombo : kKeptCombo);
      g.out[hasCand[k] ? 0 : 1][o] += 1;
      g.outPt[x.pt[k] < 10 ? 0 : (x.pt[k] < 20 ? 1 : 2)][o] += 1;
      if (o == kRight) nr += 1; else if (o == kWrong) nr -= 1;
    }
    mx = x.mx; my = x.my;
    for (int k = 0; k < x.nk; ++k) if (h[k] >= 0 && h[k] < x.nv) { mx[h[k]] -= x.px[k]; my[h[k]] -= x.py[k]; }
    double d = 0, q = 0;
    for (int v = 0; v < x.nv; ++v) {
      const double ax = x.mx[v] - mx[v], ay = x.my[v] - my[v];     // A_v
      q += ax*ax + ay*ay - 2*(ax*fx[v] + ay*fy[v]);
      const double m2 = mx[v]*mx[v] + my[v]*my[v];
      g.met2 += m2; g.metN += 2;
      if (v == x.pv) { g.met2pv += m2; g.metNpv += 2; }
      if (x.owner[v] < 0) continue;
      const double dx = mx[v] - x.gmx[v], dy = my[v] - x.gmy[v], r2 = dx*dx + dy*dy;
      g.res2 += r2; g.resN += 2;
      if (v == x.pv) { g.res2pv += r2; g.resNpv += 2; }
      d += r2 - rn2[v];
    }
    g.dSum += d; g.dSum2 += d*d; g.nrSum += nr; g.nrSum2 += nr*nr; g.nCross += 1; g.qSum += q; g.qSum2 += q*q;
  }

  struct Worker {
    std::string file; Long64_t first = 0, last = 0; Knobs fk;
    Reg reg; Cand cand; Grid grid[kNCfg]; Grid ref[kNRef];
    double nCross = 0, nCl = 0, optMis = 0, hostMis = 0, hostMisCross = 0, trkMis = 0, validated = 0;
    std::string err;
    std::atomic<long> *progress = 0;
    Worker() { Zero(reg); Zero(cand); for (int i = 0; i < kNCfg; ++i) Zero(grid[i]); for (int i = 0; i < kNRef; ++i) Zero(ref[i]); }
    void Run();
  };

  void Worker::Run() {
    std::unique_ptr<TFile> f(TFile::Open(file.c_str()));
    if (!f || f->IsZombie()) { err = "cannot open " + file; return; }
    TTreeReader r("jva", f.get());
    if (r.SetEntriesRange(first, last) != TTreeReader::kEntryValid && last > first) { err = "bad entry range"; return; }
    typedef TTreeReaderArray<float> AF; typedef TTreeReaderArray<int> AI;
    AI key(r,"vtx_key"), own(r,"vtx_owner");
    AF mx(r,"vtx_mx"), my(r,"vtx_my"), s(r,"vtx_s");
    AF gmx(r,"vtx_genmet_x"), gmy(r,"vtx_genmet_y"), gfx(r,"vtx_genfwd_x"), gfy(r,"vtx_genfwd_y"), gcx(r,"vtx_gencen_x"), gcy(r,"vtx_gencen_y");
    std::vector< std::unique_ptr<AF> > ht, sj2, hx, hy, ghx, ghy;
    for (int t = 0; t < kNT; ++t) {                   // (no Form() in a thread)
      const std::string T = std::to_string(kT[t]);
      ht.emplace_back(new AF(r, ("vtx_ht" + T).c_str()));        sj2.emplace_back(new AF(r, ("vtx_htsq" + T).c_str()));
      hx.emplace_back(new AF(r, ("vtx_mht" + T + "_x").c_str()));  hy.emplace_back(new AF(r, ("vtx_mht" + T + "_y").c_str()));
      ghx.emplace_back(new AF(r, ("vtx_genmht" + T + "_x").c_str())); ghy.emplace_back(new AF(r, ("vtx_genmht" + T + "_y").c_str()));
    }
    AF cpt(r,"cl_pt"), ceta(r,"cl_eta"), cphi(r,"cl_phi"), cpx(r,"cl_px"), cpy(r,"cl_py");
    AI nopt(r,"cl_nopt"), ncand(r,"cl_ncand"), coff(r,"cl_candoff"), cvtx(r,"cand_vtx"), nch(r,"cl_nch"), choff(r,"cl_choff"), chv(r,"ch_vtx");
    AF ctp(r,"cand_trkpt"), chp(r,"ch_pt");
    AI tru(r,"cl_true"), hoj(r,"cl_host_ojet"), hjv(r,"cl_host_jva"), htk(r,"cl_host_trk"), hlv(r,"cl_host_lv");

    Raw x; std::vector<jva::Cluster> CL[kNMode], CLf; std::vector<jva::Vertex> V[kNRec];
    std::vector<char> has[kNMode], hasF; std::vector<double> rn2, sx, sy, fvx, fvy;
    // the file's recoil among the tuple's (validation), -1 if it is not one of them
    int recF = -1;
    if (fk.ok) {
      if (fk.c.recoilMode == 0) recF = 0;
      else for (int t = 0; t < kNT; ++t) if (std::fabs(fk.c.mhtMin - kT[t]) < 1e-9) recF = 1 + t;
    }
    while (r.Next()) {
      if (r.GetEntryStatus() != TTreeReader::kEntryValid) { err = "read error"; return; }
      // ---- the crossing ------------------------------------------------------------
      const int nv = (int)key.GetSize(), nk = (int)cpt.GetSize();
      x.nv = nv; x.nk = nk; x.pv = -1;
      x.owner.assign(own.begin(), own.end());
      x.mx.assign(mx.begin(), mx.end()); x.my.assign(my.begin(), my.end()); x.s.assign(s.begin(), s.end());
      x.gmx.assign(gmx.begin(), gmx.end()); x.gmy.assign(gmy.begin(), gmy.end());
      x.gfx.assign(gfx.begin(), gfx.end()); x.gfy.assign(gfy.begin(), gfy.end());
      x.gcx.assign(gcx.begin(), gcx.end()); x.gcy.assign(gcy.begin(), gcy.end());
      for (int t = 0; t < kNT; ++t) {
        x.ht[t].assign(ht[t]->begin(), ht[t]->end()); x.sj2[t].assign(sj2[t]->begin(), sj2[t]->end());
        x.hx[t].assign(hx[t]->begin(), hx[t]->end()); x.hy[t].assign(hy[t]->begin(), hy[t]->end());
        x.ghx[t].assign(ghx[t]->begin(), ghx[t]->end()); x.ghy[t].assign(ghy[t]->begin(), ghy[t]->end());
      }
      for (int v = 0; v < nv; ++v) if (key[v] == 0) x.pv = v;
      x.px.assign(cpx.begin(), cpx.end()); x.py.assign(cpy.begin(), cpy.end()); x.pt.assign(cpt.begin(), cpt.end());
      x.eta.assign(ceta.begin(), ceta.end()); x.phi.assign(cphi.begin(), cphi.end());
      x.nopt.assign(nopt.begin(), nopt.end()); x.tru.assign(tru.begin(), tru.end());
      x.host[kRefNone].assign(nk, jva::kNull);
      x.host[kRefOjet].assign(hoj.begin(), hoj.end()); x.host[kRefJva].assign(hjv.begin(), hjv.end());
      x.host[kRefTrk].assign(htk.begin(), htk.end()); x.host[kRefLv].assign(hlv.begin(), hlv.end());
      x.tv.resize(nk); x.tp.resize(nk); x.cv.resize(nk); x.cp.resize(nk);
      for (int k = 0; k < nk; ++k) {
        x.tv[k].clear(); x.tp[k].clear(); x.cv[k].clear(); x.cp[k].clear();
        for (int j = 0; j < ncand[k]; ++j) { x.tv[k].push_back(cvtx[coff[k] + j]); x.tp[k].push_back(ctp[coff[k] + j]); }
        for (int j = 0; j < nch[k]; ++j)   { x.cv[k].push_back(chv[choff[k] + j]);  x.cp[k].push_back(chp[choff[k] + j]); }
      }
      nCross += 1; nCl += nk;

      // the recoils and the candidates the grid uses
      for (int rr = 0; rr < kNRec; ++rr) {
        V[rr].assign(nv, jva::Vertex());
        for (int v = 0; v < nv; ++v) {
          jva::Vertex &q = V[rr][v];
          q.s = x.s[v];
          if (rr == 0) { q.mx = x.mx[v]; q.my = x.my[v]; }
          else { const int t = rr - 1; q.mx = x.hx[t][v]; q.my = x.hy[t][v]; q.ht = x.ht[t][v]; q.sj2 = x.sj2[t][v]; }
        }
      }
      for (int m = 0; m < kNMode; ++m) {
        MakeClusters(x, m, kTrkMin, kChMin, CL[m]);
        has[m].resize(nk); for (int k = 0; k < nk; ++k) has[m][k] = !CL[m][k].cand.empty();
      }

      // ---- (0) the production's own knobs --------------------------------------------
      if (fk.ok) {
        MakeClusters(x, fk.c.candMode, fk.c.trkMinPt, fk.c.chMinPt, CLf);
        hasF.resize(nk);
        for (int k = 0; k < nk; ++k) {
          hasF[k] = x.nopt[k] > 0;
          if ((int)CLf[k].cand.size() != x.nopt[k]) optMis += 1;
          if (TrkChoice(CLf[k], fk.c.candMode, x.pv) != x.host[kRefTrk][k]) trkMis += 1;
        }
        if (recF >= 0) {
          const std::vector<int> h = jva::Assign(V[recF], CLf, fk.c);
          int mis = 0; for (int k = 0; k < nk; ++k) if (h[k] != x.host[kRefJva][k]) ++mis;
          hostMis += mis; if (mis) hostMisCross += 1; validated += 1;
        }
      } else { hasF.assign(nk, 0); }

      // ---- (a) the recoil against the truth, owned vertices ---------------------------
      fvx.assign(nv, 0.); fvy.assign(nv, 0.);
      for (int k = 0; k < nk; ++k) if (x.tru[k] >= 0 && x.tru[k] < nv) { fvx[x.tru[k]] += x.px[k]; fvy[x.tru[k]] += x.py[k]; }
      for (int v = 0; v < nv; ++v) {
        if (x.owner[v] < 0) continue;
        reg.nown += 1;
        for (int t = 0; t < kNT; ++t) if (x.ht[t][v] > 0) reg.nz[t] += 1;
        const double agf = std::sqrt(x.gfx[v]*x.gfx[v] + x.gfy[v]*x.gfy[v]);
        const int cats[3] = {0, 1 + Bin(agf, kGf, kNGf), 1 + kNGf + Bin(x.s[v], kSv, kNSv)};
        for (int c = 0; c < 2; ++c) {
          double P[kNP], G[kNG];
          P[0] = c ? x.my[v] : x.mx[v];
          for (int t = 0; t < kNT; ++t) { P[1+t] = c ? x.hy[t][v] : x.hx[t][v]; P[6+t] = -(c ? x.ghy[t][v] : x.ghx[t][v]); G[3+t] = P[6+t]; }
          P[11] = -(c ? x.gcy[v] : x.gcx[v]); P[12] = c ? fvy[v] : fvx[v];
          G[0] = c ? x.gfy[v] : x.gfx[v]; G[1] = c ? x.gmy[v] : x.gmx[v]; G[2] = P[11]; G[8] = P[12];
          for (int ic = 0; ic < 3; ++ic) {
            const int cc = cats[ic];
            reg.n[cc] += 1;
            for (int p = 0; p < kNP; ++p) { reg.pp[cc][p] += P[p]*P[p]; for (int g = 0; g < kNG; ++g) reg.pg[cc][p][g] += P[p]*G[g]; }
            for (int g = 0; g < kNG; ++g) reg.gg[cc][g] += G[g]*G[g];
          }
        }
      }

      // ---- (b) the candidates of single-interaction clusters -------------------------
      for (int k = 0; k < nk; ++k) {
        const int tv = x.tru[k];
        if (tv < 0 || !(x.pt[k] >= kPtB[0])) continue;
        const int b = Bin(x.pt[k], kPtB, kNPtB);
        cand.n[b] += 1;
        for (int m = 0; m < kNMode; ++m) {
          const jva::Cluster &q = CL[m][k];
          if (!q.cand.empty()) cand.has[m][b] += 1;
          if (std::find(q.cand.begin(), q.cand.end(), tv) != q.cand.end()) cand.among[m][b] += 1;
          if (!q.cand.empty() && TrkChoice(q, m, x.pv) == tv) cand.leads[m][b] += 1;
        }
        const jva::Cluster &q = CL[1][k];               // chVtx is the same in every mode
        if (!q.chVtx.empty()) cand.has[3][b] += 1;
        if (std::find(q.chVtx.begin(), q.chVtx.end(), tv) != q.chVtx.end()) cand.among[3][b] += 1;
        if (ArgMax(q.chVtx, q.chPt) == tv) cand.leads[3][b] += 1;
      }

      // ---- (c) the grid and the references -------------------------------------------
      rn2.assign(nv, 0.);
      for (int v = 0; v < nv; ++v) { const double dx = x.mx[v] - x.gmx[v], dy = x.my[v] - x.gmy[v]; rn2[v] = dx*dx + dy*dy; }
      for (int m = 0; m < kNMode; ++m) {
        std::vector<int> &hm = x.host[kRefT0 + m]; hm.assign(nk, jva::kNull);
        for (int k = 0; k < nk; ++k) hm[k] = TrkChoice(CL[m][k], m, jva::kNull);
      }
      for (int i = 0; i < kNRef; ++i) Tally(ref[i], x, x.host[i], i >= kRefT0 ? has[i - kRefT0] : hasF, rn2, fvx, fvy, sx, sy);
      for (int rr = 0; rr < kNRec; ++rr)
        for (int m = 0; m < kNMode; ++m)
          for (int p = 0; p < kNPrice; ++p)
            for (int sg = 0; sg < kNSig; ++sg) {
              jva::Config c;                           // the analyzer's defaults, then the grid point
              if (fk.ok) { c.maxCombos = fk.c.maxCombos; c.maxSweeps = fk.c.maxSweeps; }
              c.candMode = m; c.recoilMode = rr > 0 ? 1 : 0; c.mhtMin = rr > 0 ? kT[rr-1] : 10.;
              c.tauTrk = kPrice[p].tt; c.tauAll = kPrice[p].ta;
              const double k = kSigScale[sg];
              c.sigma0 = 1.0*k; c.sigmaK = 1.0*k; c.sigmaJ = 0.35*k;
              const std::vector<int> h = jva::Assign(V[rr], CL[m], c);
              Tally(grid[Cfg(rr, m, p, sg)], x, h, has[m], rn2, fvx, fvy, sx, sy);
            }
      if (progress) ++(*progress);
    }
  }

  // ---- printing ----------------------------------------------------------------------
  inline double Rms(double s2, double n) { return n > 0 ? std::sqrt(s2/n) : 0.; }
  struct NetRight { double n, eb, ec; };             // right - wrong, binomial error, error from the crossings
  inline NetRight Net(const Grid &g, double nsingle) {
    const double R = g.out[0][kRight] + g.out[1][kRight], W = g.out[0][kWrong] + g.out[1][kWrong];
    NetRight x; x.n = R - W;
    x.eb = nsingle > 0 ? std::sqrt(std::max(0., R + W - (R - W)*(R - W)/nsingle)) : 0.;
    const double N = g.nCross;
    x.ec = N > 1 ? std::sqrt(std::max(0., (g.nrSum2 - g.nrSum*g.nrSum/N))*N/(N - 1)) : 0.;
    return x;
  }
  inline void DRes(const Grid &g, double &d, double &e) {   // d<res^2> per component to none, and its error
    const double N = g.nCross, n = g.resN;
    d = n > 0 ? g.dSum/n : 0.;
    e = (N > 1 && n > 0) ? std::sqrt(std::max(0., g.dSum2 - g.dSum*g.dSum/N)*N/(N - 1))/n : 0.;
  }
  inline void DQ(const Grid &g, double &d, double &e) {     // d|F - A|^2 per crossing to none, and its error
    const double N = g.nCross;
    d = N > 0 ? g.qSum/N : 0.;
    e = N > 1 ? std::sqrt(std::max(0., g.qSum2 - g.qSum*g.qSum/N)/(N - 1)/N) : 0.;
  }
  void PrintRow(const char *name, const Grid &g, double nsingle, double ncombo) {
    double o[kNOut]; for (int i = 0; i < kNOut; ++i) o[i] = g.out[0][i] + g.out[1][i];
    const NetRight nr = Net(g, nsingle); double d, e, q, qe; DRes(g, d, e); DQ(g, q, qe);
    printf("  %-24s %6.0f %6.0f %6.0f %6.0f %6.0f | %6.0f +- %4.0f (%4.0f) %6.2f%% %7.0f | %6.3f %+7.3f +- %5.3f | %6.3f %6.3f | %6.3f %5.3f | %+7.2f +- %4.2f\n",
           name, o[kRight], o[kWrong], o[kNulledSingle], o[kNulledCombo], o[kKeptCombo],
           nr.n, nr.eb, nr.ec, nsingle > 0 ? 100.*nr.n/nsingle : 0., o[kRight] - o[kKeptCombo],
           Rms(g.res2, g.resN), d, e, Rms(g.met2, g.metN), Rms(g.met2pv, g.metNpv),
           Rms(g.res2pv, g.resNpv), ncombo > 0 ? o[kKeptCombo]/ncombo : 0., q, qe);
  }
  void PrintHeader() {
    printf("  %-24s %6s %6s %6s %6s %6s | %25s %7s %7s | %6s %16s | %6s %6s | %6s %5s | %15s\n", "", "right", "wrong", "nulS", "nulC", "keptC",
           "net right +- binom (xing)", "of sing", "R-keptC", "res", "d<res^2> +- err", "MET", "MET pv", "res pv", "keptC", "d|F-A|^2 +- err");
    printf("  %-24s %6s %6s %6s %6s %6s | %25s %7s %7s | %6s %16s | %6s %6s | %6s %5s | %15s\n", "", "", "", "", "", "",
           "", "", "", "[GeV]", "[GeV^2]", "[GeV]", "[GeV]", "[GeV]", "/comb", "[GeV^2/xing]");
  }
  std::string CfgName(int r, int m, int p, int s) {
    char b[64]; snprintf(b, sizeof(b), "%-5s c%d %-5s s%.1f", RecName(r).c_str(), m, kPrice[p].name, kSigScale[s]); return b; }
}

void scanJVA(const char *file = "rootfiles/JVA_v2.root", int nthreads = 8, double maxEntries = 0)
{
  using namespace scanjva;
  TStopwatch sw; sw.Start();
  Knobs fk; Long64_t n = 0;
  {
    std::unique_ptr<TFile> f(TFile::Open(file));
    if (!f || f->IsZombie()) { printf("scanJVA: cannot open %s\n", file); return; }
    TTree *t = (TTree*)f->Get("jva");
    if (!t) { printf("scanJVA: %s has no tuple (storeTuple 0?)\n", file); return; }
    n = t->GetEntries();
    if (maxEntries > 0) n = std::min(n, (Long64_t)maxEntries);
    fk = FileKnobs(f.get());
  }
  printf("scanJVA: %s, %lld crossings, %d threads\n", file, n, nthreads);
  if (fk.ok)
    printf("  production knobs (hopts / %d jobs): candMode %d trkMinPt %g chMinPt %g, recoilMode %d mhtMin %g, "
           "sigma0 %g sigmaK %g sigmaJ %g, tauTrk %g tauAll %g, maxCombos %g maxSweeps %d\n",
           fk.njobs, fk.c.candMode, fk.c.trkMinPt, fk.c.chMinPt, fk.c.recoilMode, fk.c.mhtMin, fk.c.sigma0, fk.c.sigmaK,
           fk.c.sigmaJ, fk.c.tauTrk, fk.c.tauAll, fk.c.maxCombos, fk.c.maxSweeps);
  else printf("  WARNING: no hist/hopts, the stored hosts are not checked\n");
  if (fk.ok && fk.c.trkMinPt > kTrkMin + 1e-9)
    printf("  WARNING: the production's trkMinPt %g is above the grid's %g: the grid's track candidates are those of %g\n",
           fk.c.trkMinPt, kTrkMin, fk.c.trkMinPt);

  ROOT::EnableThreadSafety();
  nthreads = std::max(1, nthreads);
  std::vector< std::unique_ptr<Worker> > w;
  std::atomic<long> progress(0);
  for (int i = 0; i < nthreads; ++i) {
    w.emplace_back(new Worker());
    w.back()->file = file; w.back()->fk = fk; w.back()->progress = &progress;
    w.back()->first = (Long64_t)(double(i)*n/nthreads); w.back()->last = (Long64_t)(double(i+1)*n/nthreads);
  }
  std::vector<std::thread> th;
  for (int i = 0; i < nthreads; ++i) th.emplace_back([&w, i]() { w[i]->Run(); });
  for (auto &t : th) t.join();
  Worker A; A.fk = fk;
  for (int i = 0; i < nthreads; ++i) {
    if (!w[i]->err.empty()) { printf("scanJVA: thread %d: %s\n", i, w[i]->err.c_str()); return; }
    AddTo(A.reg, w[i]->reg); AddTo(A.cand, w[i]->cand);
    for (int c = 0; c < kNCfg; ++c) AddTo(A.grid[c], w[i]->grid[c]);
    for (int c = 0; c < kNRef; ++c) AddTo(A.ref[c], w[i]->ref[c]);
    A.nCross += w[i]->nCross; A.nCl += w[i]->nCl; A.optMis += w[i]->optMis; A.hostMis += w[i]->hostMis;
    A.hostMisCross += w[i]->hostMisCross; A.trkMis += w[i]->trkMis; A.validated += w[i]->validated;
  }
  printf("  %.0f crossings, %.2f clusters per crossing, %.0f s\n", A.nCross, A.nCl/std::max(1., A.nCross), sw.RealTime());

  // ---- (0) ----------------------------------------------------------------------------
  printf("\n(0) the offline rebuild at the production's knobs against what the analyzer stored\n");
  if (fk.ok) {
    printf("  candidates (cl_nopt) mismatching: %.0f of %.0f clusters\n", A.optMis, A.nCl);
    printf("  trk choice (cl_host_trk) mismatching: %.0f of %.0f clusters%s\n", A.trkMis, A.nCl,
           fk.c.candMode == 2 ? " (candMode 2: a T_ku below trkMinPt is not in the tuple)" : "");
    if (A.validated > 0)
      printf("  jva hosts (cl_host_jva) mismatching: %.0f of %.0f clusters, in %.0f of %.0f crossings\n",
             A.hostMis, A.nCl, A.hostMisCross, A.validated);
    else printf("  jva hosts not checked: the production's recoil (mhtMin %g) is not one of the tuple's thresholds\n", fk.c.mhtMin);
  }

  // ---- (a) ----------------------------------------------------------------------------
  const Reg &R = A.reg;
  auto slope = [&](int c, int p, int g) { return R.gg[c][g] > 0 ? R.pg[c][p][g]/R.gg[c][g] : 0.; };
  auto corr  = [&](int c, int p, int g) { return R.pp[c][p] > 0 && R.gg[c][g] > 0 ? R.pg[c][p][g]/std::sqrt(R.pp[c][p]*R.gg[c][g]) : 0.; };
  auto resid = [&](int c, int p, int g) { const double b = slope(c,p,g); return R.n[c] > 0 ? std::sqrt(std::max(0., (R.pp[c][p] - 2*b*R.pg[c][p][g] + b*b*R.gg[c][g])/R.n[c])) : 0.; };
  auto reso  = [&](int c, int p, int g) {   // rms of g - (sum pg / sum pp) p
    const double a = R.pp[c][p] > 0 ? R.pg[c][p][g]/R.pp[c][p] : 0.;
    return R.n[c] > 0 ? std::sqrt(std::max(0., (R.gg[c][g] - 2*a*R.pg[c][p][g] + a*a*R.pp[c][p])/R.n[c])) : 0.; };
  auto rmsP = [&](int c, int p) { return R.n[c] > 0 ? std::sqrt(R.pp[c][p]/R.n[c]) : 0.; };
  auto rmsG = [&](int c, int g) { return R.n[c] > 0 ? std::sqrt(R.gg[c][g]/R.n[c]) : 0.; };
  printf("\n(a) the recoil R_v against the truth, %.0f owned vertices, both components pooled (fits through 0)\n", R.nown);
  printf("  sign convention: M_v and MHT_T = -(sum of reco p), 'what is missing'; genfwd = +(owner's generated pT at 2.5 < |eta| < 5);\n"
         "  genmet = MET_true(owner) = -(gencen + genfwd); a perfect recoil is R_v = genfwd + genmet = -gencen\n");
  printf("  rms per component: genfwd %.2f, genmet %.2f, gencen %.2f, F_v %.2f GeV; owned vertices with a jet above T:",
         rmsG(0,0), rmsG(0,1), rmsG(0,2), rmsG(0,8));
  for (int t = 0; t < kNT; ++t) printf(" %d: %.3f", kT[t], R.nown > 0 ? R.nz[t]/R.nown : 0.);
  printf("\n");
  const int tg[4] = {0, 1, 2, 8};
  printf("  %-22s %6s |", "R_v", "rms");
  for (int i = 0; i < 4; ++i) { char h[64]; snprintf(h, sizeof(h), "vs %s: b  r  resid  reso", GName(tg[i])); printf(" %-29s|", h); }
  printf("\n");
  for (int p = 0; p < 6; ++p) {
    printf("  %-22s %6.2f |", PName(p), rmsP(0,p));
    for (int i = 0; i < 4; ++i) printf(" %6.3f %6.3f %6.2f %6.2f |", slope(0,p,tg[i]), corr(0,p,tg[i]), resid(0,p,tg[i]), reso(0,p,tg[i]));
    printf("\n");
  }
  printf("  the truth chain (what an ideal central measurement would give):\n");
  for (int p = 6; p < kNP; ++p) {
    printf("  %-22s %6.2f |", PName(p), rmsP(0,p));
    for (int i = 0; i < 4; ++i) printf(" %6.3f %6.3f %6.2f %6.2f |", slope(0,p,tg[i]), corr(0,p,tg[i]), resid(0,p,tg[i]), reso(0,p,tg[i]));
    printf("\n");
  }
  printf("  R_v against genfwd in bins of |genfwd| [GeV] and of S_v [GeV]: slope b / correlation r / resid rms(R - b genfwd) / reso rms(genfwd - a R)\n");
  printf("  %-12s %7s |", "bin", "n/2");
  for (int p = 0; p < 6; ++p) printf(" %-27s|", p == 0 ? "M_v" : PName(p));
  printf("\n");
  for (int c = 1; c < kNCatA; ++c) {
    const bool gf = c <= kNGf; const int b = gf ? c - 1 : c - 1 - kNGf;
    const double lo = gf ? kGf[b] : kSv[b], hi = gf ? kGf[b+1] : kSv[b+1];
    char lab0[64];
    if (hi > 1e29) snprintf(lab0, sizeof(lab0), "%s > %g", gf ? "|gf|" : "S", lo);
    else snprintf(lab0, sizeof(lab0), "%s %g-%g", gf ? "|gf|" : "S", lo, hi);
    const std::string lab = lab0;
    printf("  %-12s %7.0f |", lab.c_str(), R.n[c]/2);
    for (int p = 0; p < 6; ++p) printf(" %5.2f %5.2f %6.2f %6.2f |", slope(c,p,0), corr(c,p,0), resid(c,p,0), reso(c,p,0));
    printf("\n");
  }
  printf("  R_v against F_v (the oracle's clusters) in the same bins: slope b / correlation r\n");
  printf("  %-12s %7s |", "bin", "n/2");
  for (int p = 0; p < 6; ++p) printf(" %-13s|", p == 0 ? "M_v" : PName(p));
  printf("\n");
  for (int c = 1; c < kNCatA; ++c) {
    const bool gf = c <= kNGf; const int b = gf ? c - 1 : c - 1 - kNGf;
    const double lo = gf ? kGf[b] : kSv[b], hi = gf ? kGf[b+1] : kSv[b+1];
    char lab0[64];
    if (hi > 1e29) snprintf(lab0, sizeof(lab0), "%s > %g", gf ? "|gf|" : "S", lo);
    else snprintf(lab0, sizeof(lab0), "%s %g-%g", gf ? "|gf|" : "S", lo, hi);
    printf("  %-12s %7.0f |", lab0, R.n[c]/2);
    for (int p = 0; p < 6; ++p) printf(" %6.3f %6.3f|", slope(c,p,8), corr(c,p,8));
    printf("\n");
  }
  printf("  MHT_T (reco, raw threshold) against the truth jets -genmht_T' (gen threshold): slope / correlation\n");
  printf("  %-8s", "");
  for (int t2 = 0; t2 < kNT; ++t2) { char h[32]; snprintf(h, sizeof(h), "-genmht%d", kT[t2]); printf(" %14s", h); }
  printf("\n");
  for (int t = 0; t < kNT; ++t) {
    printf("  %-8s", PName(1 + t));
    for (int t2 = 0; t2 < kNT; ++t2) printf("  %5.3f / %5.3f", slope(0,1+t,3+t2), corr(0,1+t,3+t2));
    printf("\n");
  }
  { int best = 0; for (int p = 1; p < 6; ++p) if (corr(0,p,0) > corr(0,best,0)) best = p;
    printf("  largest correlation with genfwd: %s (r = %.3f); with F_v: ", PName(best), corr(0,best,0));
    int bf = 0; for (int p = 1; p < 6; ++p) if (corr(0,p,8) > corr(0,bf,8)) bf = p;
    printf("%s (r = %.3f)\n", PName(bf), corr(0,bf,8)); }

  // ---- (b) ----------------------------------------------------------------------------
  const Cand &C = A.cand;
  printf("\n(b) single-interaction clusters (cl_true >= 0): the true vertex among the candidates / carrying the most evidence\n");
  printf("  %-10s %8s | %-21s | %-21s | %-21s | %-21s\n", "pT [GeV]", "clusters", "candMode 0 (T>=0.5)", "candMode 1 (C>=0.5)",
         "candMode 2 (either)", "charged argmax (C>0)");
  printf("  %-10s %8s |", "", "");
  for (int q = 0; q < kNRule; ++q) printf(" %6s %6s %6s |", "has", "among", "leads");
  printf("\n");
  for (int b = 0; b < kNPtB; ++b) {
    char lab0[64];
    if (b == kNPtB-1) snprintf(lab0, sizeof(lab0), "> %g", kPtB[b]); else snprintf(lab0, sizeof(lab0), "%g-%g", kPtB[b], kPtB[b+1]);
    const std::string lab = lab0;
    printf("  %-10s %8.0f |", lab.c_str(), C.n[b]);
    for (int q = 0; q < kNRule; ++q) {
      const double nn = C.n[b] > 0 ? C.n[b] : 1;
      printf(" %6.3f %6.3f %6.3f |", C.has[q][b]/nn, C.among[q][b]/nn, C.leads[q][b]/nn);
    }
    printf("\n");
  }
  { double nn = 0, am[kNRule] = {0}, ld[kNRule] = {0};
    for (int b = 1; b < kNPtB; ++b) { nn += C.n[b]; for (int q = 0; q < kNRule; ++q) { am[q] += C.among[q][b]; ld[q] += C.leads[q][b]; } }
    printf("  %-10s %8.0f |", "> 5", nn);
    for (int q = 0; q < kNRule; ++q) printf(" %6s %6.3f %6.3f |", "", nn > 0 ? am[q]/nn : 0., nn > 0 ? ld[q]/nn : 0.);
    printf("\n  (binomial error on a fraction f of n clusters: sqrt(f(1-f)/n), 0.004-0.008 above 5 GeV)\n"); }

  // ---- (c) ----------------------------------------------------------------------------
  const Grid &gn = A.ref[kRefNone];
  const double nsingle = gn.out[0][kNulledSingle] + gn.out[1][kNulledSingle];
  const double ncombo  = gn.out[0][kNulledCombo] + gn.out[1][kNulledCombo];
  printf("\n(c) the grid: clusters above %g GeV, %.0f single-interaction and %.0f combinations (none: all nulled)\n", kClMin, nsingle, ncombo);
  printf("  net right = right - wrong, error binomial (and from the spread over crossings); R-keptC = right - kept_combo, the change of\n"
         "  'correct' (right + nulled_combo) against none; res = rms per component of MET_v - MET_true(owner),\n"
         "  owned vertices; d<res^2> = paired difference to none, per component; MET = rms per component of MET_v, all vertices / PV;\n"
         "  d|F-A|^2 = sum over the vertices of |F_v - A_v|^2 - |F_v|^2 per crossing, A_v the pT the method gave v, F_v what ojet gives it\n"
         "  (the energy-weighted net right: < 0 when the forward pT handed out is more the right vertex's than noise; ojet: -sum |F_v|^2)\n"
         "  with tauTrk = 0 and tauAll = 1e9 every price is 0 or never paid, X scales as a whole with sigma^2, and s1.0 = s1.5 exactly\n");
  PrintHeader();
  for (int i = 0; i < kNRef; ++i) PrintRow(kRefName[i], A.ref[i], nsingle, ncombo);
  for (int rr = 0; rr < kNRec; ++rr) {
    printf("  --\n");
    for (int m = 0; m < kNMode; ++m)
      for (int p = 0; p < kNPrice; ++p)
        for (int sg = 0; sg < kNSig; ++sg)
          PrintRow(CfgName(rr, m, p, sg).c_str(), A.grid[Cfg(rr, m, p, sg)], nsingle, ncombo);
  }
  printf("  the same split by whether the cluster has candidates under the rule (right wrong nulS nulC keptC):\n");
  for (int i = 0; i < kNRef; ++i) {
    printf("  %-24s with:", kRefName[i]); for (int o = 0; o < kNOut; ++o) printf(" %6.0f", A.ref[i].out[0][o]);
    printf("   without:"); for (int o = 0; o < kNOut; ++o) printf(" %6.0f", A.ref[i].out[1][o]); printf("\n");
  }
  for (int rr = 0; rr < kNRec; ++rr) for (int m = 0; m < kNMode; ++m) for (int p = 0; p < kNPrice; ++p) for (int sg = 0; sg < kNSig; ++sg) {
    const Grid &g = A.grid[Cfg(rr, m, p, sg)];
    printf("  %-24s with:", CfgName(rr, m, p, sg).c_str()); for (int o = 0; o < kNOut; ++o) printf(" %6.0f", g.out[0][o]);
    printf("   without:"); for (int o = 0; o < kNOut; ++o) printf(" %6.0f", g.out[1][o]); printf("\n");
  }

  // ---- (d) ----------------------------------------------------------------------------
  const double resNone = Rms(gn.res2, gn.resN);
  int best = -1, bestAny = -1;
  for (int c = 0; c < kNCfg; ++c) {
    const double nr = Net(A.grid[c], nsingle).n;
    if (bestAny < 0 || nr > Net(A.grid[bestAny], nsingle).n) bestAny = c;
    double d, e; DRes(A.grid[c], d, e);
    if (d <= 0 && (best < 0 || nr > Net(A.grid[best], nsingle).n)) best = c;
  }
  auto name = [&](int c) { const int sg = c % kNSig, p = (c / kNSig) % kNPrice, m = (c / kNSig / kNPrice) % kNMode, rr = c / kNSig / kNPrice / kNMode;
                           return CfgName(rr, m, p, sg); };
  printf("\n(d) none: res %.3f GeV per component; the largest net right of all: %s, %.0f\n", resNone, name(bestAny).c_str(), Net(A.grid[bestAny], nsingle).n);
  if (best >= 0) {
    const Grid &g = A.grid[best]; const NetRight nr = Net(g, nsingle); double d, e; DRes(g, d, e);
    printf("  the largest net right with d<res^2> <= 0: %s: net right %.0f +- %.0f binomial (+- %.0f from the crossings) of %.0f single clusters"
           " = %.2f +- %.2f %%, %.1f sigma from none; res %.3f GeV (d<res^2> %+.3f +- %.3f GeV^2)\n",
           name(best).c_str(), nr.n, nr.eb, nr.ec, nsingle, 100.*nr.n/nsingle, 100.*nr.eb/nsingle, nr.eb > 0 ? nr.n/nr.eb : 0.,
           Rms(g.res2, g.resN), d, e);
    printf("  its outcomes by cluster pT (right wrong nulS nulC keptC): 5-10 GeV");
    for (int o = 0; o < kNOut; ++o) printf(" %.0f", g.outPt[0][o]);
    printf(", 10-20 GeV"); for (int o = 0; o < kNOut; ++o) printf(" %.0f", g.outPt[1][o]);
    printf(", > 20 GeV"); for (int o = 0; o < kNOut; ++o) printf(" %.0f", g.outPt[2][o]);
    printf("\n");
  } else printf("  no configuration keeps d<res^2> <= 0\n");
  printf("  the largest net right per recoil (any residual):\n");
  PrintHeader();
  for (int rr = 0; rr < kNRec; ++rr) {
    int b = -1;
    for (int m = 0; m < kNMode; ++m) for (int p = 0; p < kNPrice; ++p) for (int sg = 0; sg < kNSig; ++sg) {
      const int c = Cfg(rr, m, p, sg);
      if (b < 0 || Net(A.grid[c], nsingle).n > Net(A.grid[b], nsingle).n) b = c; }
    PrintRow(name(b).c_str(), A.grid[b], nsingle, ncombo);
  }
  // the ten largest net right, whatever the residual
  std::vector<int> idx(kNCfg); for (int c = 0; c < kNCfg; ++c) idx[c] = c;
  std::stable_sort(idx.begin(), idx.end(), [&](int a, int b) { return Net(A.grid[a], nsingle).n > Net(A.grid[b], nsingle).n; });
  printf("  the ten largest net right:\n");
  PrintHeader();
  for (int i = 0; i < 10; ++i) PrintRow(name(idx[i]).c_str(), A.grid[idx[i]], nsingle, ncombo);
  printf("scanJVA: done in %.0f s\n", sw.RealTime());
}
