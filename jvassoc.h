// genseed - Copyright 2026 Mikko Voutilainen - Apache License 2.0
#ifndef genseed_jvassoc_h
#define genseed_jvassoc_h
// jvassoc.h - jet-vertex association (JVA) of the forward, vertex-blind energy
// by a global minimisation of the missing transverse momentum over all vertices.
// (Planned as jva.h; the default macOS file system is case-insensitive, where
// jva.h and the analyzer's JVA.h are one file, so the ROOT-free half is
// jvassoc.h.  The namespace is jva.)
//
// THE PROBLEM.  PUPPI tells vertices apart with tracks, and past the tracker it
// has none: every PF candidate with |eta| > 2.5 gets the same weight wFwd at
// every vertex hypothesis (puppi.h, FORWARD REGIONS ARE VERTEX-BLIND).  Run per
// vertex, as genseed runs it, every vertex therefore carries the same forward
// jets.  On 300 crossings of MC24NanoV15_PU_IT_OOT (<N_PU> = 45): 3.6 jets
// above 5 GeV at 2.5 < |eta| < 3.0 and 1.7 at 3 < |eta| < 5 per vertex, 98-100%
// of their pT from |eta| > 2.5 candidates.  The part of their cone that is
// still inside the tracker holds 4.2 tracks per jet, 4.1 of them associated to
// OTHER vertices than the one under study.  Nor are they jets of one
// interaction: the dominant interaction gives less than 0.4 of the linked
// generated pT for 71% of them and more than 0.6 for 7%, with 3.2 interactions
// above 10% each.  That is the 2.5 < |eta| < 3.0 spike of the genseed
// production (reco/gen = 2-5 at 10-24 GeV), and it is what the MET sees: the
// PUPPI MET per vertex has an rms of 25 GeV at the PV and 23 GeV at the other
// vertices, against 6 GeV for the true visible imbalance of the interaction.
//
// THE IDEA.  The forward energy is not wrong, it is unassigned.  Partition it
// among the reconstructed vertices, whole clusters at a time: anti-kT R of the
// vertex-blind candidates at wFwd p_i (every one of them in exactly one
// cluster F_k), and each cluster goes to the vertex whose recoil it balances
// best - or to nobody.  The vertex-resolved part of the event (|eta| <= etaFwd
// with the vertex's own PUPPI weights) has a MET of its own,
//     M_v = -sum_{|eta_i| <= etaFwd} w_i(v) p_i ,     S_v = sum w_i(v) pT_i ,
// and a cluster that belongs to v is, to the resolution of M_v, part of the
// recoil M_v is missing.  The assignment h (h_k a vertex index or null = -1)
// minimises
//     X(h) = sum_v |M_v - sum_{k: h_k = v} p_k|^2 / sigma_v^2  +  sum_{k: h_k != null} tau_k ,
//     sigma_v^2 = sigma0^2 + sigmaK^2 S_v   [GeV^2 per component],
// over ALL vertices at once: a cluster cannot balance two vertices, and two
// clusters that each half-balance one vertex are judged together.  The
// resolution model is measured, not assumed: the vertex-resolved MET against
// the owner interaction's true central recoil has sigma ~= sqrt(1 + 1.0^2 S_v)
// GeV per component, 1.5 GeV at S = 3, 5.8 at S = 30, 9.3 at S = 90.
//
// THE PRICE tau, AND WHY TRACKS ARE NEEDED.  Moving cluster k from null to v
// changes X by (|p_k|^2 - 2 M_v.p_k)/sigma_v^2 + tau_k: it is kept only if it
// explains more of M_v than it costs.  Without a price every cluster would find
// some vertex among forty whose noise it happens to reduce.  The price alone is
// not enough, though: the true forward recoil of an interaction is 4.8 GeV rms
// per component, about the MET noise of a typical vertex, so MET balance alone
// can place only forward jets well above that noise.  Tracks narrow the
// choice.  The forward clusters at the tracker edge (2.5 < |eta| < ~2.9) have
// tracks inside their cone, and the vertices those tracks belong to,
//     T_ku = sum of pT of vertex u's charged candidates with |eta| <= etaFwd
//            within dR < R of the cluster axis,  candidate if T_ku >= trkMinPt,
// are the only vertices such a cluster may go to, at price tauTrk (4), besides
// null.  A cluster with no track at all (HF) may go to any vertex at price
// tauAll (12): it has to be a strong balance to be believed.
//
// THE CANDIDATES, candMode.  The in-cone tracks turned out a weak handle (WHAT
// IT DOES NOT FIND, below): the tracker part of a cone at the tracker edge is
// mostly other interactions' tracks, which is why they are in the cone at
// all.  The cluster has a better one inside it.  Its charged constituents at
// 2.5 < |eta| < ~3 still have a track and a vertexRef, which PUPPI ignores
// past etaTracker (they get wFwd like the neutrals), and they are part of the
// cluster's own energy, not a neighbour's:
//     C_ku = sum of the raw pT of cluster k's charged constituents with
//            vertexRef u  (raw, not wFwd pT: a track's vertex does not
//            depend on how much PUPPI believes its energy).
// On 300 crossings the leading C_ku is the true vertex of 40, 49, 31-36 and
// 21-26% of the single-interaction clusters at 0-2, 2-5, 5-10 and above 10
// GeV, the leading T_ku of 3-7%.  candMode 0 takes the u with T_ku >=
// trkMinPt as candidates (v1), 1 those with C_ku >= chMinPt, 2 either.  The
// price is tauTrk whichever evidence made the candidate, and a cluster with no
// candidate still goes anywhere at tauAll.  Cluster keeps the evidence of
// both kinds whatever the mode, and Evidence() ranks the candidates for the
// analyzer's trk method: T_ku, C_ku, or T_ku + C_ku (vertex u's charged pT in
// and around the cluster).
//
// THE RECOIL, recoilMode.  The other half of the trouble: at the true vertex
// M_v carries 0.5-0.7 of the cluster's pT along it, against ~11 GeV of noise
// per component.  The particle MET sums every vertex-resolved candidate, and
// at 45 interactions per crossing most of S_v at any vertex is soft: the
// vertex's own underlying event, and what PUPPI lets through of the others'
// neutrals and mis-associated tracks, at a fraction of their pT and in random
// directions.  That is the sqrt(S_v) noise, and it carries next to none of
// the balance of a 10 GeV forward jet; the recoil that does is in the
// vertex's jets.  recoilMode 1 therefore takes the MHT of the vertex-resolved
// jets - the jets of the event made of the vertex-resolved candidates alone,
// at w_i(v), no vertex-blind candidate in it (the analyzer's "none" event) -
// above mhtMin in pT^raw:
//     M_v = -sum_{j: pT_j > mhtMin} p_j ,   HT_v = sum_j pT_j ,
//     sigma_v^2 = sigma0^2 + sigmaJ^2 sum_j pT_j^2 + sigmaK^2 (S_v - HT_v) .
// The jets bring their own resolution, sigmaJ ~ 0.35 of the raw pT per jet
// (the raw jet resolution at 10-20 GeV, put on either component: a prior, not
// a fit); the unclustered rest is out of M_v, and with it the part of the
// true recoil it carried, which the MET must now miss by about what that rest
// weighed: the same sqrt(S) term on S_v - HT_v.  Raw on both sides, as in
// recoilMode 0: p_k is raw (at wFwd) and so is M_v.  recoilMode 0 is the
// particle MET with HT_v = sum pT_j^2 = 0, i.e. v1 bit for bit.  mhtMin
// trades noise against recoil - below ~5 GeV the soft noise is back, above ~15
// a PU vertex has no jet left - and is to be chosen on the MHT against the
// generated recoil of the owner interaction, which the analyzer's tuple
// carries at 3, 5, 7, 10 and 15 GeV.
//
// THE OPTIMISER.  Two clusters interact in X only through a shared vertex, so
// the clusters fall into connected components (union-find over "k may go to
// v"; a cluster without tracks connects to every vertex).  A component is
// solved exhaustively over the product of its options (candidates + null) when
// that has at most maxCombos combinations - an odometer, one residual update
// per step, so 2e5 combinations cost a few ms - and otherwise by coordinate
// descent: every cluster starts at its best single move against M_v alone
// (as if it were the only one), then each in turn moves to its best option
// given the others, sweep after sweep until a sweep changes nothing (or
// maxSweeps), and then a pairwise pass re-optimises every pair of clusters
// that share a vertex JOINTLY over the product of their options (which
// contains the swap of their two hosts, and also both moving elsewhere); if
// that improves anything the descent resumes.  The joint move is exact in O(1)
// per combination: with both removed and residual R_v, putting both at v costs
// |R_v - p_k - p_l|^2 = single(k) + single(l) + 2 p_k.p_l.  Pairs with
// disjoint options are skipped - their joint optimum is the product of the
// single ones, which the descent already has.  Everything runs in a fixed
// order, a move is made only for a strict improvement (> kEps), and ties go to
// null, then to the lowest vertex index, so the result is deterministic.
// Stats reports the objective, the components, how many were solved exactly,
// and the sweeps.
//
// WHAT IT CANNOT DO.  A cluster is assigned whole.  A forward jet made of three
// interactions balances none of them and goes to null together with its
// single-interaction part: combinations are dropped, not split.  Splitting
// would need something that tells interactions apart inside one calorimeter
// cluster, and there is nothing (precision timing would be).  A cluster whose
// true vertex is not among its candidates (no track in its cone in candMode
// 0, no charged constituent in 1), but which has another vertex's, can only
// go to that vertex or to null.  The descent is a local
// search and on components beyond maxCombos it may stop short of the optimum:
// on random instances of 4-9 clusters forced onto it, it reached the
// brute-force minimum in 972 of 979; on real crossings, where the clusters
// without tracks join everything into one component (24 clusters, 34
// vertices) and all but 1-3 per mille go to the descent, an iterated local
// search started from its answer found nothing lower in 4000 crossings.
// And the prices are priors, not probabilities: tauTrk and tauAll are knobs.
//
// WHAT IT DOES NOT FIND, measured on 4000 crossings at the v1 knobs
// (candMode 0, recoilMode 0), and why the two knobs above exist.
// The tracker-edge tracks in the cone are mostly other vertices': of the
// single-interaction clusters above 5 GeV that have track candidates, the
// true vertex is among them for 31% and has the largest T_ku for 21%.  And
// at the true vertex M_v carries only 0.5-0.7 of a 5-40 GeV cluster's pT
// along its direction (the raw central response, the interaction's other
// forward energy), so the chi2 gain of the right assignment, 2 M_v.p_k -
// |p_k|^2, averages about zero, against 11 GeV of MET noise per component
// across the cluster at those vertices (sigma_v ~ 9 GeV): nearly every
// cluster goes to null, and neither a balance factor on p_k nor lower prices
// make right outnumber wrong.
//
// ROOT-free, like puppi.h, genlink.h and tiledjet.h; BuildClusters uses
// tiledjet.h whatever backend the analyzer clusters its jets with.
#include "puppi.h"
#include "tiledjet.h"
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

namespace jva {

  struct Config {
    double R          = 0.4;     // forward clusters, anti-kT, and the track cone around their axis
    double etaFwd     = 2.5;     // |eta| > etaFwd is vertex-blind
    int    candMode   = 0;       // candidates: 0 T_ku >= trkMinPt, 1 C_ku >= chMinPt, 2 either
    double trkMinPt   = 0.5;     // GeV, T_ku (in-cone tracker tracks) that makes u a candidate
    double chMinPt    = 0.5;     // GeV, C_ku (in-cluster charged constituents) that makes u a candidate
    double sigma0     = 1.0;     // GeV, sigma_v^2 = sigma0^2 + sigmaJ^2 sum pT_j^2 + sigmaK^2 (S_v - HT_v)
    double sigmaK     = 1.0;     // GeV^1/2
    double sigmaJ     = 0.35;    // per recoil jet (recoilMode 1)
    int    recoilMode = 0;       // 0 particle MET of the vertex-resolved part, 1 MHT of its jets (the analyzer builds it)
    double mhtMin     = 10.;     // GeV, pT^raw of the recoil jets (recoilMode 1)
    double tauTrk     = 4.0;     // price of a candidate
    double tauAll     = 12.0;    // price of any vertex for a cluster without candidates
    double maxCombos  = 200000;  // exhaustive search up to this many combinations per component
    int    maxSweeps  = 50;      // descent: coordinate sweeps + pair passes per component
  };

  // A usable vertex: the recoil its forward clusters are to balance (GeV), the
  // scalar sum S_v of its vertex-resolved part, and for recoilMode 1 how much
  // of S_v is in the recoil jets (ht) and their sum of pT^2 (sj2); both 0 for
  // the particle MET, which makes Sigma2 v1's sigma0^2 + sigmaK^2 S_v exactly.
  struct Vertex { double mx = 0, my = 0, s = 0, ht = 0, sj2 = 0; };

  // A forward cluster: its transverse momentum (the sum of wFwd p_i), axis,
  // its candidates under candMode (vertex indices, ascending) with the T_ku
  // and C_ku of each, and the evidence of both kinds whatever candMode: the u
  // with T_ku >= trkMinPt (v1's track candidates) and the u with C_ku > 0,
  // ascending.  An empty cand means no candidate: any vertex at tauAll.
  struct Cluster {
    double px = 0, py = 0, pt = 0, eta = 0, phi = 0;
    std::vector<int>    cand;
    std::vector<double> candTrkPt, candChPt;
    std::vector<int>    trkVtx, chVtx;
    std::vector<double> trkPt, chPt;
  };

  struct Stats {
    double objective = 0;        // X(h) of the result, all vertices
    int components   = 0;        // components with at least one cluster
    int exact        = 0;        // of which solved exhaustively
    int sweeps       = 0;        // coordinate sweeps + pair passes, summed over the descent components
    int largest      = 0;        // clusters in the largest component
  };

  const int    kNull = -1;
  const double kEps  = 1e-9;     // a move must gain more than this
  const double kInf  = std::numeric_limits<double>::infinity();

  // (v1's expression is its own statement and the jet term is added only when
  // there is one, so that no fused multiply-add can move recoilMode 0 by an ulp)
  inline double Sigma2(const Vertex &v, const Config &c) {
    double s2 = c.sigma0*c.sigma0 + c.sigmaK*c.sigmaK*std::max(0., v.s - v.ht);
    if (v.sj2 > 0) s2 += c.sigmaJ*c.sigmaJ*v.sj2;
    return s2;
  }

  // recoilMode 1: the recoil of v is the MHT of its vertex-resolved jets
  // above mhtMin (pT, phi of each; any order), see THE RECOIL.  S_v stays.
  inline void SetMHT(Vertex &v, const std::vector<double> &jpt, const std::vector<double> &jphi, double mhtMin) {
    v.mx = v.my = v.ht = v.sj2 = 0;
    for (size_t j = 0; j < jpt.size(); ++j) {
      if (!(jpt[j] > mhtMin)) continue;
      v.mx -= jpt[j]*std::cos(jphi[j]); v.my -= jpt[j]*std::sin(jphi[j]);
      v.ht += jpt[j]; v.sj2 += jpt[j]*jpt[j];
    }
  }

  // The weight of candidate j of cluster k in the trk method's ranking: the
  // evidence that made it a candidate (T_ku, C_ku, or both summed).
  inline double Evidence(const Cluster &k, size_t j, const Config &c) {
    return c.candMode == 1 ? k.candChPt[j] : (c.candMode == 2 ? k.candTrkPt[j] + k.candChPt[j] : k.candTrkPt[j]); }

  // tau of cluster k at vertex v, +inf where k may not go
  inline double Price(const Cluster &k, int v, int nv, const Config &c) {
    if (v < 0 || v >= nv) return kInf;
    if (k.cand.empty()) return c.tauAll;
    for (size_t j = 0; j < k.cand.size(); ++j) if (k.cand[j] == v) return c.tauTrk;
    return kInf;
  }

  // X(h), from scratch.  h[k] < 0 is null; an assignment that is not allowed
  // (a vertex that is not among a tracked cluster's candidates, an index out
  // of range) is +inf.
  inline double Objective(const std::vector<Vertex> &vs, const std::vector<Cluster> &ks,
                          const std::vector<int> &h, const Config &c)
  {
    const int nv = (int)vs.size(), nk = (int)ks.size();
    if ((int)h.size() != nk) return kInf;
    std::vector<double> rx(nv), ry(nv);
    for (int v = 0; v < nv; ++v) { rx[v] = vs[v].mx; ry[v] = vs[v].my; }
    double x = 0;
    for (int k = 0; k < nk; ++k) {
      if (h[k] < 0) continue;
      const double t = Price(ks[k], h[k], nv, c);
      if (!(t < kInf)) return kInf;
      x += t; rx[h[k]] -= ks[k].px; ry[h[k]] -= ks[k].py;
    }
    for (int v = 0; v < nv; ++v) x += (rx[v]*rx[v] + ry[v]*ry[v])/Sigma2(vs[v], c);
    return x;
  }

  namespace detail {
    struct UnionFind {
      std::vector<int> p;
      explicit UnionFind(int n) : p(n) { for (int i = 0; i < n; ++i) p[i] = i; }
      int find(int x) { while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; } return x; }
      void unite(int a, int b) { a = find(a); b = find(b); if (a != b) { if (a < b) p[b] = a; else p[a] = b; } }
    };

    // The state of one Assign call: options, residuals, per-vertex chi2 terms.
    // Digits: d = 0 is null, d = j >= 1 is opt[k][j-1].
    struct Solver {
      const std::vector<Cluster> &ks; const Config &c;
      std::vector< std::vector<int> > opt; std::vector<double> tau;
      std::vector<double> rx, ry, isg, term;
      double cur = 0;                          // X - X(all null) over what has been moved
      Solver(const std::vector<Vertex> &vs, const std::vector<Cluster> &k, const Config &cc)
        : ks(k), c(cc)
      {
        const int nv = (int)vs.size(), nk = (int)ks.size();
        rx.resize(nv); ry.resize(nv); isg.resize(nv); term.resize(nv);
        for (int v = 0; v < nv; ++v) {
          rx[v] = vs[v].mx; ry[v] = vs[v].my; isg[v] = 1./Sigma2(vs[v], c);
          term[v] = (rx[v]*rx[v] + ry[v]*ry[v])*isg[v];
        }
        opt.resize(nk); tau.resize(nk);
        for (int q = 0; q < nk; ++q) {
          if (ks[q].cand.empty()) {
            opt[q].resize(nv); for (int v = 0; v < nv; ++v) opt[q][v] = v;
            tau[q] = c.tauAll;
          } else {
            for (size_t j = 0; j < ks[q].cand.size(); ++j)
              if (ks[q].cand[j] >= 0 && ks[q].cand[j] < nv) opt[q].push_back(ks[q].cand[j]);
            std::sort(opt[q].begin(), opt[q].end());
            opt[q].erase(std::unique(opt[q].begin(), opt[q].end()), opt[q].end());
            tau[q] = c.tauTrk;
          }
        }
      }
      int host(int q, int d) const { return d == 0 ? kNull : opt[q][d-1]; }
      void setv(int v, double dx, double dy) {
        rx[v] += dx; ry[v] += dy;
        const double t = (rx[v]*rx[v] + ry[v]*ry[v])*isg[v];
        cur += t - term[v]; term[v] = t;
      }
      // cluster q from vertex a to vertex b (either may be null)
      void move(int q, int a, int b) {
        if (a == b) return;
        if (a >= 0) { setv(a,  ks[q].px,  ks[q].py); cur -= tau[q]; }
        if (b >= 0) { setv(b, -ks[q].px, -ks[q].py); cur += tau[q]; }
      }
      // what putting q (currently nowhere) at v would change: chi2 + tau, 0 for null
      double single(int q, int v) const {
        if (v < 0) return 0.;
        const double ax = rx[v] - ks[q].px, ay = ry[v] - ks[q].py;
        return (ax*ax + ay*ay)*isg[v] - term[v] + tau[q];
      }

      // Every combination of the component's digits, in odometer order; the
      // first minimum is kept, and a later combination replaces it only if it
      // is lower by more than kEps, the rule of the descent: cur is summed
      // over up to maxCombos incremental updates, and a round-off difference
      // must not decide an exact tie (all null first, so a tie goes to null).
      // Leaves the clusters at the optimum.
      void Exhaustive(const std::vector<int> &K, std::vector<int> &h) {
        const int n = (int)K.size();
        std::vector<int> d(n, 0), best(n, 0);
        double xbest = cur;                          // all null
        for (;;) {
          int j = 0;
          for (; j < n; ++j) {
            const int q = K[j], a = host(q, d[j]);
            if (++d[j] > (int)opt[q].size()) d[j] = 0;
            move(q, a, host(q, d[j]));
            if (d[j] != 0) break;                    // no carry
          }
          if (j == n) break;                         // wrapped around to all null
          if (cur < xbest - kEps) { xbest = cur; best = d; }
        }
        for (int j = 0; j < n; ++j) { h[K[j]] = host(K[j], best[j]); move(K[j], kNull, h[K[j]]); }
      }

      // One coordinate sweep; true if anything moved.
      bool Sweep(const std::vector<int> &K, std::vector<int> &d) {
        bool changed = false;
        for (size_t j = 0; j < K.size(); ++j) {
          const int q = K[j], a = host(q, d[j]);
          move(q, a, kNull);
          const double now = single(q, a);
          int bj = 0; double bv = 0.;                // null first, then ascending vertex
          for (int o = 1; o <= (int)opt[q].size(); ++o) {
            const double x = single(q, opt[q][o-1]);
            if (x < bv) { bv = x; bj = o; }
          }
          if (bv < now - kEps) { d[j] = bj; changed = true; }
          move(q, kNull, host(q, d[j]));
        }
        return changed;
      }

      // Every pair sharing a vertex, jointly over the product of its options.
      bool PairPass(const std::vector<int> &K, std::vector<int> &d) {
        bool changed = false;
        std::vector<double> sa, sb;
        for (size_t j1 = 0; j1 < K.size(); ++j1)
          for (size_t j2 = j1 + 1; j2 < K.size(); ++j2) {
            const int qa = K[j1], qb = K[j2];
            const std::vector<int> &oa = opt[qa], &ob = opt[qb];
            bool share = false;                      // sorted lists: merge walk
            for (size_t x = 0, y = 0; x < oa.size() && y < ob.size() && !share; ) {
              if (oa[x] == ob[y]) share = true; else if (oa[x] < ob[y]) ++x; else ++y;
            }
            if (!share) continue;
            const int ha = host(qa, d[j1]), hb = host(qb, d[j2]);
            move(qa, ha, kNull); move(qb, hb, kNull);
            const int na = (int)oa.size(), nb = (int)ob.size();
            sa.assign(na + 1, 0.); sb.assign(nb + 1, 0.);
            for (int o = 1; o <= na; ++o) sa[o] = single(qa, oa[o-1]);
            for (int o = 1; o <= nb; ++o) sb[o] = single(qb, ob[o-1]);
            const double pp = 2.*(ks[qa].px*ks[qb].px + ks[qa].py*ks[qb].py);
            auto joint = [&](int ia, int ib) {
              double x = sa[ia] + sb[ib];
              if (ia && ib && oa[ia-1] == ob[ib-1]) x += pp*isg[oa[ia-1]];
              return x;
            };
            const double now = joint(d[j1], d[j2]);
            int ba = 0, bb = 0; double bv = joint(0, 0);
            for (int ia = 0; ia <= na; ++ia)
              for (int ib = 0; ib <= nb; ++ib) {
                const double x = joint(ia, ib);
                if (x < bv) { bv = x; ba = ia; bb = ib; }
              }
            if (bv < now - kEps) { d[j1] = ba; d[j2] = bb; changed = true; }
            move(qa, kNull, host(qa, d[j1])); move(qb, kNull, host(qb, d[j2]));
          }
        return changed;
      }

      // Start at the best single move of every cluster against M_v alone,
      // then sweeps and pair passes as described at the top.
      int Descent(const std::vector<int> &K, std::vector<int> &h) {
        const int n = (int)K.size();
        std::vector<int> d(n, 0);
        for (int j = 0; j < n; ++j) {                // all null here: residuals are M_v
          const int q = K[j]; double bv = 0.;
          for (int o = 1; o <= (int)opt[q].size(); ++o) {
            const double x = single(q, opt[q][o-1]);
            if (x < bv) { bv = x; d[j] = o; }
          }
        }
        for (int j = 0; j < n; ++j) move(K[j], kNull, host(K[j], d[j]));
        int nsw = 0;
        while (nsw < c.maxSweeps) {
          const bool moved = Sweep(K, d); ++nsw;
          if (moved) continue;
          if (nsw >= c.maxSweeps) break;
          const bool paired = PairPass(K, d); ++nsw;
          if (!paired) break;
        }
        for (int j = 0; j < n; ++j) h[K[j]] = host(K[j], d[j]);
        return nsw;
      }
    };
  } // namespace detail

  // The assignment minimising X: h[k] = vertex index or kNull.
  inline std::vector<int> Assign(const std::vector<Vertex> &vs, const std::vector<Cluster> &ks,
                                 const Config &c, Stats *st = 0)
  {
    const int nv = (int)vs.size(), nk = (int)ks.size();
    std::vector<int> h(nk, kNull);
    Stats s;
    if (nk > 0 && nv > 0) {
      detail::Solver sv(vs, ks, c);
      detail::UnionFind uf(nk + nv);
      for (int k = 0; k < nk; ++k)
        for (size_t j = 0; j < sv.opt[k].size(); ++j) uf.unite(k, nk + sv.opt[k][j]);
      // components in the order of their lowest cluster, clusters ascending;
      // a cluster with no option at all (every candidate out of range) stays null
      std::vector<int> compOf(nk + nv, -1);
      std::vector< std::vector<int> > comps;
      for (int k = 0; k < nk; ++k) {
        if (sv.opt[k].empty()) continue;
        const int r = uf.find(k);
        if (compOf[r] < 0) { compOf[r] = (int)comps.size(); comps.push_back(std::vector<int>()); }
        comps[compOf[r]].push_back(k);
      }
      for (size_t ic = 0; ic < comps.size(); ++ic) {
        const std::vector<int> &K = comps[ic];
        double combos = 1;
        for (size_t j = 0; j < K.size() && combos <= c.maxCombos; ++j) combos *= double(sv.opt[K[j]].size() + 1);
        if (combos <= c.maxCombos) { sv.Exhaustive(K, h); ++s.exact; }
        else s.sweeps += sv.Descent(K, h);
        ++s.components;
        s.largest = std::max(s.largest, (int)K.size());
      }
    }
    s.objective = Objective(vs, ks, h, c);
    if (st) *st = s;
    return h;
  }

  // The forward clusters of one crossing: anti-kT R (ptmin 0) of the candidates
  // with |eta| > etaFwd and wFwd > 0 at wFwd p_i (mass scaled likewise), so that
  // each of them is in exactly one cluster, with the evidence of both kinds
  // (T_ku of the tracker tracks in the cone, C_ku of the charged constituents)
  // and the candidates candMode makes of it.  vtxIndex maps a vertex key
  // (vertexRef) to its index in the usable list, -1 for an unusable vertex;
  // tracks and constituents on no usable vertex are no evidence.  members, if
  // given, gets the PF indices of every cluster, ascending.  Clusters are in
  // pT order.
  inline std::vector<Cluster> BuildClusters(const puppi::Event &pf, const std::vector<float> &wFwd,
                                            const std::vector<int> &vtxIndex, const Config &c,
                                            std::vector< std::vector<int> > *members = 0)
  {
    std::vector<Cluster> out;
    if (members) members->clear();
    const int np = (int)pf.size();
    std::vector<tiledjet::PseudoJet> in; std::vector<int> idx;
    for (int i = 0; i < np; ++i) {
      if (std::fabs(pf.eta[i]) <= c.etaFwd || !(wFwd[i] > 0)) continue;
      const double w = wFwd[i];
      in.push_back(tiledjet::PseudoJet::PtEtaPhiM(w*pf.pt[i], pf.eta[i], pf.phi[i], w*pf.mass[i], (int)in.size()));
      idx.push_back(i);
    }
    if (in.empty()) return out;
    tiledjet::Config tc; tc.R = c.R; tc.p = -1; tc.ptmin = 0;
    const std::vector<tiledjet::Jet> js = tiledjet::cluster(in, tc);
    const int nk = (int)js.size();
    int nv = 0;
    for (size_t u = 0; u < vtxIndex.size(); ++u) nv = std::max(nv, vtxIndex[u] + 1);
    out.resize(nk);
    if (members) members->resize(nk);
    // C_ku: the charged constituents of every cluster, by usable vertex, raw pT
    std::vector<double> C((size_t)nk*std::max(nv, 1), 0.);
    for (int k = 0; k < nk; ++k) {
      Cluster &q = out[k];
      q.px = js[k].p.px; q.py = js[k].p.py; q.pt = js[k].p.pt();
      q.eta = js[k].p.eta(); q.phi = js[k].p.phi_std();
      if (members) {
        std::vector<int> &m = (*members)[k];
        for (size_t j = 0; j < js[k].constituents.size(); ++j) m.push_back(idx[js[k].constituents[j]]);
        std::sort(m.begin(), m.end());
      }
      for (size_t j = 0; j < js[k].constituents.size() && nv > 0; ++j) {
        const int i = idx[js[k].constituents[j]];
        if (!pf.charged[i]) continue;
        const int key = pf.vref[i];
        if (key < 0 || key >= (int)vtxIndex.size() || vtxIndex[key] < 0) continue;
        C[(size_t)k*nv + vtxIndex[key]] += pf.pt[i];
      }
    }
    // T_ku: the tracker-edge tracks of every usable vertex inside each cone
    std::vector<double> T((size_t)nk*std::max(nv, 1), 0.);
    const double R2 = c.R*c.R;
    for (int i = 0; i < np && nv > 0; ++i) {
      if (!pf.charged[i] || std::fabs(pf.eta[i]) > c.etaFwd) continue;
      const int key = pf.vref[i];
      if (key < 0 || key >= (int)vtxIndex.size() || vtxIndex[key] < 0) continue;
      const int u = vtxIndex[key];
      for (int k = 0; k < nk; ++k) {
        const double de = pf.eta[i] - out[k].eta;
        if (de*de >= R2) continue;
        const double df = puppi::dphi(pf.phi[i], out[k].phi);
        if (de*de + df*df < R2) T[(size_t)k*nv + u] += pf.pt[i];
      }
    }
    for (int k = 0; k < nk; ++k)
      for (int u = 0; u < nv; ++u) {
        Cluster &q = out[k];
        const double t = T[(size_t)k*nv + u], h = C[(size_t)k*nv + u];
        const bool bt = t >= c.trkMinPt && t > 0, bh = h >= c.chMinPt && h > 0;
        if (bt) { q.trkVtx.push_back(u); q.trkPt.push_back(t); }
        if (h > 0) { q.chVtx.push_back(u); q.chPt.push_back(h); }
        if (c.candMode == 1 ? bh : (c.candMode == 2 ? (bt || bh) : bt)) {
          q.cand.push_back(u); q.candTrkPt.push_back(t); q.candChPt.push_back(h); }
      }
    return out;
  }

} // namespace jva
#endif
