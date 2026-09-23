#ifndef tiledjet_h
#define tiledjet_h
// genseed - Copyright 2026 Mikko Voutilainen - Apache License 2.0
//
// tiledjet.h - sequential-recombination jet finding (anti-kT by default; kT
// and Cambridge/Aachen through the exponent p) with the tiled strategy, lazy
// nearest-neighbour updates and a min-heap.  A clean-room implementation:
// written from the published algorithm descriptions (Cacciari and Salam,
// hep-ph/0512210, Sec. 3-4; Cacciari, Salam and Soyez, arXiv:0802.1189), not
// from any existing code.  No external dependency.
//
// WHY.  genseed clusters every event several times per vertex hypothesis - the
// PF candidates plus thousands of soft ghosts, 9000 particles and more - and it
// has to run where FastJet is not installed.  The textbook algorithm is N^3,
// the N^2 nearest-neighbour cache in cluster_n2 below is 0.3 s per event at this
// size, and neither is usable for 25000 events x 50 vertices.  The tiled
// strategy makes each step cost the occupancy of a 3x3 block of tiles instead
// of N, and the clustering ends up close to N log N in practice.  Keeping the
// implementation in one header, GPL-free, is what lets the package be built
// anywhere with a C++17 compiler; FastJet stays an optional backend and the
// benchmark (bench_tiledjet.C) says which one the driver should prefer.
//
// THE ALGORITHM.  Every particle i has kt_i^2p (p = -1 anti-kT, 0 C/A, 1 kT)
// and distances d_ij = min(kt_i^2p, kt_j^2p) dR_ij^2 / R^2 and d_iB = kt_i^2p,
// dR^2 = dy^2 + dphi^2 in rapidity and azimuth.  Repeatedly take the smallest:
// d_ij merges i and j by four-vector addition (E-scheme), d_iB declares i a
// jet.  Two observations make it fast.  (1) The smallest d_ij over all pairs
// is always attained by some i together with its GEOMETRIC nearest neighbour:
// if (i,j) is the minimum pair and kt_i^2p <= kt_j^2p, a particle k closer to
// i than j would give d_ik <= kt_i^2p dR_ik^2/R^2 < d_ij.  Likewise the pair
// can only win over d_iB if dR_ij < R.  So each particle needs to know only
// its nearest neighbour within R, and the next step is the minimum over i of
// key_i = min(kt_i^2p, kt_NN^2p) dR_NN^2/R^2, or kt_i^2p if it has no
// neighbour within R.  (2) With the (y, phi) plane cut into tiles of side >= R,
// that neighbour lies in the 3x3 block of tiles around the particle, so a
// nearest-neighbour search costs the block occupancy rather than N.  After a
// merge or a jet declaration only two kinds of particles change: those whose
// neighbour WAS one of the removed pair (they rescan their 3x3 block in full),
// and those that find the new merged particle closer than their current
// neighbour (one distance each, over the merged particle's 3x3 block).
// Everybody else keeps its neighbour untouched - the lazy update.  A binary
// min-heap over the particle keys, with a position table so that a changed
// key is re-sifted in log N, picks the next step.
//
// DETAILS THAT MATTER.
// - Rapidity, not pseudorapidity: y = sign(pz) ln((E+|pz|)/sqrt(pt^2+m^2)),
//   which is exact for massive particles at zero pT and does not cancel
//   catastrophically at |y| ~ 10.  A massless particle exactly along the beam
//   (E == |pz|, pt = 0) gets y = +-kMaxRap; it never lies within R of anything
//   and ends as a zero-pT jet, kept only if ptmin <= 0.
// - Phi is kept in [0, 2pi) and the difference is folded to [0, pi]; the tile
//   grid is periodic in phi, so the 3x3 blocks wrap.
// - The rapidity range of the tiling is that of the data, clipped to
//   +-kRapClip; particles outside fall into the edge tiles, which is safe
//   because tile membership is only an acceleration structure: a particle
//   beyond the last edge can be within R only of particles in the last two
//   rows.  The number of tiles is also capped at a few per particle so that a
//   handful of far-forward particles does not make the setup O(range/R).
// - Tile boundaries are computed in two ways that round differently (the
//   tile index of a particle, the edge distance used to skip a tile), so the
//   tile side is R + kEdgeEps and every edge distance is shortened by
//   kEdgeEps = 1e-12: a skip is then never wrong, not even for a particle
//   within an ulp of an edge.  Any R > 0 works; for R > pi the phi grid is a
//   single tile and the block is the own row and the two next to it.
// - Zero-pT input with p < 0 would give kt^2p = inf and inf * 0 = NaN, so the
//   pT^2 is floored at 1/kInfKt.  Ghosts at pT 1e-9 give kt^-2 = 1e18, far
//   from any overflow, and their keys sort after all real d_ij, as they must.
// - Recombination is the plain sum of the four-vectors; the merged particle
//   is a NEW slot (indices never reused), so 'my neighbour was i' is an
//   unambiguous test, and constituent lists are concatenated in O(1) through
//   per-slot head/tail pointers into a single next-array over the inputs.
// - Ties are broken by rank: slot 0..N-1 is the input index, slot N+n the
//   n-th merged object, and an exactly equal distance or key goes to the
//   lower slot (a nearest neighbour is the closest, then the lowest slot;
//   the heap orders by key, then slot).  Real PF input HAS exact ties - HF
//   towers and ECAL crystals are a grid, two towers at the same float eta
//   and mirrored phi are equidistant to the bit - and without the rule the
//   result would depend on the tiling.  cluster_n2 applies the same rule, so
//   the two agree bit for bit; FastJet has its own, see bench_tiledjet.C.
// - Output: jets with pT >= ptmin, sorted by pT descending, with their
//   constituents as indices into the input vector (in recombination order).
// - Ghosts: scale the MASS with the pT (PtEtaPhiM(pt*1e-9, eta, phi, m*1e-9)).
//   A ghost with pT 1e-9 and the real mass has pz << m and sits at rapidity
//   0 whatever its eta; a whole event of them lands in one tile row and the
//   block scans go back to N^2.  The clustering stays correct, only slow.
//
// EXACTNESS.  With identical input four-vectors the merge sequence is the same
// as in any other correct implementation unless two distances tie, so the
// jets agree with the N^2 reference cluster_n2 (same header, same tie rule,
// and tiles that only ever skip what provably cannot win - see kEdgeEps)
// bit for bit, and with FastJet jet by jet except where real PF input ties:
// the calorimeter grid puts HF towers and ECAL photons at exactly equal float
// eta or phi, two d_ij then agree exactly or to 1e-13 - the last bit of a
// rapidity, which differs between any two formulas - and which pair merges
// first is a convention, not physics.  On 200 events that is 505 jets in 27
// events for the unweighted PF candidates, and none with PUPPI weights, with
// ghosts, or with a 1e-7 jitter of the same candidates (bench_tiledjet.C
// shows all four).  The verdict of the last run on this machine is recorded
// in kAtLeastAsFastAsFastJet / kBenchmark below.  runGenSeed.C chooses the
// default backend from kAtLeastAsFastAsFastJet, which it reads from this
// file AS TEXT (the header need not be interpretable by cling): FastJet is
// used only if it is found and the verdict is false.  kPreferFastJet is
// derived from the verdict for C++ users; editing it changes nothing.
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstddef>

namespace tiledjet {

  // ---- benchmark verdict (written by hand from bench_tiledjet.C output) ---
  // Every input (PF x PUPPI N ~ 120, PF all N ~ 2100, PF + gen ghosts
  // N ~ 10000) within 1.0x of FastJet's Best strategy on this machine -> true.
  // THE switch: runGenSeed.C parses the next line as text, so keep it in
  // this form, '= true;' or '= false;' with nothing after it (no comment).
  inline constexpr bool kAtLeastAsFastAsFastJet = true;
  // derived, informational only - do not edit, runGenSeed.C never reads it
  inline constexpr bool kPreferFastJet = !kAtLeastAsFastAsFastJet;
  inline const char *const kBenchmark =
    "bench_tiledjet.C 2026-09-23, 200 events of MC24NanoV15_PU_IT_OOT/NANOAODSIM_1.root, "
    "anti-kT R=0.4 ptmin=1, Apple M2 (4P+4E) 24 GB, macOS 13.5 (Darwin 22.6), Apple clang 15, "
    "ACLiC -O3, FastJet 3.5.1 (Best): ms/event mean+-rms tiledjet / FastJet / N^2 ref: "
    "PF x PUPPI w>0 N=123: 0.061+-0.024 / 0.070+-0.034 (N2Tiled) / 0.14, ratio 0.87; "
    "PF all N=2094: 1.20+-0.30 / 1.61+-0.39 (N2MinHeapTiled) / 12.9, ratio 0.74; "
    "PF all jittered 1e-7 N=2094: 1.18+-0.27 / 1.60+-0.38 (N2MinHeapTiled) / 12.9, ratio 0.74; "
    "PF all + gen ghosts N=10149: 6.54+-2.19 / 9.25+-3.22 (N2MHTLazy9) / 340, ratio 0.71 "
    "(machine loaded to ~4.5 by other jobs; the other of two runs that day gave ratios "
    "0.79 / 0.74 / 0.73 / 0.71). Exactness: tiledjet == cluster_n2 on every jet of every input; "
    "== FastJet on every jet of the PUPPI, ghost and jittered inputs; 505 jets in 27 "
    "unweighted-PF events differ by grid-tie resolution";

  // ---- constants ----------------------------------------------------------
  constexpr double kPi      = 3.14159265358979323846;
  constexpr double kTwoPi   = 2.0*kPi;
  constexpr double kMaxRap  = 1e5;    // rapidity given to E == |pz|, pt == 0
  constexpr double kRapClip = 50.0;   // the tiling covers at most [-50, 50]
  constexpr double kInfKt   = 1e200;  // kt^2p of a zero-pT particle for p < 0
  constexpr double kTileScale = 1.0;  // tile side = kTileScale * R (>= R) + kEdgeEps
  constexpr int    kMaxTilesPerParticle = 1;   // caps the tile count at small N
  constexpr double kEdgeEps = 1e-12;  // tile-edge safety margin, see edges() and setup()

  // ---- the four-vector ----------------------------------------------------
  struct PseudoJet {
    double px = 0, py = 0, pz = 0, E = 0;
    int user_index = -1;
    PseudoJet() = default;
    PseudoJet(double x, double y, double z, double e, int idx = -1)
      : px(x), py(y), pz(z), E(e), user_index(idx) {}
    // From (pT, eta, phi, m): the conversion every reader needs.  A negative
    // m^2 (float round-off) is treated as zero.
    static PseudoJet PtEtaPhiM(double pt, double eta, double phi, double m,
                               int idx = -1) {
      const double x = pt*std::cos(phi), y = pt*std::sin(phi), z = pt*std::sinh(eta);
      const double m2 = m > 0 ? m*m : 0.0;
      return PseudoJet(x, y, z, std::sqrt(x*x + y*y + z*z + m2), idx);
    }
    double pt2() const { return px*px + py*py; }
    double pt()  const { return std::sqrt(pt2()); }
    double modp2() const { return px*px + py*py + pz*pz; }
    double m2()  const { return E*E - modp2(); }
    double m()   const { const double x = m2(); return x > 0 ? std::sqrt(x) : -std::sqrt(-x); }
    double e()   const { return E; }
    // phi in [0, 2pi)
    double phi() const {
      double f = std::atan2(py, px);
      if (f < 0) f += kTwoPi;
      if (f >= kTwoPi) f = 0;
      return f;
    }
    // phi in (-pi, pi]
    double phi_std() const { return std::atan2(py, px); }
    // rapidity, exact for massive zero-pT particles, guarded at E == |pz|
    double rap() const {
      const double apz = std::fabs(pz), pt2v = pt2();
      double m2v = E*E - pt2v - pz*pz;
      if (m2v < 0) m2v = 0;
      const double den = pt2v + m2v;                 // (E+|pz|)(E-|pz|)
      if (den <= 0) return E + apz <= 0 ? 0.0 : (pz < 0 ? -kMaxRap : kMaxRap);
      const double y = std::log((E + apz)/std::sqrt(den));
      return pz < 0 ? -y : y;
    }
    // pseudorapidity, +-kMaxRap along the beam
    double eta() const {
      const double pt1 = pt();
      if (pt1 <= 0) return pz < 0 ? -kMaxRap : (pz > 0 ? kMaxRap : 0.0);
      return std::asinh(pz/pt1);
    }
    PseudoJet &operator+=(const PseudoJet &o) {
      px += o.px; py += o.py; pz += o.pz; E += o.E; return *this;
    }
    friend PseudoJet operator+(PseudoJet a, const PseudoJet &b) { a += b; a.user_index = -1; return a; }
  };

  struct Config {
    double R = 0.4;
    double p = -1;     // -1 anti-kT, 0 Cambridge/Aachen, 1 kT
    double ptmin = 0;  // jets with pT >= ptmin are returned
  };

  struct Jet {
    PseudoJet p;
    std::vector<int> constituents;   // indices into the input vector
  };

  inline bool ByPt(const Jet &a, const Jet &b) { return a.p.pt2() > b.p.pt2(); }

  // kt^2p with the zero-pT guard
  inline double Kt2P(double pt2, double p) {
    if (p == -1) return pt2 > 1.0/kInfKt ? 1.0/pt2 : kInfKt;
    if (p ==  1) return pt2;
    if (p ==  0) return 1.0;
    if (pt2 <= 1.0/kInfKt) return p < 0 ? kInfKt : 0.0;
    return std::pow(pt2, p);
  }

  namespace detail {

    // Binary min-heap of slot ids over an external key array, with a position
    // table so any slot's key can be changed and re-sifted in O(log n).
    struct MinHeap {
      struct Node { double key; int slot; };
      std::vector<Node> h; std::vector<int> pos;
      void init(int cap) { h.clear(); h.reserve(cap); pos.assign(cap, -1); }
      int  size() const { return (int)h.size(); }
      int  top()  const { return h[0].slot; }
      // (key, slot) lexicographic: exact ties go to the lower slot
      static bool less(const Node &a, const Node &b) { return a.key < b.key || (a.key == b.key && a.slot < b.slot); }
      void up(int idx) {
        const Node n = h[idx];
        while (idx > 0) {
          const int par = (idx - 1) >> 1;
          if (!less(n, h[par])) break;
          h[idx] = h[par]; pos[h[idx].slot] = idx; idx = par;
        }
        h[idx] = n; pos[n.slot] = idx;
      }
      void down(int idx) {
        const int sz = (int)h.size(); const Node n = h[idx];
        for (;;) {
          int c = 2*idx + 1;
          if (c >= sz) break;
          if (c + 1 < sz && less(h[c+1], h[c])) ++c;
          if (!less(h[c], n)) break;
          h[idx] = h[c]; pos[h[idx].slot] = idx; idx = c;
        }
        h[idx] = n; pos[n.slot] = idx;
      }
      void push(int s, double k) { h.push_back(Node{k, s}); pos[s] = (int)h.size() - 1; up(pos[s]); }
      void build() { for (int i = (int)h.size()/2 - 1; i >= 0; --i) down(i); }
      void erase(int s) {
        const int idx = pos[s]; pos[s] = -1;
        const Node last = h.back(); h.pop_back();
        if (idx == (int)h.size()) return;
        h[idx] = last; pos[last.slot] = idx;
        up(idx); down(pos[last.slot]);
      }
      void pop() { erase(h[0].slot); }
      void update(int s, double k) { h[pos[s]].key = k; up(pos[s]); down(pos[s]); }
      // k takes over j's place: one sift instead of an erase and a push
      void replace(int j, int k, double key) { const int idx = pos[j]; pos[j] = -1; h[idx].slot = k; pos[k] = idx; update(k, key); }
    };

    // One record per live particle, stored CONTIGUOUSLY per tile: a block
    // scan is a linear walk over 32-byte records, not a pointer chase.  The
    // slot id is the particle's identity (heap key, momentum, kt^2p, lists).
    struct Rec { double rap, phi, dnn; int slot, nn; };

    class Engine {
    public:
      Engine() = default;
      // The engine is reused between calls (a thread_local in cluster()), so
      // the vectors keep their capacity and an event costs no allocation.
      void reset(const std::vector<PseudoJet> &in, const Config &c) {
        cfg = c; N = (int)in.size(); R2 = c.R*c.R; invR2 = 1.0/R2; setup(in);
      }

      std::vector<Jet> run() {
        std::vector<int> jets; jets.reserve(N);
        std::vector<int> changed; changed.reserve(128);
        while (heap.size()) {
          const int i = heap.top();
          const int j = pool[rpos[i]].nn;
          // whoever had i (or j) as neighbour must look again: the reverse
          // lists know exactly who, no block has to be scanned for it
          // (j may be in i's list and i in j's: both die, neither is rescanned)
          for (int m = rhead[i]; m >= 0; m = rnext[m]) if (m != j) changed.push_back(m);
          if (j >= 0) for (int m = rhead[j]; m >= 0; m = rnext[m]) if (m != i) changed.push_back(m);
          if (j < 0) {
            // ---- i becomes a jet -------------------------------------------
            heap.pop(); remove(i); jets.push_back(i);
          } else {
            // ---- merge i and j into the new slot k --------------------------
            // (the dead leave the reverse lists of their own neighbours)
            setNN(i, j, -1); setNN(j, pool[rpos[j]].nn, -1);
            heap.pop(); remove(i); remove(j);
            const int k = nslot++;
            mom[k] = mom[i]; mom[k] += mom[j];
            setkin(k); insert(k);
            chead[k] = chead[i]; cnext[ctail[i]] = chead[j]; ctail[k] = ctail[j];
            rhead[k] = -1;
            // k's block: find k's neighbour and offer k to everybody there.
            // Own tile first; a neighbouring tile is skipped when k is
            // farther from its edge than from its current neighbour AND
            // than any occupant of that tile is from its own (tmax).
            // (k's best is kept in locals: written through rk it would be
            // reloaded every iteration, the compiler having to assume that
            // rk and *r may alias)
            const Rec rk = pool[rpos[k]];
            double bestd = R2; int bestn = -1;
            const int tk = tile[k], *L = &nb[9*tk];
            double mxself = 0;
            {
              Rec *r = &pool[tbeg[tk]], *e = r + tcnt[tk];
              for (; r != e; ++r) {
                if (r->slot == k) continue;
                const double d = dist2(rk, *r);
                if (d < bestd || (d == bestd && r->slot < bestn)) { bestd = d; bestn = r->slot; }
                if (d < r->dnn) {   // k is the newest slot: on a tie it never replaces
                  const int old = r->nn;
                  r->dnn = d; setNN(r->slot, old, k); r->nn = k;
                  if (old != i && old != j) changed.push_back(r->slot);
                }
                if (r->dnn > mxself) mxself = r->dnn;
              }
            }
            double edge[9]; edges(rk, tk, edge);
            for (int o = 0; o < 8; ++o) {
              const int q = kOrder[o], t = L[q];
              if (t < 0 || t == tk || (edge[q] > bestd && edge[q] >= tmax[t])) continue;   // (t == tk: see findNN)
              double mx = 0;
              Rec *r = &pool[tbeg[t]], *e = r + tcnt[t];
              for (; r != e; ++r) {
                const double d = dist2(rk, *r);
                if (d < bestd || (d == bestd && r->slot < bestn)) { bestd = d; bestn = r->slot; }
                if (d < r->dnn) {
                  const int old = r->nn;
                  r->dnn = d; setNN(r->slot, old, k); r->nn = k;
                  if (old != i && old != j) changed.push_back(r->slot);
                }
                if (r->dnn > mx) mx = r->dnn;
              }
              tmax[t] = mx;
            }
            pool[rpos[k]].dnn = bestd; pool[rpos[k]].nn = bestn;
            tmax[tk] = mxself > bestd ? mxself : bestd;
            setNN(k, -1, bestn);
            heap.replace(j, k, keyof(k));
          }
          for (size_t q = 0; q < changed.size(); ++q) {
            const int m = changed[q];
            Rec &rm = pool[rpos[m]];
            if (rm.nn == i || rm.nn == j) {
              const int old = rm.nn; findNN(m, rm); setNN(m, old, rm.nn);
              if (rm.dnn > tmax[tile[m]]) tmax[tile[m]] = rm.dnn;
            }
            heap.update(m, keyof(m));
          }
          changed.clear();
        }
        // ---- output ---------------------------------------------------------
        std::vector<Jet> out; out.reserve(jets.size());
        for (size_t q = 0; q < jets.size(); ++q) {
          const int s = jets[q];
          if (mom[s].pt() < cfg.ptmin) continue;
          Jet J; J.p = mom[s]; J.p.user_index = -1;
          for (int c = chead[s]; c >= 0; c = cnext[c]) J.constituents.push_back(c);
          out.push_back(std::move(J));
        }
        std::stable_sort(out.begin(), out.end(), ByPt);
        return out;
      }

    private:
      Config cfg;
      int N = 0, nslot = 0;
      double R2 = 0, invR2 = 0;
      // tiling: tile t = ir*nphi + ip owns pool[tbeg[t] .. tbeg[t]+tcnt[t])
      int nrap = 1, nphi = 1, ntiles = 1;
      bool skipOK = false;
      double rap0 = 0, drap = 1, drapInv = 0, dphi = kTwoPi, dphiInv = 0;
      std::vector<int> tbeg, tcnt, tcap, nb;   // nb: 9 per tile, (dr+1)*3+(dp+1), -1 absent
      std::vector<double> tmax;                // per tile: upper bound of dnn over its records
      static constexpr int kOrder[8] = {1, 3, 5, 7, 0, 2, 6, 8};   // edge-adjacent first, corners last
      std::vector<Rec> pool;
      // per slot (2N): momentum, kt^2p, record position, tile,
      // reverse neighbour lists (who has me as NN), constituent chains
      std::vector<PseudoJet> mom;
      std::vector<double> kt2p;
      std::vector<int> rpos, tile, rhead, rnext, rprev, chead, ctail;
      std::vector<int> cnext;        // N: constituent chain over the inputs
      MinHeap heap;

      // branch-free: fabs and the fold compile to selects, and the loops
      // over a tile's records then run without mispredictions
      static double dist2(const Rec &a, const Rec &b) {
        const double dr = a.rap - b.rap;
        const double da = std::fabs(a.phi - b.phi);
        const double dp = da > kPi ? kTwoPi - da : da;
        return dr*dr + dp*dp;
      }
      int tileOf(double r, double f) const {
        int ir = (int)((r - rap0)*drapInv);
        if (ir < 0) ir = 0; else if (ir >= nrap) ir = nrap - 1;
        int ip = (int)(f*dphiInv);
        if (ip >= nphi) ip = nphi - 1; else if (ip < 0) ip = 0;
        return ir*nphi + ip;
      }
      void setkin(int s) {
        const PseudoJet &q = mom[s];
        const double r = q.rap(), f = q.phi();
        kt2p[s] = Kt2P(q.pt2(), cfg.p);
        tile[s] = tileOf(r, f);
        rec.rap = r; rec.phi = f; rec.slot = s;   // staged, insert() copies it
      }
      Rec rec;
      // tile segments: append; on overflow move the tile to the end of the pool
      void insert(int s) {
        const int t = tile[s];
        if (tcnt[t] == tcap[t]) {
          const int nc = 2*tcap[t] + 4, nbg = (int)pool.size();
          pool.resize(nbg + nc);
          for (int q = 0; q < tcnt[t]; ++q) { pool[nbg + q] = pool[tbeg[t] + q]; rpos[pool[nbg + q].slot] = nbg + q; }
          tbeg[t] = nbg; tcap[t] = nc;
        }
        const int idx = tbeg[t] + tcnt[t]++;
        pool[idx] = rec; rpos[s] = idx;
      }
      void remove(int s) {
        const int t = tile[s], idx = rpos[s], last = tbeg[t] + --tcnt[t];
        if (idx != last) { pool[idx] = pool[last]; rpos[pool[idx].slot] = idx; }
        rpos[s] = -1;
      }
      // reverse lists: m's neighbour changes from a to b
      void setNN(int m, int a, int b) {
        if (a == b) return;
        if (a >= 0) {
          if (rprev[m] >= 0) rnext[rprev[m]] = rnext[m]; else rhead[a] = rnext[m];
          if (rnext[m] >= 0) rprev[rnext[m]] = rprev[m];
        }
        if (b >= 0) {
          rnext[m] = rhead[b]; rprev[m] = -1;
          if (rhead[b] >= 0) rprev[rhead[b]] = m;
          rhead[b] = m;
        }
      }
      double keyof(int s) const {
        const Rec &r = pool[rpos[s]];
        if (r.nn < 0) return kt2p[s];
        const double k = kt2p[s] < kt2p[r.nn] ? kt2p[s] : kt2p[r.nn];
        return k*r.dnn*invR2;
      }
      void scanTile(int t, const Rec &rm, double &bestd, int &bestn) const {
        const Rec *r = &pool[tbeg[t]], *e = r + tcnt[t];
        double bd = bestd; int bn = bestn;
        for (; r != e; ++r) {
          const double d = dist2(rm, *r);
          if (d < bd || (d == bd && r->slot < bn)) { bd = d; bn = r->slot; }
        }
        bestd = bd; bestn = bn;
      }
      // Squared distance from a particle in tile t to the edge shared with
      // each of the 8 neighbouring tiles (0 for its own, and for everything
      // when the phi grid is too coarse to tell -1 from +1).  Particles
      // clipped into an edge tile get 0 too: never skipped, always correct.
      // Every linear distance is shortened by kEdgeEps before it is squared,
      // so the bound is a LOWER bound for sure: tileOf() puts a particle in
      // a tile with (r - rap0)*drapInv and f*dphiInv, the edges here are
      // rap0 + ir*drap and ip*dphi, and the two round differently within a
      // few ulps of a boundary (1e-13 at worst for |y| <= kRapClip; the phi
      // wrap, nphi*dphi against 2pi, is the same).  Without the margin a
      // particle that tileOf() had put in the next tile, 1e-16 closer than
      // the edge computed here, was skipped, and a near-tie merged in another
      // order than in cluster_n2.  The margin changes a decision only when a
      // particle is within ~1e-12 of the edge distance, so it costs nothing
      // measurable, and none of the skips below can then be wrong.
      void edges(const Rec &rm, int t, double *edge) const {
        if (!skipOK) { for (int q = 0; q < 9; ++q) edge[q] = 0; return; }
        const int ir = t/nphi, ip = t - ir*nphi;
        double dl = rm.rap - (rap0 + ir*drap) - kEdgeEps, dh = (rap0 + (ir + 1)*drap) - rm.rap - kEdgeEps;
        double fl = rm.phi - ip*dphi - kEdgeEps, fh = (ip + 1)*dphi - rm.phi - kEdgeEps;
        dl = dl > 0 ? dl*dl : 0; dh = dh > 0 ? dh*dh : 0; fl = fl > 0 ? fl*fl : 0; fh = fh > 0 ? fh*fh : 0;
        edge[0] = dl + fl; edge[1] = dl; edge[2] = dl + fh;
        edge[3] = fl;      edge[4] = 0;  edge[5] = fh;
        edge[6] = dh + fl; edge[7] = dh; edge[8] = dh + fh;
      }
      // Nearest neighbour of m within R: own tile first, then the four
      // edge-adjacent and the four corner tiles, each skipped when m is
      // already closer to its current neighbour than to that tile's edge.
      // Only the own-tile loop skips m itself, so the own tile must never be
      // scanned again from an off-centre slot: setup() guarantees that for
      // every phi grid, and the 'u == t' test (eight integer compares per
      // search, nothing per record) keeps it true whatever the tiling does.
      void findNN(int m, Rec &rec) const {
        const Rec rm = rec;
        double bestd = R2; int bestn = -1;
        const int t = tile[m], *L = &nb[9*t];
        {
          const Rec *r = &pool[tbeg[t]], *e = r + tcnt[t];
          for (; r != e; ++r) {
            if (r->slot == m) continue;
            const double d = dist2(rm, *r);
            if (d < bestd || (d == bestd && r->slot < bestn)) { bestd = d; bestn = r->slot; }
          }
        }
        double edge[9]; edges(rm, t, edge);
        for (int o = 0; o < 8; ++o) {
          const int q = kOrder[o], u = L[q];
          if (u < 0 || u == t || edge[q] > bestd) continue;   // (>: a tie on the edge is still scanned)
          scanTile(u, rm, bestd, bestn);
        }
        rec.dnn = bestd; rec.nn = bestn;
      }

      void setup(const std::vector<PseudoJet> &in) {
        const int cap = 2*N + 1;
        mom.resize(cap); kt2p.resize(cap); rpos.assign(cap, -1); tile.resize(cap);
        rhead.assign(cap, -1); rnext.assign(cap, -1); rprev.assign(cap, -1);
        chead.resize(cap); ctail.resize(cap); cnext.assign(N, -1);
        // rapidities first: they fix the tiling
        std::vector<double> r(N), f(N);
        double rmin = kRapClip, rmax = -kRapClip;
        for (int i = 0; i < N; ++i) {
          r[i] = in[i].rap(); f[i] = in[i].phi();
          const double rc = r[i] < -kRapClip ? -kRapClip : (r[i] > kRapClip ? kRapClip : r[i]);
          if (rc < rmin) rmin = rc;
          if (rc > rmax) rmax = rc;
        }
        if (N == 0) { rmin = rmax = 0; }
        // The side is R plus kEdgeEps for the reason edges() gives: with a
        // side of exactly R (range/R an integer, e.g. the clipped [-50, 50]
        // at R = 0.4 once N is large enough not to cap the tile count) the
        // tileOf() rounding could put two particles closer than R two rows
        // apart, outside each other's 3x3 block.
        const double side = (kTileScale > 1 ? kTileScale : 1.0)*cfg.R + kEdgeEps;
        nphi = std::max(1, (int)std::floor(kTwoPi/side));
        const double range = rmax - rmin;
        nrap = std::max(1, (int)std::floor(range/side));
        const int maxTiles = std::max(nphi, kMaxTilesPerParticle*std::max(N, 16));
        if (nrap*nphi > maxTiles) nrap = std::max(1, maxTiles/nphi);
        ntiles = nrap*nphi;
        rap0 = rmin;
        drap = range > 0 ? range/nrap : 1.0; drapInv = range > 0 ? nrap/range : 0.0;
        dphi = kTwoPi/nphi; dphiInv = nphi/kTwoPi;
        skipOK = nphi >= 3;
        nb.assign(9*ntiles, -1);
        for (int ir = 0; ir < nrap; ++ir)
          for (int ip = 0; ip < nphi; ++ip) {
            // The own tile goes into the centre slot FIRST.  With one phi
            // tile (R > pi) all three phi columns are the own column, and
            // filled in dr/dp order the own tile landed in slot 3, which the
            // neighbour loops scan without skipping the particle itself: it
            // became its own nearest neighbour at d = 0, merged with itself
            // and wrote pool[-1].  Seeded here, the dedup below drops every
            // alias of it, and slot 4 is the only place it can be.
            int *L = &nb[9*(ir*nphi + ip)];
            L[4] = ir*nphi + ip;
            const int jp[3] = {ip == 0 ? nphi - 1 : ip - 1, ip, ip == nphi - 1 ? 0 : ip + 1};
            for (int dr = -1; dr <= 1; ++dr) {
              const int jr = ir + dr;
              if (jr < 0 || jr >= nrap) continue;
              for (int dp = 0; dp < 3; ++dp) {
                const int id = jr*nphi + jp[dp];
                bool dup = false;   // only a phi grid of 1 or 2 tiles can wrap onto itself
                if (!skipOK) for (int q = 0; q < 9; ++q) if (L[q] == id) dup = true;
                if (!dup) L[(dr + 1)*3 + dp] = id;
              }
            }
          }
        // slot = input index; the records are ordered by tile (a counting
        // sort), each tile owning a contiguous segment of the pool with a
        // little slack for merged objects that move in
        tcnt.assign(ntiles, 0); tcap.resize(ntiles); tbeg.resize(ntiles);
        for (int i = 0; i < N; ++i) { tile[i] = tileOf(r[i], f[i]); ++tcnt[tile[i]]; }
        int off = 0;
        for (int t = 0; t < ntiles; ++t) { tbeg[t] = off; tcap[t] = tcnt[t] + 4; off += tcap[t]; tcnt[t] = 0; }
        if ((int)pool.capacity() < off + 4*N) pool.reserve(off + 4*N);
        pool.resize(off);
        for (int i = 0; i < N; ++i) {
          const int t = tile[i], idx = tbeg[t] + tcnt[t]++;
          mom[i] = in[i];
          kt2p[i] = Kt2P(in[i].pt2(), cfg.p);
          chead[i] = ctail[i] = i;
          Rec &q = pool[idx]; q.rap = r[i]; q.phi = f[i]; q.dnn = R2; q.slot = i; q.nn = -1;
          rpos[i] = idx;
        }
        nslot = N;
        // initial neighbours: every particle searches its own block, with
        // the edge-distance skipping of findNN; cheaper than the symmetric
        // pair scan because most particles never leave their own tile
        for (int s = 0; s < N; ++s) findNN(s, pool[rpos[s]]);
        tmax.assign(ntiles, 0);
        for (int t = 0; t < ntiles; ++t)
          for (int q = tbeg[t]; q < tbeg[t] + tcnt[t]; ++q) if (pool[q].dnn > tmax[t]) tmax[t] = pool[q].dnn;
        heap.init(cap);
        for (int s = 0; s < N; ++s) {
          setNN(s, -1, pool[rpos[s]].nn);
          heap.h.push_back(MinHeap::Node{keyof(s), s}); heap.pos[s] = s;
        }
        heap.build();
      }
    };
  } // namespace detail

  // The clusterer.  Returns the jets with pT >= c.ptmin, pT-ordered, with
  // constituents as indices into 'in'.
  inline std::vector<Jet> cluster(const std::vector<PseudoJet> &in, const Config &c = Config()) {
    if (in.empty()) return std::vector<Jet>();
    static thread_local detail::Engine e;
    e.reset(in, c);
    return e.run();
  }

  // Plain N^2 reference: the same algorithm with a nearest-neighbour cache
  // over ALL particles and a linear scan for the minimum.  No tiles, no heap,
  // nothing clever; it exists so that bench_tiledjet.C can prove the fast
  // one right.  About 0.34 s per event at N ~ 10000 (50x the tiled one).
  inline std::vector<Jet> cluster_n2(const std::vector<PseudoJet> &in, const Config &c = Config()) {
    const int N = (int)in.size();
    std::vector<Jet> out;
    if (N == 0) return out;
    const double R2 = c.R*c.R, invR2 = 1.0/R2;
    const int cap = 2*N;
    std::vector<PseudoJet> mom(cap);
    std::vector<double> rap(cap), phi(cap), kt2p(cap), dnn(cap);
    std::vector<int> nn(cap), chead(cap), ctail(cap), cnext(N, -1), alive;
    alive.reserve(N);
    auto setkin = [&](int s) {
      rap[s] = mom[s].rap(); phi[s] = mom[s].phi(); kt2p[s] = Kt2P(mom[s].pt2(), c.p);
    };
    auto dist2 = [&](int a, int b) {
      const double dr = rap[a] - rap[b];
      double dp = std::fabs(phi[a] - phi[b]);
      if (dp > kPi) dp = kTwoPi - dp;
      return dr*dr + dp*dp;
    };
    auto keyof = [&](int s) {
      if (nn[s] < 0) return kt2p[s];
      return std::min(kt2p[s], kt2p[nn[s]])*dnn[s]*invR2;
    };
    // ties: the closest, then the lowest slot - the rule of cluster()
    auto rescan = [&](int m) {
      nn[m] = -1; dnn[m] = R2;
      for (size_t q = 0; q < alive.size(); ++q) {
        const int n = alive[q];
        if (n == m) continue;
        const double d = dist2(m, n);
        if (d < dnn[m] || (d == dnn[m] && n < nn[m])) { dnn[m] = d; nn[m] = n; }
      }
    };
    for (int i = 0; i < N; ++i) {
      mom[i] = in[i]; setkin(i); nn[i] = -1; dnn[i] = R2; chead[i] = ctail[i] = i; alive.push_back(i);
    }
    for (int i = 0; i < N; ++i)
      for (int j = i + 1; j < N; ++j) {
        const double d = dist2(i, j);
        if (d < dnn[i] || (d == dnn[i] && j < nn[i])) { dnn[i] = d; nn[i] = j; }
        if (d < dnn[j] || (d == dnn[j] && i < nn[j])) { dnn[j] = d; nn[j] = i; }
      }
    int nslot = N;
    std::vector<int> jets;
    while (!alive.empty()) {
      int best = 0; double kbest = keyof(alive[0]);
      for (size_t q = 1; q < alive.size(); ++q) {
        const double k = keyof(alive[q]);
        if (k < kbest || (k == kbest && alive[q] < alive[best])) { kbest = k; best = (int)q; }
      }
      const int i = alive[best], j = nn[i];
      int k = -1;
      if (j < 0) {
        jets.push_back(i);
        alive[best] = alive.back(); alive.pop_back();
      } else {
        k = nslot++;
        mom[k] = mom[i]; mom[k] += mom[j]; setkin(k);
        chead[k] = chead[i]; cnext[ctail[i]] = chead[j]; ctail[k] = ctail[j];
        alive[best] = alive.back(); alive.pop_back();
        for (size_t q = 0; q < alive.size(); ++q) if (alive[q] == j) { alive[q] = alive.back(); alive.pop_back(); break; }
        alive.push_back(k);
        rescan(k);
      }
      for (size_t q = 0; q < alive.size(); ++q) {
        const int m = alive[q];
        if (m == k) continue;
        if (nn[m] == i || nn[m] == j) rescan(m);
        else if (k >= 0) { const double d = dist2(m, k); if (d < dnn[m]) { dnn[m] = d; nn[m] = k; } }   // k newest: no tie
      }
    }
    out.reserve(jets.size());
    for (size_t q = 0; q < jets.size(); ++q) {
      const int s = jets[q];
      if (mom[s].pt() < c.ptmin) continue;
      Jet J; J.p = mom[s]; J.p.user_index = -1;
      for (int x = chead[s]; x >= 0; x = cnext[x]) J.constituents.push_back(x);
      out.push_back(std::move(J));
    }
    std::stable_sort(out.begin(), out.end(), ByPt);
    return out;
  }

} // namespace tiledjet
#endif
