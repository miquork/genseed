// genseed - Copyright 2026 Mikko Voutilainen - Apache License 2.0
#ifndef genseed_puppi_h
#define genseed_puppi_h
// puppi.h - PUPPI weights recomputed from the stored particle-flow candidates,
// for an arbitrary reconstructed vertex rather than only for the primary one.
//
// genseed port of lowptjets/puppi.h (v4), without ROOT: the only ROOT call was
// ROOT::Math::chisquared_cdf(x, 1), and the chi2 distribution with one degree
// of freedom is the square of a standard normal, so its cdf is
// P(|N| < sqrt(x)) = erf(sqrt(x/2)), which <cmath> provides.  Verified against
// the v4 header on real events: max |dw| < 1e-6 over every candidate for the
// primary and for vertex 5.  Numerically identical to v4 wherever v4's chi2 is
// finite.  The one place it is not is a forward region with a zero rms (at
// most four candidates, see Prepare): v4 divided by it and got NaN weights,
// this port takes rms = 1 there, as Weights() always did in the tracker.
//
// WHY.  NanoAOD gives one PUPPI weight per candidate, computed against the
// vertex the reconstruction chose as primary.  v1 of this analysis showed that
// inside the tracker that single choice is what decides whether a jet exists
// at all: 90% of the generated jets from the chosen vertex are reconstructed
// against 2-9% of the jets from any other.  To ask what the other forty-four
// interactions would have looked like, the weights have to be recomputed with
// a different vertex playing the role of the primary one, and that is what
// this does.
//
// THE ALGORITHM (Bertolini, Harris, Low, Tran, arXiv:1407.6013), as CMS runs
// it.  For every candidate i,
//
//     alpha_i = log sum_{j != i, dR_ij < R0} ( pT_j / dR_ij )^2
//
// where the sum runs over the charged candidates assigned to the vertex under
// test when i is inside the tracker, and over every candidate when it is not -
// there being no tracks in HF to assign.  Charged candidates in the tracker
// are then given weight 1 if they belong to the vertex and 0 if they do not;
// everything else gets
//
//     chi2_i = (alpha_i - med) |alpha_i - med| / rms^2 ,   w_i = F_chi2,1(chi2_i)
//
// with med and rms taken from the charged pileup candidates in the tracker,
// and from all candidates in each forward region, in both cases only those
// above rmsPtMin = 0.1 GeV.  med is the median in the tracker but the 20%
// quantile in each forward region (medQuantileFwd = 0.2), and rms is not the
// plain rms: it is 0.7 times the RMS about med of the alphas below med only
// (rmsScale, leftRMS; see Config).  Dropping the 0.7 alone would halve every
// chi2.  A neutral candidate whose weighted pT falls below a pileup-dependent
// floor is dropped outright.
//
// WHAT IS AVAILABLE AND WHAT IS NOT.  PFCandV2_vertexRef gives the reco vertex
// of every charged candidate and the number of distinct keys per event is
// exactly PV_npvs, so every vertex is represented.  Only four vertex z
// positions are stored, though (PV_z and three OtherPV_z against 34 vertices),
// so VertexZ() below recovers z for all of them from the median of
// PV_z + PFCandV2_dz over the tracks of each vertex.  That is good to 25 um
// for the primary vertex and 49 um for the next three, checked against the
// stored positions in puppiprobe.C.  Those z positions are what the gen
// interactions are matched to; the weights themselves never need them,
// because the candidate-to-vertex assignment is the stored vertexRef.
//
// FORWARD REGIONS ARE VERTEX-BLIND BY CONSTRUCTION.  Past |eta| = 2.5 the
// alpha sum and the median both run over all candidates, so the weights come
// out identical whichever vertex is under test.  That is not an approximation
// in this code, it is what the algorithm does, and it is the same statement
// that v1 measured as a matched-jet fraction: PUPPI cannot separate vertices
// where it has no tracks.  Anything that merges the jets of several vertices
// has to deal with the resulting duplication in HF.
//
// THE FREE PARAMETERS ARE MEASURED, NOT GUESSED.  CMS sets them in
// PuppiProducer_cfi, which is not in the tuple, so they were read back off the
// stored weights (see the note in puppiprobe.C):
//
//  * the charged-to-vertex rule.  Of the tracker charged candidates,
//    vertexRef == 0 together with fromPVvertexRef == UsedInFit, or PVTight
//    with pvAssociationQuality >= UsedInFitLoose, selects a set that is 100%
//    pure in stored weight 1 - but only 90% efficient, and it carries 89.6%
//    of the pT that CMS put weight 1 on, which shows up directly as a 10%
//    deficit in the jets.  CMS also lets a track in on |dz| alone, and adding
//    |z_track - z_vertex| < 0.03 cm brings the carried pT to 0.99 of the
//    reference at 92% purity.  Both halves of the rule are defined relative
//    to the candidate's own vertex, so they carry over to every other vertex
//    unchanged.
//  * the neutral floor.  It is a hard cut on w*pT rising with the number of
//    vertices, so the lower edge of w*pT against N_PV reads it off directly -
//    on NEUTRAL candidates only, because the charged ones past the tracker
//    carry no floor and hide it.  A straight-line fit to the 0.5% quantile
//    gives 0.184 + 0.0131 N_PV below |eta| = 2.5 and about 1.85 + 0.070 N_PV
//    beyond it, the second being CMS's 2.0 + 0.07 N_PV to the precision this
//    can be measured at.
//  * who the floor applies to.  Charged candidates between |eta| = 2.5 and 3
//    survive at weights far below it, so it is a cut on charge == 0 and not
//    on "everything without a vertex".  Those candidates are outside the
//    tracker, get their weight from the chi2 like a neutral, and are then
//    left alone.
#include <vector>
#include <algorithm>
#include <cmath>
#include <type_traits>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace puppi {

  // Version stamp of this implementation.  hopts is numeric and has no slot
  // for it, so the analyzer prints it at start-up and writes it into the
  // histogram file as the TNamed hist/puppiVersion.
  inline constexpr const char *kVersion = "genseed-1.0";

  // F_chi2,1(x): the cdf of chi2 with one degree of freedom.  The v4 header
  // asked ROOT::Math for it; erf(sqrt(x/2)) is the same function.
  inline double Chi2Cdf1(double x) { return x > 0 ? std::erf(std::sqrt(0.5*x)) : 0.; }

  struct Config {
    double cone        = 0.4;    // R0 of the alpha sum
    double dr2min      = 1e-4;   // pairs closer than this are skipped
    double etaTracker  = 2.5;    // where the alpha sum switches to all-candidate
    double etaVtxAssoc = 2.5;    // how far the charged-to-vertex rule reaches
    double dzCut       = 0.03;   // cm, the |z| path into the same rule
    double rmsPtMin    = 0.1;    // candidates below this do not set med/rms
    double rmsScale    = 0.7;
    bool   leftRMS     = true;   // rms from the low side of the median only
    // Forward regions get their median from every candidate in them, which
    // includes the ones from the vertex under test.  CMS shifts the quantile
    // down by the fraction those make up; medQuantileFwd < 0 asks for that,
    // a positive value fixes the quantile instead.
    double medQuantile    = 0.5;
    double medQuantileFwd = 0.2;
    // The floor, applied only to what the algorithm calls neutral, i.e. to
    // candidates with no vertex of their own.  Region 0 is |eta| < 2.5,
    // 1 is 2.5 to 3, 2 is beyond; 1 and 2 share the measured forward values
    // but keep their own median.
    double etaBound[2]    = {2.5, 3.0};
    double neutralPt[3]   = {0.184, 1.85, 1.85};
    double neutralSlope[3]= {0.0131, 0.070, 0.070};
  };

  // One event's candidates, filled once and reused for every vertex.
  struct Event {
    std::vector<float> pt, eta, phi, mass, ztrk;   // ztrk = PV_z + dz
    std::vector<char>  charged;
    // pdgId and the sign of the charge, not needed by PUPPI itself but by
    // the gen-to-PF linker of v4 (genlink.h), which has to know a track from
    // a photon from a neutral hadron and can veto a track link on the sign.
    // q is signed char, not char: plain char is unsigned on Linux aarch64 and
    // ppc64le, where -1 would read back as 255 and the linker would bend
    // every negative particle the wrong way and 255 times as far.
    std::vector<int>   pdg;
    std::vector<signed char> q;  // -1, 0, +1
    std::vector<int>   vref;     // PFCandV2_vertexRef
    std::vector<char>  fpvref;   // PFCandV2_fromPVvertexRef
    std::vector<char>  paq;      // PFCandV2_pvAssocQuality
    int   nPV = 0;
    // filled by Prepare()
    std::vector<float> alphaFwd;                   // vertex independent
    std::vector<float> wFwd;                       // vertex independent
    bool  prepared = false;
    size_t size() const { return pt.size(); }
    void clear() { pt.clear(); eta.clear(); phi.clear(); mass.clear();
                   ztrk.clear(); charged.clear(); pdg.clear(); q.clear();
                   vref.clear(); fpvref.clear();
                   paq.clear(); alphaFwd.clear(); wFwd.clear();
                   prepared = false; }
    // Is candidate i one of vertex vkey's own charged tracks?  This is the
    // whole per-vertex handle; see the note at the top.  Either it was used
    // in that vertex's fit, or it points at it in z.
    bool fromVertex(size_t i, int vkey, double zv, double dzCut) const {
      if (!charged[i]) return false;
      if (vref[i] == vkey && (fpvref[i] == 3 || (fpvref[i] == 2 && paq[i] >= 6)))
        return true;
      return fabs(ztrk[i] - zv) < dzCut;
    }
  };
  static_assert(std::is_signed<decltype(Event::q)::value_type>::value,
                "puppi::Event::q holds -1 and must be signed");

  inline double dphi(double a, double b) {
    double d = a - b;
    if (d >  M_PI) d -= 2*M_PI;
    else if (d < -M_PI) d += 2*M_PI;
    return d;
  }

  inline int Region(double aeta, const Config &c) {
    if (aeta < c.etaBound[0]) return 0;
    if (aeta < c.etaBound[1]) return 1;
    return 2;
  }

  // Quantile of an unsorted copy; v is sorted in place.
  inline double Quantile(std::vector<float> &v, double q) {
    if (v.empty()) return 0;
    std::sort(v.begin(), v.end());
    const int i = std::min((int)v.size()-1,
                           std::max(0, int(std::floor(q*v.size()))));
    return v[i];
  }

  // RMS about med, optionally from the low side only, as CMS does when
  // applyLowPUCorr is on: the high side is where the signal sits and would
  // inflate the width that the signal is then tested against.
  inline double Rms(const std::vector<float> &v, double med, bool leftOnly) {
    double s = 0; int n = 0;
    for (size_t i = 0; i < v.size(); ++i) {
      if (leftOnly && v[i] > med) continue;
      s += (v[i]-med)*(v[i]-med); ++n;
    }
    return n > 0 ? sqrt(s/n) : 1.;
  }

  // A coarse eta-phi grid, so the all-candidate alpha sum is local rather
  // than quadratic.  2100 candidates and a 0.4 cone would otherwise be
  // 4.4 million pairs an event before the per-vertex loop even starts.
  struct Grid {
    double cell; int neta, nphi; double etamin;
    std::vector< std::vector<int> > cells;
    void build(const std::vector<float> &eta, const std::vector<float> &phi,
               double c, double emin = -6., double emax = 6.) {
      cell = c; etamin = emin;
      neta = std::max(1, int((emax-emin)/cell));
      nphi = std::max(1, int(2*M_PI/cell));
      cells.assign(neta*nphi, std::vector<int>());
      for (size_t i = 0; i < eta.size(); ++i) {
        const int ie = idxEta(eta[i]), ip = idxPhi(phi[i]);
        cells[ie*nphi+ip].push_back(i);
      }
    }
    int idxEta(double e) const {
      return std::min(neta-1, std::max(0, int((e-etamin)/cell))); }
    int idxPhi(double p) const {
      int i = int((p+M_PI)/cell); return ((i%nphi)+nphi)%nphi; }
    // neighbours of (eta,phi) in the 3x3 block of cells
    template <class F> void forNear(double e, double p, F fn) const {
      const int ie = idxEta(e), ip = idxPhi(p);
      for (int de = -1; de <= 1; ++de) {
        const int je = ie+de; if (je < 0 || je >= neta) continue;
        for (int dp = -1; dp <= 1; ++dp) {
          const int jp = ((ip+dp)%nphi+nphi)%nphi;
          const std::vector<int> &c = cells[je*nphi+jp];
          for (size_t k = 0; k < c.size(); ++k) fn(c[k]);
        }
      }
    }
  };

  // Everything that does not depend on which vertex is under test: the alpha
  // and the weight of every candidate outside the tracker.
  inline void Prepare(Event &e, const Config &c)
  {
    const size_t n = e.size();
    e.alphaFwd.assign(n, 0.f);
    e.wFwd.assign(n, 0.f);
    Grid g; g.build(e.eta, e.phi, c.cone);
    const double cone2 = c.cone*c.cone;

    // alpha over every candidate, needed for |eta| > etaTracker
    for (size_t i = 0; i < n; ++i) {
      if (fabs(e.eta[i]) <= c.etaTracker) continue;
      double var = 0;
      const double ei = e.eta[i], pi = e.phi[i];
      g.forNear(ei, pi, [&](int j){
        if ((size_t)j == i) return;
        const double de = ei - e.eta[j], df = dphi(pi, e.phi[j]);
        const double dr2 = de*de + df*df;
        if (dr2 < c.dr2min || dr2 >= cone2) return;
        var += e.pt[j]*e.pt[j]/dr2;
      });
      e.alphaFwd[i] = var > 0 ? log(var) : -99.f;
    }

    // median and rms per forward region, from every candidate in it
    for (int r = 1; r <= 2; ++r) {
      std::vector<float> v;
      for (size_t i = 0; i < n; ++i) {
        const double ae = fabs(e.eta[i]);
        if (ae <= c.etaTracker || Region(ae,c) != r) continue;
        if (e.pt[i] < c.rmsPtMin || e.alphaFwd[i] < -90) continue;
        v.push_back(e.alphaFwd[i]);
      }
      if (v.empty()) continue;
      const double q = c.medQuantileFwd > 0 ? c.medQuantileFwd : c.medQuantile;
      std::vector<float> s(v);
      const double med = Quantile(s, q);
      double rms = Rms(v, med, c.leftRMS) * c.rmsScale;
      // With at most four candidates in the region the 0.2 quantile is the
      // lowest of them, nothing lies below it and the left-side rms is exactly
      // 0; chi2 would be +-inf or NaN.  v4 took that and wrote NaN into the
      // weights.  Take rms = 1 instead, the guard Weights() has in the tracker.
      // It is a guard for low-pileup and special samples: in the 25000 events
      // of NANOAODSIM_1.root the sparsest region still has six candidates.
      if (!(rms > 0)) rms = 1;
      for (size_t i = 0; i < n; ++i) {
        const double ae = fabs(e.eta[i]);
        if (ae <= c.etaTracker || Region(ae,c) != r) continue;
        if (e.alphaFwd[i] < -90) { e.wFwd[i] = 0; continue; }
        const double d = e.alphaFwd[i] - med;
        const double chi2 = d*fabs(d)/(rms*rms);
        e.wFwd[i] = Chi2Cdf1(chi2);
      }
    }
    e.prepared = true;
  }

  // The weights with vertex vkey playing the role of the primary one.
  // Prepare() must have run first.
  inline void Weights(const Event &e, int vkey, double zv, const Config &c,
                      std::vector<float> &w)
  {
    const size_t n = e.size();
    w.assign(n, 0.f);
    const double cone2 = c.cone*c.cone;

    // 0 = neutral or outside the tracker, 1 = charged from this vertex,
    // 2 = charged from somewhere else
    std::vector<char> id(n, 0);
    std::vector<int>  fromV;
    for (size_t i = 0; i < n; ++i) {
      if (!e.charged[i] || fabs(e.eta[i]) > c.etaVtxAssoc) continue;
      const bool own = e.fromVertex(i, vkey, zv, c.dzCut);
      id[i] = own ? 1 : 2;
      // Only the ones inside the tracker feed the alpha sum: past 2.5 the
      // algorithm is track-blind whatever the candidate's own vertex is.
      if (own && fabs(e.eta[i]) <= c.etaTracker) fromV.push_back(i);
    }

    // alpha inside the tracker, summed over this vertex's charged candidates
    std::vector<float> alpha(n, -99.f);
    for (size_t i = 0; i < n; ++i) {
      if (fabs(e.eta[i]) > c.etaTracker) continue;
      double var = 0;
      const double ei = e.eta[i], pi = e.phi[i];
      for (size_t k = 0; k < fromV.size(); ++k) {
        const int j = fromV[k];
        if ((size_t)j == i) continue;
        const double de = ei - e.eta[j], df = dphi(pi, e.phi[j]);
        const double dr2 = de*de + df*df;
        if (dr2 < c.dr2min || dr2 >= cone2) continue;
        var += e.pt[j]*e.pt[j]/dr2;
      }
      alpha[i] = var > 0 ? log(var) : -99.f;
    }

    // median and rms from the charged candidates of the other vertices
    std::vector<float> v;
    for (size_t i = 0; i < n; ++i)
      if (id[i] == 2 && e.pt[i] >= c.rmsPtMin && alpha[i] > -90)
        v.push_back(alpha[i]);
    double med = 0, rms = 1;
    if (!v.empty()) {
      std::vector<float> s(v);
      med = Quantile(s, c.medQuantile);
      rms = Rms(v, med, c.leftRMS) * c.rmsScale;
      if (rms <= 0) rms = 1;
    }

    for (size_t i = 0; i < n; ++i) {
      // A candidate with a vertex of its own is decided by that vertex, and
      // the neutral floor does not apply to it.
      if (id[i] == 1) { w[i] = 1; continue; }
      if (id[i] == 2) { w[i] = 0; continue; }
      const double ae = fabs(e.eta[i]);
      double x;
      if (ae > c.etaTracker) x = e.wFwd[i];
      else if (alpha[i] < -90 || v.empty()) x = 0;
      else {
        const double d = alpha[i] - med;
        const double chi2 = d*fabs(d)/(rms*rms);
        x = Chi2Cdf1(chi2);
      }
      // The floor is a cut on genuinely neutral candidates.  A charged one
      // past the tracker has no vertex to be judged by, so it takes the chi2
      // weight and keeps it.
      if (!e.charged[i]) {
        const int r = Region(ae, c);
        if (x*e.pt[i] < c.neutralPt[r] + c.neutralSlope[r]*e.nPV) x = 0;
      }
      w[i] = x;
    }
  }

  // z of every reconstructed vertex, from the median of PV_z + dz over its
  // tracks.  Only four are stored in NanoAOD and there are ~34 of them.
  // Returns z indexed by vertexRef key; keys with fewer than minTrk tracks
  // get 1e9 and should be skipped.
  inline void VertexZ(const std::vector<float> &ztrk,
                      const std::vector<int> &vref,
                      const std::vector<char> &charged,
                      const std::vector<float> &eta,
                      double etaTracker, int minTrk,
                      std::vector<double> &zv, std::vector<int> &ntrk)
  {
    int kmax = -1;
    for (size_t i = 0; i < vref.size(); ++i)
      if (charged[i] && fabs(eta[i]) <= etaTracker && vref[i] > kmax) kmax = vref[i];
    zv.assign(kmax+1, 1e9);
    ntrk.assign(kmax+1, 0);
    if (kmax < 0) return;
    std::vector< std::vector<float> > z(kmax+1);
    for (size_t i = 0; i < vref.size(); ++i) {
      if (!charged[i] || fabs(eta[i]) > etaTracker) continue;
      z[vref[i]].push_back(ztrk[i]);
    }
    for (int k = 0; k <= kmax; ++k) {
      ntrk[k] = (int)z[k].size();
      if (ntrk[k] < minTrk) continue;
      std::sort(z[k].begin(), z[k].end());
      zv[k] = z[k][z[k].size()/2];
    }
  }

} // namespace puppi
#endif
