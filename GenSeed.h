// genseed - Copyright 2026 Mikko Voutilainen - Apache License 2.0
#ifndef genseed_GenSeed_h
#define genseed_GenSeed_h
// GenSeed.h - generated jets seeded by the reconstruction, and reconstructed
// jets seeded by the generator, through the particle-level link between the
// generated particles and the particle-flow candidates.
//
// THE PROBLEM.  A response matrix needs a pairing between generated and
// reconstructed jets.  In a bunch crossing with 45 interactions, PUPPI run for
// one vertex hypothesis reconstructs the jets of that interaction with part of
// the pileup left in and part of the interaction's own soft neutrals floored
// out.  Pairing a generated jet with the nearest reconstructed one inside
// dR < 0.2 truncates the left tail of R = pT_raw/pT_pure (a badly reconstructed
// jet moves out of the cone) and leaves the right tail in (the pileup the seed
// knows nothing about); ghost association makes the pairing complete but keeps
// the right tail, because the seed still does not know what PUPPI kept.  Both
// are artefacts of the pairing, not of the detector, and both were being
// unfolded as if they were response.
//
// THE IDEA.  Pair where the ambiguity is smallest: between generated particles
// and PF candidates.  A track is one particle seen twice; a calorimeter cluster
// is a small, geometrically predictable set of particles.  genlink.h gives
// every generated particle g its PF partner(s) i with fractions f_gi, once per
// event and independent of any vertex hypothesis, and from that two objects
// follow with no jet-level pairing at all:
//
//   GenJetSeed   the generated energy inside a reconstructed jet, weighted the
//                way PUPPI weighted the jet's candidates;
//   RecoJetSeed  the reconstructed energy that a generated jet's own particles
//                produced, weighted and clustered the way the reconstruction
//                would have done it.
//
// THE GenJetSeed (pT^link).  For vertex hypothesis v the PF candidates with
// PUPPI weight w_i(v) > 0 are clustered at pT_i w_i (anti-kT, R, from
// ptMinCluster) together with ghosts (pT x 1e-9) for the UNLINKED particles of
// the interaction v owns.  The seed of jet k is
//     S_k = sum_{i in k} w_i(v) G_i  +  sum_{unlinked own ghosts g in k} p_g ,
// with G_i = sum_g f_gi p_g the routed generated four-vector of candidate i,
// at the GENERATED pT with the generated mass.  A linked particle is credited
// to the jet its PF partner ended up in ("routed", not clustered), with the
// partner's PUPPI weight; an unlinked own particle lands in whichever jet's
// catchment area its direction falls in, at weight 1; an unlinked far particle
// is nothing.  pT^linkz is the routed sum alone, without the ghosts.  Because
// pT_i never enters, R = pT_raw/pT^link is a detector response - calorimeter
// non-linearity, thresholds, tracking inefficiency, the neutral floor - and a
// correction derived from it transfers to data.  Had the seed used pT_i, R
// would be 1 by construction for everything PUPPI kept.  The seed is not the
// truth either: everything PUPPI decides moves out of R and into the bridge
// seed -> pure, which is measured as a matrix (bridge/hbridge).
//
// THE RecoJetSeed (pT^rs).  For every pure generated jet J of interaction j,
// at the usable vertex v nearest to j if that is within dzRecoSeed: every PF
// candidate i linked to a particle of J enters with the share
//     s_iJ = sum_{g in J} f_gi pT_g / S_i ,   S_i = sum_g f_gi pT_g ,
// the fraction of i's attributed (scalar) generated pT that came from J, in
// [0, 1] and summing to at most 1 over the pure jets, so that a candidate
// shared with another jet or with pileup is never credited twice; weighted by
// s_iJ w_i(v) with the PUPPI weight at v, the set is clustered with anti-kT R
// (no threshold) and the hardest jet is the RecoJetSeed of J.  With it come
// nPF (candidates in that jet), fLead (its pT over the scalar sum of the
// weighted candidates: how much of J's reconstructed energy one jet holds),
// fLinked (the linked fraction of J's generated scalar pT), fPuppi
// (sum s w pT / sum s pT: what PUPPI kept of J's candidates) and pT^rs0, the
// hardest jet of the same candidates clustered at weight s_iJ alone: J as
// the detector saw it before PUPPI decided anything.  pT^rs/pT^pure is what
// the detector made of THIS generated jet - no dR, no share threshold, no
// dedup, no other interaction's particles - but it is not efficiency-free:
// where PUPPI keeps none of J's candidates (most forward jets of a few GeV)
// pT^rs = 0 and the jet is not in hresp_rs.  The response is filled for
// pT^rs > 0, the fraction that has one is eff/hgmat_rs, and the two are to be
// read together; hresp_rs0 is filled wherever anything of J was
// reconstructed.  The reco jet then factorises as
//     reco/pure = (rs/pure) x (reco/rs),
// reco/rs (hresp_reco_rs: pT^raw of the kept copy of J's mutual partner over
// pT^rs) being what pileup and the clustering of the whole event add.  No
// correction is applied to pT^rs: the reco-jet table maps pT^raw to the link
// scale, which a RecoJetSeed is not on.
//
// OWNERSHIP AND ROUTING.  Only four vertex z are stored, so puppi::VertexZ
// recovers every reconstructed vertex from the median PV_z + dz of its tracks
// (25-49 um); a vertex with fewer than minTrkPerVertex tracks is unusable.  A
// vertex OWNS the generated interaction nearest to it if that is within
// dzVertexGen; ownership is unique, so no seed is ever divided among near
// interactions.  An interaction is looked at on the generated side once, at
// the usable vertex nearest to it (nearVtx), if that is within dzRecoSeed;
// its unlinked particles are ghosts there only if that vertex owns it.  The
// seed of a jet is decomposed into own / far (by interaction) / unlinked
// parts (fown + ffar = 1, funl <= fown), and the reconstructed pT into
// attributed / orphan parts (forph: PF candidates with no generated partner -
// the calorimeter excess of a tracked hadron, out-of-time pileup, noise).
//
// TOPOLOGY.  Reco side, per jet k at vertex v: the pure jets giving more than
// shareFrac of pT^link are its contributors (n_in); p* is the largest; n_out
// is the number of jets at v receiving more than shareFrac of p*.  Classes:
// none (n_in = 0), 1to1 (1, 1), split (1, n_out != 1), merge (>= 2, 1),
// tangle (>= 2, >= 2); they partition the reconstructed jets.  Gen side, per
// pure jet J at its vertex: n_r jets receive more than shareFrac of pT^pure;
// kmain is the one receiving most; lost (n_r = 0), 1to1 (1 and kmain has at
// most one contributor), split (>= 2, kmain has at most one), merge (1, kmain
// has >= 2), tangle (>= 2, >= 2); they partition the pure jets with a vertex
// (a jet whose interaction has none keeps gj_topo = lost in the tuple but is
// not in topo/).  The count is kmain's n_in whether or not J is among them:
// a kmain that is mostly far pileup (n_in = 0, reco side none) still makes J
// 1to1 or split.
// MUTUAL DOMINANCE: k and J are a pair when J is k's largest contributor and
// gave it more than shareFrac of pT^link, k is J's main jet and k took more
// than shareFrac of J - one relation, bpair on the reco side and dom on the
// gen side, so eff/hgmat_dom and the pairing behind fake/hrunpaired and
// bridge/ differ only by the dedup (which keeps one copy of k).  The bridge
// seed -> pure is a matrix over the mutual pairs, with the unpaired seeds
// (hbfake) and the unpaired pure jets (hbmiss; also those whose interaction
// has no vertex) as additive vectors, so that hseed = hbfake + proj_x and
// hgen = hbmiss + proj_y hold exactly, per |y| bin and half.  The classic dR < 0.2 partner is kept
// alongside (gj_dr02reco, resp/hresp_dr).
//
// DEDUP.  Past |eta| = 2.5 PUPPI is vertex-blind, so the same forward jet comes
// out of every vertex hypothesis; a neutral-only barrel jet does too.  The
// jets of all vertices are collected, copies (dR < dupDR, |dpT| < dupPtFrac)
// are grouped, the copy with the largest own fraction of pT^link is kept
// (pT order breaks ties) and every response histogram is filled ONCE per kept
// jet.  rj_ncopies says how many hypotheses produced the jet.
//
// THRESHOLDS.  Reco jets are clustered from ptMinCluster = 1 GeV raw and pure
// gen jets from ptMinGen = 1 GeV, far below the 5 GeV the measurement is
// reported from, so that the 2-5 GeV bins of the extended axis are populated
// and threshold migrations are inside the matrix rather than in miss and
// fake.  The tuple stores reco jets above ptStoreReco and gen jets above
// ptStoreGen (3 GeV each); an index to a partner below the threshold is -1.
// Generated particles beyond etaMaxPart = 5.5 and neutrinos are dropped.
//
// KNOBS (all through SetOpt, all written to hist/hopts as value x jobs):
//   R              0.4    jet radius of all four collections
//   ptMinCluster   1      reco jets clustered from this PUPPI-weighted raw pT [GeV]
//   ptMinGen       1      pure gen jets clustered from this [GeV]
//   ptStoreReco    3      tuple: reco jets stored above this raw pT [GeV]
//   ptStoreGen     3      tuple: pure gen jets stored above this [GeV]
//   etaMaxPart     5.5    generated particles kept
//   shareFrac      1/3    topology and mutual-dominance share
//   dzVertexGen    0.05   cm, a vertex owns the nearest interaction within this
//   dzRecoSeed     0.2    cm, an interaction is processed at nearVtx within this
//   minTrkPerVertex 3     vertices with fewer tracks are unusable
//   dupDR          0.15   dedup across hypotheses: dR
//   dupPtFrac      0.10   dedup across hypotheses: relative pT
//   kernelScale    1      linker: scales every pass-2 sigma and cone
//   useDz          1      linker: pass-1 dz cut on
//   useCharge      1      linker: pass-1 same charge sign required
//   parallaxRef    1      linker: calo directions seen from PV_z (0: from z = 0)
//   bendSign       +1     linker: sign of the charged bend (link/hbend checks it)
//   jecFile        ""     correction table from pass 1 (SetJEC), applied to reco jets
//                         only; jecLoaded in hopts, the name in hist/jecFile
//                         ("none" in pass 1); a table given that does not load
//                         stops the job before anything is written
//   maxVertices    -1     stop after this many vertices (timing tests)
//   firstEntry, lastEntry  the job's [first, last), even-aligned by the driver;
//                         last <= first is an empty job: a valid file, no events
//   progressEvery  200    print a progress line every this many events
//   useFastJet     compile time (-DGENSEED_USE_FASTJET), recorded in hopts
//
// OUTPUT (names are the contract with drawGenSeed.C; <y> = y00..y09, <e> =
// e00..e09, <h> = a|b by entry parity, pT axis kPtBinsExt, R axis 0-5/500):
//   genseed                    TTree, one entry per event, std::vector branches
//                              (gj_rs_* and gj_rs0_pt: -1 / 0 / > 0, see RSeed)
//   hist/hcount hopts hnghost hlinkms hncopies
//   hist/puppiVersion jecFile  TNamed: puppi::kVersion, the table or "none" (hadd: a cycle per job)
//   spec/hgen_<y>_<h> (pure) hreco_<y>_<h> (corr, kept) hseed_<y>_<h> (pT^link)
//   spec/hrseed_<y>_<h> (pT^rs > 0, in the pure jet's |y|, no cut on pT^pure)
//   eff/hgall_<e> hgmat_dr_<e> hgmat_dom_<e> hgmat_rs_<e> hgmat_link_<e>      vs pT^pure
//   fake/hrall_<e> hrnone_<e> hrunpaired_<e>                                 vs pT^corr
//   resp/hresp_link_<e> hresp_linkz_<e> (vs pT^link, pT^linkz) hresp_dr_<e> (vs pT^pure)
//   resp/hresp_rs_<e> (pT^rs/pT^pure, pT^rs > 0) hresp_rs0_<e> (pT^rs0/pT^pure, pT^rs0 > 0)
//   resp/hresp_reco_rs_<e> (pT^raw of J's kept mutual partner / pT^rs, pT^rs > 0), vs pT^pure
//   resp/hrtopo_link_<topo>_<e> hresppu_link_pu<0,1,2>_<e> hfar_<e> hunl_<e> horph_<e>
//   topo/hall_<y> h<topo>_<y> (gen: 1to1 split merge tangle lost, vs pT^pure, with a vertex)
//   topo/hr<topo>_<y> (reco: 1to1 split merge tangle none, vs pT^corr)
//   bridge/hbridge_<y>_<h> (pT^link, pT^pure) hbfake_<y>_<h> hbmiss_<y>_<h>
//   map/hjet_reco hjet_gen hjet_rs presp_link presp_rs porph pflink pfar
//   link/hlinkfrac_<class>_<reg> hparallax hbend horphan_<pfclass>_<reg> hlinkdz_<reg> hpfratio_<pfclass>
//
// CLUSTERING BACKEND.  tiledjet.h (clean-room anti-kT, tiled strategy) by
// default; FastJet under -DGENSEED_USE_FASTJET, which runGenSeed.C sets when
// FastJet is found and tiledjet::kAtLeastAsFastAsFastJet is false.  Both give
// E-scheme anti-kT with the same constituent sets; genseed::ClusterBackend()
// says which one a file was made with.
#include "puppi.h"
#include "genlink.h"
#include "jec.h"
#include "nanoreader.h"

#ifdef GENSEED_USE_FASTJET
#include "fastjet/ClusterSequence.hh"
#else
#include "tiledjet.h"
#endif

#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <map>
#include <string>
#include <vector>
#include <cmath>

namespace genseed {

  // ---- axes -------------------------------------------------------------
  // The CMS inclusive-jet pT edges extended down to 2 GeV: the reco jets are
  // clustered from 1 GeV and the seeds have no threshold, so nothing forces
  // a histogram to start where the measurement is reported.
  const double kPtBinsExt[] = {
      2,   3,   4,   5,   6,   8,  10,  12,  15,  18,  21,  24,  28,
     32,  37,  43,  49,  56,  64,  74,  84,  97, 114, 133, 153, 174,
    196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592,
    638, 686, 737, 790 };
  const int kNPtExt = sizeof(kPtBinsExt)/sizeof(kPtBinsExt[0]) - 1;   // 42
  inline int PtBinExt(double pt) {
    if (pt < kPtBinsExt[0] || pt >= kPtBinsExt[kNPtExt]) return -1;
    int lo = 0, hi = kNPtExt;
    while (hi - lo > 1) { const int m = (lo+hi)/2; if (pt < kPtBinsExt[m]) hi = m; else lo = m; }
    return lo;
  }
  const int    kNY = 10;                 // |y| bins 0.5 wide to 5.0
  const double kYWidth = 0.5;
  inline int YBin(double y) { const int i = int(std::fabs(y)/kYWidth); return (i >= 0 && i < kNY) ? i : -1; }
  inline std::string YTag(int i)   { char b[8];  snprintf(b, sizeof(b), "y%02d", i); return b; }
  inline std::string YLabel(int i) { char b[48]; snprintf(b, sizeof(b), "%.1f < |y| < %.1f", i*kYWidth, (i+1)*kYWidth); return b; }
  const int    kNR = 500; const double kRMax = 5.;   // response axis
  const double kDRMatch = 0.2;           // the classic association cone
  const double kGhost = 1e-9;            // ghost scale: too small to move an axis, far from underflow
  // No cross section to read off a NeutrinoGun: the interactions ARE the
  // events, L = sum(nPU)/sigma_inel with the Run 3 value 80 mb; 1 mb = 1e9 pb.
  const double kSigmaInelMb = 80.0;
  inline double LumiPbInv(double nInter) { return nInter/(kSigmaInelMb*1e9); }

  // ---- the clustering, one call over either backend ----------------------
  struct Jet { double pt, y, eta, phi, m; int ncon; };

  inline const char *ClusterBackend() {
#ifdef GENSEED_USE_FASTJET
    return "FastJet";
#else
    return "tiledjet";
#endif
  }

  // E-scheme anti-kT of the Cartesian four-vectors, jets above ptmin in pT
  // order, cons[k] the input indices of jet k.  Kinematics are computed here
  // from the jet four-vector so both backends give identical conventions
  // (phi in (-pi, pi], y the true rapidity, m from E^2 - p^2 clipped at 0).
  inline std::vector<Jet> Cluster(const std::vector<double> &px, const std::vector<double> &py,
                                  const std::vector<double> &pz, const std::vector<double> &E,
                                  double R, double ptmin, std::vector< std::vector<int> > &cons)
  {
    std::vector<Jet> out;
    cons.clear();
    const int N = (int)px.size();
    if (N == 0) return out;
    auto finish = [&](double jx, double jy, double jz, double jE, std::vector<int> &ix) {
      Jet j;
      const double pt2 = jx*jx + jy*jy, p = sqrt(pt2 + jz*jz);
      j.pt = sqrt(pt2); j.phi = atan2(jy, jx);
      j.eta = (p > std::fabs(jz)) ? 0.5*log((p+jz)/(p-jz)) : (jz >= 0 ? 1e3 : -1e3);
      j.y   = (jE > std::fabs(jz)) ? 0.5*log((jE+jz)/(jE-jz)) : j.eta;
      const double m2 = jE*jE - p*p; j.m = m2 > 0 ? sqrt(m2) : 0.;
      j.ncon = (int)ix.size();
      out.push_back(j); cons.push_back(ix);
    };
#ifdef GENSEED_USE_FASTJET
    std::vector<fastjet::PseudoJet> in; in.reserve(N);
    for (int i = 0; i < N; ++i) { in.push_back(fastjet::PseudoJet(px[i],py[i],pz[i],E[i])); in.back().set_user_index(i); }
    fastjet::JetDefinition jd(fastjet::antikt_algorithm, R);
    fastjet::ClusterSequence cs(in, jd);
    std::vector<fastjet::PseudoJet> js = fastjet::sorted_by_pt(cs.inclusive_jets(ptmin));
    for (size_t k = 0; k < js.size(); ++k) {
      const std::vector<fastjet::PseudoJet> cv = js[k].constituents();
      std::vector<int> ix; ix.reserve(cv.size());
      for (size_t m = 0; m < cv.size(); ++m) ix.push_back(cv[m].user_index());
      finish(js[k].px(), js[k].py(), js[k].pz(), js[k].E(), ix);
    }
#else
    std::vector<tiledjet::PseudoJet> in(N);
    for (int i = 0; i < N; ++i) { in[i].px = px[i]; in[i].py = py[i]; in[i].pz = pz[i]; in[i].E = E[i]; in[i].user_index = i; }
    tiledjet::Config c; c.R = R; c.p = -1; c.ptmin = ptmin;
    std::vector<tiledjet::Jet> js = tiledjet::cluster(in, c);
    for (size_t k = 0; k < js.size(); ++k) {
      std::vector<int> ix(js[k].constituents);
      finish(js[k].p.px, js[k].p.py, js[k].p.pz, js[k].p.E, ix);
    }
#endif
    return out;
  }

} // namespace genseed

class GenSeed {
public:
  // Same enum, two tag sets: gen side (lost = no jet took a third of it),
  // reco side (none = no pure jet gave a third of the seed).
  enum Topo { kOneOne = 0, kSplit, kMerge, kTangle, kLost, kNTopo };
  static const char *TopoTag(int t)  { static const char *s[kNTopo] = {"1to1","split","merge","tangle","lost"}; return (t>=0&&t<kNTopo)?s[t]:"?"; }
  static const char *RTopoTag(int t) { static const char *s[kNTopo] = {"1to1","split","merge","tangle","none"}; return (t>=0&&t<kNTopo)?s[t]:"?"; }
  static const char *TopoLabel(int t) {
    static const char *s[kNTopo] = {"1 #leftrightarrow 1","1 gen #rightarrow 2+ reco (split)",
      "2+ gen #rightarrow 1 reco (merge)","2+ #leftrightarrow 2+ (tangle)","no partner"};
    return (t>=0&&t<kNTopo)?s[t]:"?"; }

  // One reco jet of one vertex hypothesis with its GenJetSeed.
  struct RJet {
    double pt, ptcorr, y, eta, phi, m;
    int    vtx, gen;               // vertex key, its own interaction (-1)
    double link, linkz;            // seed pT, routed + own ghosts / routed only
    double fown, ffar, funl, forph;
    int    topo, pstar;            // reco-side topology, dominant pure jet (-1)
    double pshare;                 // what p* gave this seed [GeV]
    char   bpair;                  // mutual dominance with p*
    int    dompure;                // pure jet mutually dominant with this jet (-1)
    int    ncopies, rep;           // dedup: copies in its group, index of the kept copy in fReco
  };
  // The RecoJetSeed of a pure gen jet, in three states (the tuple's gj_rs_*):
  //   -1 everywhere  no usable vertex within dzRecoSeed (vtx = -1, gj_vtx < 0);
  //   pt = 0         a vertex, but no linked candidate, or none PUPPI kept:
  //                  eta, phi, m, npf, flead are 0 and mean nothing, flinked,
  //                  fpuppi and pt0 are valid;
  //   pt > 0         a RecoJetSeed.  "Has one" is pt > 0 (eff/hgmat_rs), not pt >= 0.
  // pt0 is the hardest jet of the same candidates at weight s alone (no
  // PUPPI): -1 without a vertex, 0 if nothing of J was reconstructed.
  struct RSeed { double pt, eta, phi, m; int npf; double flead, flinked, fpuppi; int vtx; double pt0; };
  // One pure generated jet with everything found for it at its vertex.
  struct PJet {
    double pt, y, eta, phi, m; int ncon, inter;
    int    vtx;                    // the vertex (key) it was looked at, -1 none
    int    dr02, kmain, dom;       // fReco indices: dR < 0.2 partner, main receiver, mutual partner (-1)
    int    topo; double share;     // gen-side topology, pT given to kmain
    RSeed  rs;                     // -1 / 0 / > 0 as above; rs.vtx = vtx
  };

  GenSeed(TTree *tree, const char *outname = "rootfiles/GenSeed.root");
  ~GenSeed();
  void   SetOpt(const std::string &key, double val) { fOpt[key] = val; }
  double Opt(const std::string &key, double def) const {
    std::map<std::string,double>::const_iterator it = fOpt.find(key); return it == fOpt.end() ? def : it->second; }
  void   SetJEC(const char *f) { fJECFile = f; }
  void   Loop();

private:
  void Book();
  void BookTree();
  void ReadGen();                 // GenList, four-vectors, pure gen jets per interaction
  void ReconstructVertices();     // z per vertex, ownership, nearest vertex per interaction
  void LinkDiagnostics();
  void ProcessVertex(int v);      // PUPPI, clustering with own ghosts, seeds, topology, gen side, RecoJetSeeds
  void RecoJetSeed(int p, int v);   // with the weights fW of vertex v
  void FillEvent(int half);       // dedup, histograms, tuple
  void CheckIdentities();
  void Write();
  static double DeltaR2(double e1, double p1, double e2, double p2);

  std::map<std::string,double> fOpt;
  std::string fOutName, fJECFile;
  TTree *fTree; TFile *fOut; TTree *fTup;
  nanoaod::Reader *fReader;
  jec::Correction fJEC;
  puppi::Config   fPupCfg;
  genlink::Config fLinkCfg;
  double fR, fShare, fPtMinCluster, fPtMinGen, fPtStoreReco, fPtStoreGen, fDzRecoSeed;

  // generated side, flattened over interactions
  genlink::GenList    fGen;
  std::vector<double> fGPx, fGPy, fGPz, fGE;
  std::vector<int>    fGPure;               // pure jet index per particle, -1
  std::vector<PJet>   fPure;
  std::vector< std::vector<int> > fPureCons; // gen indices per pure jet
  std::vector<int>    fNInter;              // particles kept per interaction
  std::vector<double> fVz;                  // z per interaction
  int fNPU, fNSlot;
  // reconstructed side
  puppi::Event fPF;
  std::vector<float>  fW;
  std::vector<double> fZv;
  std::vector<int>    fNtrk, fGenOf, fNearVtx, fVtxIdx;   // fVtxIdx: key -> index in the usable list
  std::vector<int>    fUsable;              // keys of the usable vertices
  genlink::Links fLinks;
  std::vector<double> fSpt;                 // per PF: S_i = sum_g f_gi pT_g, the scalar attributed gen pT
  std::vector<RJet>   fReco;
  double fLinkMs; Long64_t fNEvents;

  // tuple branches
  UInt_t bRun, bLumi; ULong64_t bEvent; Int_t bNPU, bNVtx; Float_t bPVz, bLumiW;
  std::vector<float> bVtxZ; std::vector<int> bVtxNtrk, bVtxGen;
  std::vector<float> bGjPt, bGjEta, bGjPhi, bGjMass, bGjShare;
  std::vector<int>   bGjInter, bGjNcon, bGjVtx, bGjDr02, bGjTopo, bGjKmain;
  std::vector<float> bRsPt, bRs0Pt, bRsEta, bRsPhi, bRsMass, bRsFlead, bRsFlinked, bRsFpuppi;
  std::vector<int>   bRsNpf;
  std::vector<float> bRjPt, bRjPtcorr, bRjEta, bRjPhi, bRjMass, bRjLink, bRjLinkz,
                     bRjFown, bRjFfar, bRjFunl, bRjForph, bRjPshare;
  std::vector<int>   bRjVtx, bRjGen, bRjTopo, bRjPstar, bRjBpair, bRjNcopies;

  // histograms
  TH1D *hcount, *hopts, *hnghost, *hlinkms, *hncopies;
  TH1D *hspec_gen[genseed::kNY][2], *hspec_reco[genseed::kNY][2], *hspec_seed[genseed::kNY][2], *hspec_rs[genseed::kNY][2];
  TH1D *hgall[jec::kNEta], *hgmat_dr[jec::kNEta], *hgmat_dom[jec::kNEta], *hgmat_rs[jec::kNEta], *hgmat_link[jec::kNEta];
  TH1D *hrall[jec::kNEta], *hrnone[jec::kNEta], *hrunpaired[jec::kNEta];
  TH2D *hresp_link[jec::kNEta], *hresp_linkz[jec::kNEta], *hresp_dr[jec::kNEta], *hresp_rs[jec::kNEta];
  TH2D *hresp_rs0[jec::kNEta], *hresp_reco_rs[jec::kNEta];
  TH2D *hrtopo[kNTopo][jec::kNEta], *hresppu[3][jec::kNEta];
  TProfile *hfar[jec::kNEta], *hunl[jec::kNEta], *horph[jec::kNEta];
  TH1D *htall[genseed::kNY], *htopo[kNTopo][genseed::kNY], *hrtopo1[kNTopo][genseed::kNY];
  TH2D *hbridge[genseed::kNY][2]; TH1D *hbfake[genseed::kNY][2], *hbmiss[genseed::kNY][2];
  TH2D *hjet_reco, *hjet_gen, *hjet_rs;
  TProfile2D *presp_link, *presp_rs, *porph, *pflink, *pfar;
  TProfile *hlinkfrac[genlink::kNDiagClass][3];
  TH1D *hlinkdz[3]; TH2D *hparallax, *hbend;
  TH1D *hpfratio[genlink::kNPFClass]; TProfile *horphan[genlink::kNPFClass][3];
};

#endif
