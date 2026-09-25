// genseed - Copyright 2026 Mikko Voutilainen - Apache License 2.0
#ifndef genseed_JVA_h
#define genseed_JVA_h
// JVA.h - the forward jets of every vertex, associated by global MET
// minimisation (jvassoc.h), against the ways of doing without: the analyzer
// that measures what the association buys in the MET, in the 2.5 < |eta| <
// 3.0 jet spike, and in the dijet balance and resolution down to a few GeV.
//
// THE QUESTION.  genseed runs PUPPI once per reconstructed vertex, and past the
// tracker PUPPI is vertex-blind: every vertex gets the same forward jets.  The
// PV-only view (CMS) hides this by looking at one vertex; genseed, looking at
// all of them, sees each forward jet forty times, sees it made of several
// interactions (the dominant one gives less than 0.4 of the linked pT for 71%
// of them), and sees the MET of every vertex at 23-25 GeV (sqrt <|MET|^2>)
// against 6 GeV of true imbalance.  The forward energy has to be given to
// somebody.  jvassoc.h gives each forward cluster to the vertex whose
// vertex-resolved MET it balances best over all vertices at once, or to
// nobody; this file asks whether that is better than the obvious
// alternatives, and by how much.
//
// SEVEN WAYS TO GIVE OUT THE FORWARD ENERGY, each a complete per-vertex event
// for every usable vertex: the vertex-resolved candidates (|eta| <= etaFwd)
// with that vertex's PUPPI weights, plus the vertex-blind candidates the
// method gives it, at wFwd:
//   dup    every vertex gets all of it: genseed per-vertex PUPPI as it is, and
//          CMS PUPPI at the PV;
//   lv     only the PV (vertex key 0) gets it: the leading-vertex approach,
//          extended to the other vertices, which get no forward energy;
//   trk    per forward cluster, its candidate with the most evidence
//          (jva::Evidence: the in-cone track pT T_ku in candMode 0, the
//          in-cluster charged pT C_ku in 1, their sum in 2; the lowest index on
//          a tie), else the PV: tracks alone;
//   jva    per cluster, jva::Assign: this work;
//   ojet   per cluster, the oracle: the nearest usable vertex of the
//          cluster's dominant interaction, if that gives more than half of
//          the cluster's linked generated pT and the vertex is within
//          dzOracle of it, else null - the best any whole-cluster rule can do;
//   opart  per CANDIDATE, the same oracle: the upper bound that splitting the
//          clusters would reach;
//   none   nobody gets any of it: every vertex is its vertex-resolved part
//          alone.  The reference any association has to beat - at the
//          default knobs jva came out close to it, and "correct" (right
//          vertex or combination dropped) of none is the fraction of
//          combinations, which a method must exceed to have placed anything.
// Without a usable PV (key 0 below minTrkPerVertex tracks; counted and
// printed), lv and the fallback of trk give to nobody, and wFwd is read at the
// first usable vertex.  The dominant interaction comes from the genseed links: the linked
// generated pT wFwd_i f_gi pT_g summed by interaction slot; share = top /
// total; unlinked -> null.  Each per-vertex event is clustered (anti-kT R from
// ptMinCluster, the genseed jets), corrected with the genseed table (jec.h,
// text/jec_v1.txt: pT^corr on the link scale), and its MET is -(sum of the
// weighted pT of everything in it), type-1 corrected with the jets above
// t1Min.  Truth: every jet has its GenJetSeed without ghosts, pT^linkz =
// |sum over its PF constituents of w_i G_i| with the weight the constituent
// has in THAT event, and every owned vertex has the true MET of its owner
// interaction, MET_true(i) = -(sum of the generated pT of i's particles with
// |eta| < 5.0), neutrinos already gone - what a perfect |eta| < 5 detector
// would see.
//
// THE RECOIL jva balances (recoilMode, jvassoc.h THE RECOIL): 0 the particle
// MET M_v of the vertex-resolved part, 1 the MHT of the vertex-resolved jets
// above mhtMin.  Those jets are the jets of the none event, clustered once
// per vertex before the assignment (VertexRecoil) and reused for every
// method that gives the vertex no forward candidate - the same input in the
// same order, so the same jets; it halves the clustering, since most
// vertices get nothing under most methods.
//
// WHAT IS MEASURED (the contract with drawJVA.C; <m> = dup lv trk jva ojet
// opart none, <c> = all pv pu, pT axis jva::kPtBinsExt = genseed::kPtBinsExt):
//   hist/hcount        events_a events_b npu_a npu_b entries jobs vertices owned
//   hist/hopts         every knob, value x jobs; puppiVersion jecFile backend (TNamed)
//   met/<m>/hmetxy_<c> both MET components (raw), 400 bins -100..100
//          hmet_<c>    |MET| 200 bins 0..200
//          hmetres_<c> components of MET - MET_true(owner), owned vertices
//          hmetxy_npu  N_PU x component, all vertices
//          hmetxy_s    S_v x component;  hmetres_s  S_v x residual, owned
//   jets/<m>/ (owned vertices) and jets/gen/ (pure GenJets of the owner of every
//          owned vertex, once per owned vertex):
//          hjeteta_pt5 hjeteta_pt10 hjeteta_pt20   |eta| of jets above pT^raw (gen: pT)
//          hjetpt_r0 r1 r2    pT^corr (gen: pT) at |eta| < 2.5, 2.5-3.0, 3.0-5.0
//          pown_eta           owner fraction of the jet's linked generated pT, pT^raw > 5
//          hdom_r0 r1 r2      dominant-interaction share of it, pT^raw > 5
//                             (a jet with nothing linked counts 0 in both; a share of
//                             exactly 1 is filled at 1 - 1e-9, inside the last bin)
//   fwd/hclpt_trk hclpt_notrk hshare hcand;  fwd/<m>/hcat_trk hcat_notrk (m = lv trk
//          jva ojet none): cluster pT x outcome right wrong nulled_single nulled_combo
//          kept_combo.  "trk" / "notrk" and hcand are the candidates under candMode
//          (in candMode 0 the track candidates of v1)
//   dijet/<m>/<c>/ (c = all pv): hadb_e<ie>_a<ia> hampf_ hrtrue_ (tag-probe, probe |eta|
//          in the 18 CMS L2Res bins e00..e17, a0 alpha < 0.3, a1 none), hdbpar_e<ie>
//          hdbperp_ hmpfpar_ hmpfperp_ (same-bin dijets, bisector axis), hrlink_e<ie>
//          (pT^corr/pT^linkz of both jets), hsel (selection counts: vertices, >=2 jets,
//          dphi, tag = a barrel jet among the two, alpha<0.3 with a tag, samebin =
//          both in one |eta| bin with alpha < 0.3, tag or not)
//   jva                TTree, one entry per crossing (storeTuple), floats and ints;
//                      indices are into the usable-vertex list of the entry, -1 null:
//     per crossing     run lumi event nPU nVtxUsable jva_ncomp jva_nexact jva_objective
//     per vertex       vtx_key vtx_z vtx_owner (slot, -1)
//                      vtx_mx vtx_my vtx_s: the particle MET M_v and S_v (whatever recoilMode)
//                      vtx_ht<T> vtx_htsq<T> vtx_mht<T>_x vtx_mht<T>_y, T = 3 5 7 10 15:
//                        the vertex-resolved jets above T in pT^raw, sum pT, sum pT^2 and
//                        MHT = -sum p (the sign of M_v); recoilMode 1 balances
//                        vtx_mht<mhtMin> with sigma^2 from vtx_s, vtx_ht, vtx_htsq
//                      vtx_genmet_x/_y   MET_true(owner) = -(generated pT of the owner, |eta| < 5)
//                      vtx_genfwd_x/_y   +(generated pT of the owner, etaFwd < |eta| < 5): what
//                                        that interaction's forward clusters carry
//                      vtx_gencen_x/_y   +(generated pT of the owner, |eta| <= etaFwd), so that
//                                        genmet = -(gencen + genfwd)
//                      vtx_genmht<T>_x/_y  +(pT of the owner's pure GenJets, |eta| < etaFwd, pT > T)
//                      (the gen ones 0 without an owner; MHT and gencen/genmht have
//                      opposite signs by definition: a perfect MHT is -genmht)
//     per cluster      cl_pt cl_eta cl_phi cl_px cl_py; cl_nopt: jva's candidates under candMode
//                      (0 = none: any vertex at tauAll);
//                      cl_ncand cl_candoff cand_vtx cand_trkpt: the in-cone track evidence,
//                        u with T_ku >= trkMinPt (v1's track candidates, whatever candMode);
//                      cl_nch cl_choff ch_vtx ch_pt: the in-cluster charged evidence, u with
//                        C_ku > 0 (raw pT), both flattened, cluster k's at [off_k, off_k + n_k);
//                      cl_dom cl_share cl_true: dominant slot, its share, true vertex (-1);
//                      cl_top3_slot cl_top3_share: the three leading slots of the linked pT
//                        and their shares, 3 per cluster (-1 and 0 where fewer);
//                      cl_host_lv cl_host_trk cl_host_jva cl_host_ojet cl_host_none
//
// THE DIJETS.  Every vertex of every method is a dijet event if its two
// leading jets (pT^corr > ptMinJet, |eta| < 5.191) are back to back (dphi >
// dphiMin).  The tag is the jet with |eta| < 1.305, the probe the other; with
// both in the barrel, both assignments are filled at weight 0.5.  A_DB =
// (pT_probe - pT_tag)/(pT_probe + pT_tag) and A_MPF = MET_T1.u_tag/(2 pT_avg)
// give R = (1 + <A>)/(1 - <A>); r_true is the same ratio on the truth scale
// (pT^corr/pT^linkz of the probe over that of the tag).  The resolution comes
// from the dijets with both jets in one |eta| bin: the bisector n = (u1 -
// u2)/|u1 - u2| makes the same angle with both jets, so (p1 + p2).n carries
// their pT imbalance and (p1 + p2).n_perp only the radiation and the angular
// smearing; JER = sqrt(RMS_par^2 - RMS_perp^2)/sqrt(2), from the jets (DB) or
// from MET_T1 (MPFX).  The labels 1 and 2 are arbitrary for two jets in the
// same bin, and the leading-first order would fold the parallel projection
// onto one side, so each same-bin dijet is filled at +x and -x with weight
// 0.5: the distributions are symmetric by construction and their RMS about
// zero is the width.  hrlink holds the truth, pT^corr/pT^linkz of both jets.
//
// CONVENTIONS, the genseed ones exactly: usable vertices from puppi::VertexZ
// with minTrkPerVertex tracks; a vertex owns the nearest interaction within
// dzVertexGen (0.05 cm); nearVtx(i) is the usable vertex nearest to
// interaction i; generated particles to |eta| < 5.5 without neutrinos; halves
// a/b by global entry parity; the clustering backend of GenSeed.h (a copy of
// its wrapper below, so that this file does not drag class GenSeed into the
// dictionary).  The forward clusters are always tiledjet (jvassoc.h is ROOT-
// and FastJet-free); the two backends agree jet by jet on weighted input.
//
// IDENTITIES, printed at the end, worst violation must be 0: under dup every
// vertex gets every forward candidate; under every other method each
// vertex-blind candidate with wFwd > 0 goes to at most one vertex; under none
// no vertex gets any, and every cluster is null; every such candidate is in
// exactly one cluster; the weight of every vertex-blind candidate is the same
// at every vertex; the dup and lv events at the PV are the same event (their
// MET is compared crossing by crossing); ojet and opart agree with the truth
// recomputed from the gen side of the links, and ojet is never wrong, never
// nulls a single-interaction cluster, never keeps a combination; sum over
// outcomes of hcat = hclpt, bin by bin; the in-cluster charged evidence of a
// cluster adds up to the raw pT of its charged constituents on usable
// vertices, recounted here, and its candidates are the candMode rule applied
// to the evidence; MET_true = -(gencen + genfwd) for every interaction.
//
// KNOBS (SetOpt; all written to hist/hopts as value x jobs):
//   R 0.4, ptMinCluster 1, etaFwd 2.5, trkMinPt 0.5, sigma0 1, sigmaK 1,
//   tauTrk 4, tauAll 12, maxCombos 200000, maxSweeps 50 (jvassoc.h),
//   ptMinJet 3, t1Min 3, dphiMin 2.7, dzVertexGen 0.05, dzOracle 0.2,
//   minTrkPerVertex 3, storeTuple 1, jecLoaded (1 if a table was loaded),
//   firstEntry, lastEntry, progressEvery 200, useFastJet (compile time),
//   and appended in v2 (the defaults are v1): candMode 0, chMinPt 0.5,
//   recoilMode 0, mhtMin 10, sigmaJ 0.35 (jvassoc.h).
//   The generated-particle acceptance (|eta| < 5.5) and the MET_true one
//   (|eta| < 5.0) are constants, as are the dijet cuts |eta| < 5.191, tag
//   |eta| < 1.305 and alpha < 0.3.
//
// PRINTED AT THE END: the MET rms per component for every method and vertex
// class (and of MET - MET_true for the owned vertices, and of MET_true
// itself), the forward clusters and components per crossing, the fraction of
// components solved exactly, the recoil and its mean sigma_v, how often the
// true vertex of a single-interaction cluster above 5 GeV is among the
// candidates and leads the evidence (for either kind of evidence, whatever
// candMode), the outcome fractions above 5 GeV, the dijet selection, the
// identities, ms per crossing.
#include "jvassoc.h"
#include "puppi.h"
#include "genlink.h"
#include "jec.h"
#include "nanoreader.h"

#ifdef GENSEED_USE_FASTJET
#include "fastjet/ClusterSequence.hh"
#endif

#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TProfile.h>
#include <map>
#include <string>
#include <vector>
#include <cmath>

namespace jva {

  // genseed::kPtBinsExt, copied: the CMS inclusive-jet edges extended to 2 GeV.
  const double kPtBinsExt[] = {
      2,   3,   4,   5,   6,   8,  10,  12,  15,  18,  21,  24,  28,
     32,  37,  43,  49,  56,  64,  74,  84,  97, 114, 133, 153, 174,
    196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592,
    638, 686, 737, 790 };
  const int kNPtExt = sizeof(kPtBinsExt)/sizeof(kPtBinsExt[0]) - 1;   // 42

  // The CMS L2Res |eta| bins of the probe, e00..e17.
  const double kEtaDj[] = { 0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.930, 2.172,
                            2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191 };
  const int kNEtaDj = sizeof(kEtaDj)/sizeof(kEtaDj[0]) - 1;          // 18
  inline int EtaBinDj(double aeta) {
    if (!(aeta >= kEtaDj[0]) || aeta >= kEtaDj[kNEtaDj]) return -1;
    int i = 0; while (i < kNEtaDj - 1 && aeta >= kEtaDj[i+1]) ++i;
    return i;
  }

  // ---- the clustering, one call over either backend (GenSeed.h's wrapper) --
  struct Jet { double pt, y, eta, phi, m; int ncon; };

  inline const char *ClusterBackend() {
#ifdef GENSEED_USE_FASTJET
    return "FastJet";
#else
    return "tiledjet";
#endif
  }

  // E-scheme anti-kT of the Cartesian four-vectors, jets above ptmin in pT
  // order, cons[k] the input indices of jet k; the conventions of
  // genseed::Cluster (phi in (-pi, pi], y the true rapidity, m clipped at 0).
  inline std::vector<Jet> ClusterJets(const std::vector<double> &px, const std::vector<double> &py,
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

} // namespace jva

class JVA {
public:
  enum Method  { kDup = 0, kLv, kTrk, kJva, kOjet, kOpart, kNobody, kNMeth };
  enum VClass  { kAll = 0, kPV, kPU, kNClass };
  enum Outcome { kRight = 0, kWrong, kNulledSingle, kNulledCombo, kKeptCombo, kNOut };
  static const int kNCat = 5;                       // cluster-level methods in fwd/: lv trk jva ojet none
  static const int kNDjC = 2;                       // dijet classes: all pv
  static const int kNT   = 5;                       // the MHT thresholds of the tuple, 3 5 7 10 15 GeV
  static const char *MethTag(int m) { static const char *s[kNMeth] = {"dup","lv","trk","jva","ojet","opart","none"}; return (m>=0&&m<kNMeth)?s[m]:"?"; }
  static const char *ClassTag(int c) { static const char *s[kNClass] = {"all","pv","pu"}; return (c>=0&&c<kNClass)?s[c]:"?"; }
  static const char *OutTag(int o) { static const char *s[kNOut] = {"right","wrong","nulled_single","nulled_combo","kept_combo"}; return (o>=0&&o<kNOut)?s[o]:"?"; }
  static int CatMeth(int c) { static const int m[kNCat] = {kLv, kTrk, kJva, kOjet, kNobody}; return m[c]; }

  // One jet of one per-vertex event.
  struct RJet { double pt, ptcorr, eta, phi, linkz, fown, dshare; };

  JVA(TTree *tree, const char *outname = "rootfiles/JVA.root");
  ~JVA();
  void   SetOpt(const std::string &key, double val) { fOpt[key] = val; }
  double Opt(const std::string &key, double def) const {
    std::map<std::string,double>::const_iterator it = fOpt.find(key); return it == fOpt.end() ? def : it->second; }
  void   SetJEC(const char *f) { fJECFile = f; }
  void   Loop();

private:
  void Book();
  void BookTree();
  void ReadGen();                 // generated list, per-interaction ranges, true MET
  void ReconstructVertices();     // z, ownership, nearVtx (GenSeed's rules)
  void VertexWeights();           // PUPPI per usable vertex, wFwd, M_v and S_v
  void VertexRecoil();            // the vertex-resolved jets of every vertex, MHT, the recoil jva balances
  void Truth();                   // dominant interaction per cluster and per candidate
  void CheckEvidence();           // C_ku recounted, the candidates against the candMode rule
  void Methods();                 // the hosts of every method, the forward list per vertex
  void CheckTruth();              // ojet and opart recomputed from the gen side
  void FillFwd();
  void ProcessVertex(int u);      // the seven events of usable vertex u
  void FillDijet(int m, int c, const std::vector<RJet> &jets, double t1x, double t1y);
  const std::vector<jva::Jet> &GenJets(int inter);
  void FillTuple();
  void CheckIdentities();
  void PrintSummary(double seconds);
  void Write();

  std::map<std::string,double> fOpt;
  std::string fOutName, fJECFile;
  TTree *fTree; TFile *fOut; TTree *fTup;
  nanoaod::Reader *fReader;
  jec::Correction fJEC;
  puppi::Config   fPupCfg;
  genlink::Config fLinkCfg;
  jva::Config     fCfg;
  double fPtMinCluster, fPtMinJet, fT1Min, fDphiMin, fDzOracle;
  bool   fStoreTuple;

  // generated side
  genlink::GenList    fGen;
  std::vector<double> fVz, fTx, fTy;                // z and true MET per interaction
  std::vector<double> fCx, fCy, fFx, fFy;           // +generated pT at |eta| <= etaFwd, etaFwd < |eta| < 5
  std::vector<int>    fNInter, fIBase, fIEnd;       // particles and their range per interaction
  std::map<int, std::vector<jva::Jet> > fGenJets;   // pure jets per interaction, on demand
  int fNPU, fNSlot;
  // reconstructed side
  puppi::Event fPF;
  std::vector<double> fZv;
  std::vector<int>    fNtrk, fGenOf, fNearVtx, fVtxIdx, fUsable;
  int fPVIdx;                                       // usable index of key 0, -1
  genlink::Links fLinks;
  std::vector< std::vector<float> > fWv;            // PUPPI weights per usable vertex
  std::vector<float>  fWFwd;                        // vertex-blind weight, 0 at |eta| <= etaFwd
  std::vector<int>    fFwd;                         // vertex-blind candidates with wFwd > 0
  std::vector<jva::Vertex>  fVtx;                   // what jva balances (recoilMode)
  std::vector<jva::Vertex>  fPMet;                  // the particle MET M_v and S_v, whatever recoilMode
  // the vertex-resolved (none) event of every usable vertex: its candidates,
  // jets with their constituents (input positions), MET; and per threshold
  // kT[t] of the tuple, [u*kNT + t]: HT, sum pT^2, MHT of its jets
  std::vector< std::vector<int> > fBase;
  std::vector< std::vector<jva::Jet> > fBaseJets;
  std::vector< std::vector< std::vector<int> > > fBaseCons;
  std::vector<double> fBaseMx, fBaseMy, fHtT, fHtSqT, fMhtXT, fMhtYT;
  std::vector<jva::Cluster> fCl;
  std::vector< std::vector<int> > fClMem;
  std::vector<int>    fClOf;                        // cluster of every PF candidate, -1
  std::vector<int>    fClDom, fClTrue; std::vector<double> fClShare;
  std::vector<int>    fClTop3; std::vector<double> fClTop3Share;   // 3 per cluster
  std::vector<int>    fPartDom, fPartHost;          // per PF candidate (opart)
  std::vector<int>    fHost[kNMeth];                // per cluster (lv trk jva ojet)
  std::vector< std::vector<int> > fFwdTo[kNMeth];   // forward candidates given to each usable vertex
  jva::Stats fStats;
  double fMetPV[2][2]; bool fHavePV[2]; int fNJetPV[2];   // dup and lv at the PV, for the check

  // run totals
  Long64_t fNEvents, fNVtx, fNOwned, fNCl, fNClTrk, fNComp, fNExact, fNSweeps, fNFwdCand, fNNoPV, fLargest, fNCand, fNClusterings;
  double fAssignMs, fObjSum, fSig2Sum, fRecoilS2;
  // single-interaction clusters above 5 GeV: count; true vertex among the
  // candidates, leading the trk ranking; among / leading the T and C evidence
  double fS5, fS5In, fS5Top, fS5InT, fS5TopT, fS5InC, fS5TopC;
  double fMetN[kNMeth][kNClass], fMetS2[kNMeth][kNClass], fResN[kNMeth][kNClass], fResS2[kNMeth][kNClass];
  double fTrueN[kNClass], fTrueS2[kNClass];
  double fIdDup, fIdOnce, fIdCover, fIdWfwd, fIdPV, fIdTruth, fIdNone, fIdCh, fIdCand, fIdGen;

  // tuple branches
  UInt_t bRun, bLumi; ULong64_t bEvent; Int_t bNPU, bNVtx, bNComp, bNExact; Float_t bObj;
  std::vector<int>   bVtxKey, bVtxOwner;
  std::vector<float> bVtxZ, bVtxMx, bVtxMy, bVtxS;
  std::vector<float> bVtxHt[kNT], bVtxHtSq[kNT], bVtxMhtX[kNT], bVtxMhtY[kNT], bVtxGenMhtX[kNT], bVtxGenMhtY[kNT];
  std::vector<float> bVtxGenMetX, bVtxGenMetY, bVtxGenFwdX, bVtxGenFwdY, bVtxGenCenX, bVtxGenCenY;
  std::vector<float> bClPt, bClEta, bClPhi, bClPx, bClPy, bClShare, bCandTrkPt, bChPt, bClTop3Share;
  std::vector<int>   bClNcand, bClOff, bCandVtx, bClDom, bClTrue, bClLv, bClTrk, bClJva, bClOjet, bClNone;
  std::vector<int>   bClNopt, bClNch, bClChOff, bChVtx, bClTop3;

  // histograms
  TH1D *hcount, *hopts;
  TH1D *hmetxy[kNMeth][kNClass], *hmet[kNMeth][kNClass], *hmetres[kNMeth][kNClass];
  TH2D *hmetxy_npu[kNMeth], *hmetxy_s[kNMeth], *hmetres_s[kNMeth];
  TH1D *hjeteta[kNMeth+1][3], *hjetpt[kNMeth+1][3], *hdom[kNMeth][3];   // [kNMeth] = gen
  TProfile *pown[kNMeth];
  TH1D *hclpt[2], *hcand; TH2D *hshare, *hcat[kNCat][2];               // [0] trk, [1] notrk
  TH2F *hadb[kNMeth][kNDjC][jva::kNEtaDj][2], *hampf[kNMeth][kNDjC][jva::kNEtaDj][2], *hrtrue[kNMeth][kNDjC][jva::kNEtaDj][2];
  TH2F *hdbpar[kNMeth][kNDjC][jva::kNEtaDj], *hdbperp[kNMeth][kNDjC][jva::kNEtaDj];
  TH2F *hmpfpar[kNMeth][kNDjC][jva::kNEtaDj], *hmpfperp[kNMeth][kNDjC][jva::kNEtaDj], *hrlink[kNMeth][kNDjC][jva::kNEtaDj];
  TH1D *hsel[kNMeth][kNDjC];
};

#endif
