// genseed - Copyright 2026 Mikko Voutilainen - Apache License 2.0
// JVA.C - see JVA.h for what this measures and why, jvassoc.h for the
// association itself.
//
//   root -l -b -q 'runJVA.C("files.txt","_test",0,1,300)'      // 300 crossings, text/jec_v1.txt
//   ./runjva.sh v1 8                                           // everything, 8 jobs
//
// The per-crossing order is the one the dependencies dictate: the generated
// list and the true MET of every interaction; the PF event and
// puppi::Prepare; the vertices and their ownership; the linker (once, vertex
// independent); the PUPPI weights of every usable vertex, and from them wFwd,
// M_v and S_v; the forward clusters; the truth of every cluster and every
// vertex-blind candidate; the hosts of every method and the forward list of
// every vertex; the vertex-resolved event of every vertex, clustered once,
// and from its jets the recoil jva balances; the forward clusters and their
// evidence; the truth of every cluster and every vertex-blind candidate; the
// hosts of every method and the forward list of every vertex; then each
// usable vertex in turn with all seven methods - the event, its jets and
// their truth, its MET, every fill.  Nothing is deduplicated across vertices,
// unlike GenSeed: every per-vertex event is filled once, and the duplication
// of the forward jets is what dup shows.  A method that gives a vertex no
// forward candidate gives it the vertex-resolved event, whose jets are
// already there (the same input in the same order: the same jets).
#include "JVA.h"

#include <TStopwatch.h>
#include <TSystem.h>
#include <TDirectory.h>
#include <TNamed.h>
#include <algorithm>
#include <cstdio>
#include <cstring>

using namespace jva;

namespace {
  inline double dphi(double a, double b) {
    double d = a - b;
    if (d > M_PI) d -= 2*M_PI; else if (d < -M_PI) d += 2*M_PI;
    return d;
  }
  const double kEtaMaxPart = 5.5;    // generated particles kept: GenSeed's etaMaxPart
  const double kEtaTrue    = 5.0;    // MET_true: what a perfect |eta| < 5 detector would see
  const double kEtaJetMax  = 5.191;  // dijet jets
  const double kEtaTag     = 1.305;  // tag jet
  const double kAlphaMax   = 0.3;
  const double kPtX[3]     = {5, 10, 20};
  const double kShareMax   = 1. - 1e-9;   // a share of exactly 1 lands in the last bin of a 0..1 axis
  const double kSBins[]    = {0,3,5,7,10,15,20,30,40,55,70,90,110,140,200,300};
  const int    kNS         = sizeof(kSBins)/sizeof(kSBins[0]) - 1;
  inline int Reg(double aeta) { return aeta < 2.5 ? 0 : (aeta < 3.0 ? 1 : (aeta < 5.0 ? 2 : -1)); }
  const char *kRegLabel[3] = {"|#eta| < 2.5","2.5 < |#eta| < 3.0","3.0 < |#eta| < 5.0"};
  // the knobs written to hist/hopts, in bin order; v2's appended, so that v1's
  // bins stay where they were
  const char *kKnob[] = {"R","ptMinCluster","etaFwd","trkMinPt","sigma0","sigmaK","tauTrk","tauAll",
    "maxCombos","maxSweeps","ptMinJet","t1Min","dphiMin","dzVertexGen","dzOracle","minTrkPerVertex",
    "storeTuple","jecLoaded","firstEntry","lastEntry","progressEvery","useFastJet",
    "candMode","chMinPt","recoilMode","mhtMin","sigmaJ",
    "pupEtaTracker","pupEtaVtxAssoc","pupFloor1Pt","pupFloor1Slope"};
  const int kNKnob = sizeof(kKnob)/sizeof(kKnob[0]);
  const char *kSel[] = {"vertices",">=2 jets","dphi","tag","alpha<0.3","samebin"};
  const int kNSel = sizeof(kSel)/sizeof(kSel[0]);
  const int kT[JVA::kNT] = {3, 5, 7, 10, 15};   // GeV, the MHT thresholds of the tuple
  const double kVtxEta = 2.5;                    // tracks for the vertex z (puppi::Config's default etaTracker)
  template <class T> void Vec(TTree *t, const char *n, std::vector<T> &v) { t->Branch(n, &v); }

  // One per-vertex event as the clustering takes it: the vertex-resolved
  // candidates at w_i(v), then the vertex-blind ones given to v at wFwd; the
  // PF index and weight of every input, and the MET, -(sum of the input pT)
  // in input order.  The one place an event is built, so that the cached
  // vertex-resolved event and a rebuilt one are the same numbers.
  struct Ev {
    std::vector<double> px, py, pz, E; std::vector<int> idx; std::vector<float> wt; double mx = 0, my = 0;
    void build(const puppi::Event &pf, const std::vector<int> &base, const std::vector<float> &w,
               const std::vector<int> &fw, const std::vector<float> &wf, bool kin) {
      px.clear(); py.clear(); pz.clear(); E.clear(); idx.clear(); wt.clear(); mx = my = 0;
      auto push = [&](int i, float wi) {
        idx.push_back(i); wt.push_back(wi);
        if (!kin) return;
        const double s = wi*pf.pt[i], mm = wi*pf.mass[i];
        const double x = s*cos(pf.phi[i]), y = s*sin(pf.phi[i]), z = s*sinh(pf.eta[i]);
        px.push_back(x); py.push_back(y); pz.push_back(z); E.push_back(sqrt(x*x + y*y + z*z + mm*mm));
      };
      for (size_t j = 0; j < base.size(); ++j) push(base[j], w[base[j]]);
      for (size_t j = 0; j < fw.size(); ++j) push(fw[j], wf[fw[j]]);
      for (size_t n = 0; n < px.size(); ++n) { mx -= px[n]; my -= py[n]; }
    }
  };

  // Linked generated pT by interaction slot: dense over the slots, the touched
  // ones listed, so that clearing costs what was added.
  struct SlotSum {
    std::vector<double> v; std::vector<char> on; std::vector<int> touched; double total = 0;
    void reset(int n) { v.assign(n, 0.); on.assign(n, 0); touched.clear(); total = 0; }
    void clear() { for (size_t j = 0; j < touched.size(); ++j) { v[touched[j]] = 0.; on[touched[j]] = 0; } touched.clear(); total = 0; }
    void add(int s, double x) { if (!on[s]) { on[s] = 1; touched.push_back(s); } v[s] += x; total += x; }
    // the dominant slot (the lowest on an exact tie) and its share; -1 and 0 if nothing
    int dominant(double &share) const {
      int top = -1; double tv = 0;
      for (size_t j = 0; j < touched.size(); ++j) {
        const int s = touched[j];
        if (v[s] > tv || (v[s] == tv && top >= 0 && s < top)) { tv = v[s]; top = s; }
      }
      share = (top >= 0 && total > 0) ? tv/total : 0.;
      return top;
    }
    // the n leading slots by the same order (the first is dominant()), -1 and
    // 0 where fewer than n slots have anything
    void leading(int n, int *slot, double *share) const {
      for (int j = 0; j < n; ++j) { slot[j] = -1; share[j] = 0.; }
      for (size_t j = 0; j < touched.size(); ++j) {
        const int s = touched[j];
        if (!(v[s] > 0)) continue;
        int p = n;                                   // insertion into the sorted top n
        while (p > 0 && (slot[p-1] < 0 || v[s] > v[slot[p-1]] || (v[s] == v[slot[p-1]] && s < slot[p-1]))) --p;
        if (p >= n) continue;
        for (int q = n - 1; q > p; --q) slot[q] = slot[q-1];
        slot[p] = s;
      }
      for (int j = 0; j < n; ++j) share[j] = (slot[j] >= 0 && total > 0) ? v[slot[j]]/total : 0.;
    }
  };
  // the same on a dense row, slots ascending
  inline int DominantRow(const double *row, int n, double &share) {
    int top = -1; double tv = 0, tot = 0;
    for (int s = 0; s < n; ++s) { tot += row[s]; if (row[s] > tv) { tv = row[s]; top = s; } }
    share = (top >= 0 && tot > 0) ? tv/tot : 0.;
    return top;
  }
}

JVA::JVA(TTree *tree, const char *outname)
  : fOutName(outname), fJECFile(""), fTree(tree), fOut(0), fTup(0), fReader(0),
    fPtMinCluster(1.), fPtMinJet(3.), fT1Min(3.), fDphiMin(2.7), fDzOracle(0.2), fStoreTuple(true),
    fNPU(0), fNSlot(0), fPVIdx(-1),
    fNEvents(0), fNVtx(0), fNOwned(0), fNCl(0), fNClTrk(0), fNComp(0), fNExact(0), fNSweeps(0),
    fNFwdCand(0), fNNoPV(0), fLargest(0), fNCand(0), fNClusterings(0), fAssignMs(0), fObjSum(0), fSig2Sum(0), fRecoilS2(0),
    fS5(0), fS5In(0), fS5Top(0), fS5InT(0), fS5TopT(0), fS5InC(0), fS5TopC(0),
    fIdDup(0), fIdOnce(0), fIdCover(0), fIdWfwd(0), fIdPV(0), fIdTruth(0), fIdNone(0), fIdCh(0), fIdCand(0), fIdGen(0)
{
  for (int m = 0; m < kNMeth; ++m) for (int c = 0; c < kNClass; ++c)
    fMetN[m][c] = fMetS2[m][c] = fResN[m][c] = fResS2[m][c] = 0;
  for (int c = 0; c < kNClass; ++c) fTrueN[c] = fTrueS2[c] = 0;
}
JVA::~JVA() { delete fReader; }

////////////////////////////////////////////////////////////////////////////
void JVA::Book()
{
  fOut = new TFile(fOutName.c_str(), "RECREATE");
  if (fOut->IsZombie()) { printf("JVA: cannot open %s\n", fOutName.c_str()); delete fOut; fOut = 0; return; }
  const int nb = kNPtExt; const double *xb = kPtBinsExt;
  const char *T = fJEC.IsLoaded() ? "corr" : "raw";   // pT^corr is pT^raw without a table

  TDirectory *dh = fOut->mkdir("hist"); dh->cd();
  hcount = new TH1D("hcount","counters",8,0,8);
  const char *cl[8] = {"events_a","events_b","npu_a","npu_b","entries","jobs","vertices","owned"};
  for (int i = 0; i < 8; ++i) hcount->GetXaxis()->SetBinLabel(i+1, cl[i]);
  hopts = new TH1D("hopts","knobs, value x jobs",kNKnob,0,kNKnob);
  for (int i = 0; i < kNKnob; ++i) hopts->GetXaxis()->SetBinLabel(i+1, kKnob[i]);
  // Text next to hopts; TNamed has no Merge, so hadd keeps a cycle per job.
  TNamed("puppiVersion", puppi::kVersion).Write();
  TNamed("jecFile", fJEC.IsLoaded() ? fJECFile.c_str() : "none").Write();
  TNamed("backend", ClusterBackend()).Write();

  // ---- MET ------------------------------------------------------------------
  TDirectory *dmet = fOut->mkdir("met");
  for (int m = 0; m < kNMeth; ++m) {
    TDirectory *d = dmet->mkdir(MethTag(m)); d->cd();
    for (int c = 0; c < kNClass; ++c) {
      hmetxy[m][c]  = new TH1D(Form("hmetxy_%s",ClassTag(c)), Form("%s, %s vertices;MET_{x}, MET_{y} [GeV];components",MethTag(m),ClassTag(c)),400,-100,100);
      hmet[m][c]    = new TH1D(Form("hmet_%s",ClassTag(c)),   Form("%s, %s vertices;|MET| [GeV];vertices",MethTag(m),ClassTag(c)),200,0,200);
      hmetres[m][c] = new TH1D(Form("hmetres_%s",ClassTag(c)),Form("%s, %s owned vertices;MET - MET_{true}(owner), both components [GeV];components",MethTag(m),ClassTag(c)),400,-100,100);
    }
    hmetxy_npu[m] = new TH2D("hmetxy_npu",Form("%s, all vertices;N_{PU};MET component [GeV]",MethTag(m)),100,0,100,200,-100,100);
    hmetxy_s[m]   = new TH2D("hmetxy_s",  Form("%s, all vertices;S_{v} [GeV];MET component [GeV]",MethTag(m)),kNS,kSBins,200,-100,100);
    hmetres_s[m]  = new TH2D("hmetres_s", Form("%s, owned vertices;S_{v} [GeV];MET - MET_{true}(owner) component [GeV]",MethTag(m)),kNS,kSBins,200,-100,100);
  }

  // ---- the forward spike: jets at owned vertices, and the gen reference -------
  TDirectory *dj = fOut->mkdir("jets");
  for (int m = 0; m <= kNMeth; ++m) {
    const bool gen = (m == kNMeth);
    const char *mt = gen ? "gen" : MethTag(m);
    TDirectory *d = dj->mkdir(mt); d->cd();
    for (int x = 0; x < 3; ++x)
      hjeteta[m][x] = new TH1D(Form("hjeteta_pt%.0f",kPtX[x]),
        gen ? Form("pure GenJets of the owner, p_{T} > %.0f GeV;|#eta|;jets",kPtX[x])
            : Form("%s, owned vertices, p_{T}^{raw} > %.0f GeV;|#eta|;jets",mt,kPtX[x]),50,0,5);
    for (int r = 0; r < 3; ++r)
      hjetpt[m][r] = new TH1D(Form("hjetpt_r%d",r),
        gen ? Form("pure GenJets of the owner, %s;p_{T} [GeV];jets",kRegLabel[r])
            : Form("%s, owned vertices, %s;p_{T}^{%s} [GeV];jets",mt,kRegLabel[r],T),nb,xb);
    if (gen) continue;
    pown[m] = new TProfile("pown_eta",Form("%s, owned vertices, p_{T}^{raw} > 5 GeV;|#eta|;owner fraction of the linked generated p_{T}",mt),50,0,5);
    for (int r = 0; r < 3; ++r)
      hdom[m][r] = new TH1D(Form("hdom_r%d",r),Form("%s, owned vertices, p_{T}^{raw} > 5 GeV, %s;dominant-interaction share of the linked p_{T};jets",mt,kRegLabel[r]),20,0,1);
  }

  // ---- the forward clusters and what every method did with them ----------------
  TDirectory *df = fOut->mkdir("fwd"); df->cd();
  const char *cw = fCfg.candMode == 1 ? "in-cluster charged" : (fCfg.candMode == 2 ? "track or in-cluster charged" : "track");
  hclpt[0] = new TH1D("hclpt_trk",  Form("forward clusters with %s candidates;p_{T}^{cluster} [GeV];clusters",cw),nb,xb);
  hclpt[1] = new TH1D("hclpt_notrk",Form("forward clusters without %s candidates;p_{T}^{cluster} [GeV];clusters",cw),nb,xb);
  hshare = new TH2D("hshare","forward clusters;p_{T}^{cluster} [GeV];dominant-interaction share (0: unlinked)",nb,xb,20,0,1);
  hcand  = new TH1D("hcand",Form("forward clusters;%s candidates;clusters",cw),21,-0.5,20.5);
  for (int c = 0; c < kNCat; ++c) {
    TDirectory *d = df->mkdir(MethTag(CatMeth(c))); d->cd();
    for (int t = 0; t < 2; ++t) {
      hcat[c][t] = new TH2D(t == 0 ? "hcat_trk" : "hcat_notrk",
        Form("%s, clusters %s %s candidates;p_{T}^{cluster} [GeV];outcome",MethTag(CatMeth(c)),t == 0 ? "with" : "without",cw),
        nb,xb,kNOut,0,kNOut);
      for (int o = 0; o < kNOut; ++o) hcat[c][t]->GetYaxis()->SetBinLabel(o+1, OutTag(o));
    }
  }

  // ---- dijets from every vertex -------------------------------------------------
  TDirectory *dd = fOut->mkdir("dijet");
  for (int m = 0; m < kNMeth; ++m) {
    TDirectory *dm = dd->mkdir(MethTag(m));
    for (int c = 0; c < kNDjC; ++c) {
      TDirectory *d = dm->mkdir(ClassTag(c)); d->cd();
      for (int ie = 0; ie < kNEtaDj; ++ie) {
        const TString el = Form("%s %s, %.3f < |#eta| < %.3f", MethTag(m), ClassTag(c), kEtaDj[ie], kEtaDj[ie+1]);
        for (int ia = 0; ia < 2; ++ia) {
          const char *al = ia == 0 ? "#alpha < 0.3" : "no #alpha cut";
          hadb[m][c][ie][ia]   = new TH2F(Form("hadb_e%02d_a%d",ie,ia),  Form("%s (probe), %s;p_{T}^{avg} [GeV];A_{DB}",el.Data(),al),nb,xb,100,-1,1);
          hampf[m][c][ie][ia]  = new TH2F(Form("hampf_e%02d_a%d",ie,ia), Form("%s (probe), %s;p_{T}^{avg} [GeV];A_{MPF}",el.Data(),al),nb,xb,100,-1,1);
          hrtrue[m][c][ie][ia] = new TH2F(Form("hrtrue_e%02d_a%d",ie,ia),Form("%s (probe), %s;p_{T}^{avg} [GeV];r_{true} = (p_{T}^{%s}/p_{T}^{linkz})_{probe} / (p_{T}^{%s}/p_{T}^{linkz})_{tag}",el.Data(),al,T,T),nb,xb,100,0,3);
        }
        hdbpar[m][c][ie]   = new TH2F(Form("hdbpar_e%02d",ie),  Form("%s, both jets, #alpha < 0.3;p_{T}^{avg} [GeV];(p_{1}+p_{2}) #upoint n / p_{T}^{avg}",el.Data()),nb,xb,100,-1,1);
        hdbperp[m][c][ie]  = new TH2F(Form("hdbperp_e%02d",ie), Form("%s, both jets, #alpha < 0.3;p_{T}^{avg} [GeV];(p_{1}+p_{2}) #upoint n_{#perp} / p_{T}^{avg}",el.Data()),nb,xb,100,-1,1);
        hmpfpar[m][c][ie]  = new TH2F(Form("hmpfpar_e%02d",ie), Form("%s, both jets, #alpha < 0.3;p_{T}^{avg} [GeV];MET_{T1} #upoint n / p_{T}^{avg}",el.Data()),nb,xb,100,-1,1);
        hmpfperp[m][c][ie] = new TH2F(Form("hmpfperp_e%02d",ie),Form("%s, both jets, #alpha < 0.3;p_{T}^{avg} [GeV];MET_{T1} #upoint n_{#perp} / p_{T}^{avg}",el.Data()),nb,xb,100,-1,1);
        hrlink[m][c][ie]   = new TH2F(Form("hrlink_e%02d",ie),  Form("%s, both jets, #alpha < 0.3;p_{T}^{avg} [GeV];p_{T}^{%s} / p_{T}^{linkz} (both jets)",el.Data(),T),nb,xb,150,0,3);
      }
      hsel[m][c] = new TH1D("hsel",Form("%s %s vertices;selection;vertex events",MethTag(m),ClassTag(c)),kNSel,0,kNSel);
      for (int s = 0; s < kNSel; ++s) hsel[m][c]->GetXaxis()->SetBinLabel(s+1, kSel[s]);
    }
  }
  fOut->cd();
}

// The tuple: one entry per crossing, enough to re-run the association offline
// with other knobs (M_v, S_v, the MHT of the vertex-resolved jets at five
// thresholds, the clusters and their evidence of both kinds) and to judge it
// (the generated recoil of every vertex's owner, the dominant interactions and
// shares, true vertex, the host of every method).  The branches are listed in
// JVA.h.  Indices are into the usable-vertex list of the same entry; -1 is null.
void JVA::BookTree()
{
  fOut->cd();
  fTup = new TTree("jva","jva: vertices and their recoil, forward clusters and their evidence, truth, hosts");
  fTup->Branch("run",&bRun,"run/i"); fTup->Branch("lumi",&bLumi,"lumi/i");
  fTup->Branch("event",&bEvent,"event/l"); fTup->Branch("nPU",&bNPU,"nPU/I");
  fTup->Branch("nVtxUsable",&bNVtx,"nVtxUsable/I");
  fTup->Branch("jva_ncomp",&bNComp,"jva_ncomp/I"); fTup->Branch("jva_nexact",&bNExact,"jva_nexact/I");
  fTup->Branch("jva_objective",&bObj,"jva_objective/F");
  Vec(fTup,"vtx_key",bVtxKey); Vec(fTup,"vtx_z",bVtxZ); Vec(fTup,"vtx_owner",bVtxOwner);
  Vec(fTup,"vtx_mx",bVtxMx); Vec(fTup,"vtx_my",bVtxMy); Vec(fTup,"vtx_s",bVtxS);
  for (int t = 0; t < kNT; ++t) {
    Vec(fTup,Form("vtx_ht%d",kT[t]),bVtxHt[t]); Vec(fTup,Form("vtx_htsq%d",kT[t]),bVtxHtSq[t]);
    Vec(fTup,Form("vtx_mht%d_x",kT[t]),bVtxMhtX[t]); Vec(fTup,Form("vtx_mht%d_y",kT[t]),bVtxMhtY[t]);
  }
  Vec(fTup,"vtx_genmet_x",bVtxGenMetX); Vec(fTup,"vtx_genmet_y",bVtxGenMetY);
  Vec(fTup,"vtx_genfwd_x",bVtxGenFwdX); Vec(fTup,"vtx_genfwd_y",bVtxGenFwdY);
  Vec(fTup,"vtx_gencen_x",bVtxGenCenX); Vec(fTup,"vtx_gencen_y",bVtxGenCenY);
  for (int t = 0; t < kNT; ++t) {
    Vec(fTup,Form("vtx_genmht%d_x",kT[t]),bVtxGenMhtX[t]); Vec(fTup,Form("vtx_genmht%d_y",kT[t]),bVtxGenMhtY[t]); }
  Vec(fTup,"cl_pt",bClPt); Vec(fTup,"cl_eta",bClEta); Vec(fTup,"cl_phi",bClPhi);
  Vec(fTup,"cl_px",bClPx); Vec(fTup,"cl_py",bClPy); Vec(fTup,"cl_nopt",bClNopt);
  Vec(fTup,"cl_ncand",bClNcand); Vec(fTup,"cl_candoff",bClOff);
  Vec(fTup,"cand_vtx",bCandVtx); Vec(fTup,"cand_trkpt",bCandTrkPt);
  Vec(fTup,"cl_nch",bClNch); Vec(fTup,"cl_choff",bClChOff);
  Vec(fTup,"ch_vtx",bChVtx); Vec(fTup,"ch_pt",bChPt);
  Vec(fTup,"cl_dom",bClDom); Vec(fTup,"cl_share",bClShare); Vec(fTup,"cl_true",bClTrue);
  Vec(fTup,"cl_top3_slot",bClTop3); Vec(fTup,"cl_top3_share",bClTop3Share);
  Vec(fTup,"cl_host_lv",bClLv); Vec(fTup,"cl_host_trk",bClTrk); Vec(fTup,"cl_host_jva",bClJva); Vec(fTup,"cl_host_ojet",bClOjet);
  Vec(fTup,"cl_host_none",bClNone);
}

////////////////////////////////////////////////////////////////////////////
// The generated list (each interaction contiguous), the range of every
// interaction, and its true MET: -(sum of generated pT with |eta| < 5); and
// the two parts of that sum with their own sign, +(sum at |eta| <= etaFwd)
// and +(sum at etaFwd < |eta| < 5), so that MET_true = -(central + forward).
void JVA::ReadGen()
{
  fNPU = fReader->NPU();
  fNSlot = fReader->FillGen(fGen, fVz, kEtaMaxPart, true);
  const int ng = (int)fGen.size();
  fNInter.assign(fNSlot, 0); fIBase.assign(fNSlot, 0); fIEnd.assign(fNSlot, 0);
  fTx.assign(fNSlot, 0.); fTy.assign(fNSlot, 0.);
  fCx.assign(fNSlot, 0.); fCy.assign(fNSlot, 0.); fFx.assign(fNSlot, 0.); fFy.assign(fNSlot, 0.);
  for (int g = 0; g < ng; ++g) {
    const int i = fGen.inter[g];
    if (fNInter[i]++ == 0) fIBase[i] = g;
    fIEnd[i] = g + 1;
    const double ae = fabs(fGen.eta[g]);
    if (ae < kEtaTrue) {
      // v1's statements as they were (a contracted multiply-subtract must not
      // move MET_true by an ulp), the two parts beside them
      fTx[i] -= fGen.pt[g]*cos(fGen.phi[g]); fTy[i] -= fGen.pt[g]*sin(fGen.phi[g]);
      const double x = fGen.pt[g]*cos(fGen.phi[g]), y = fGen.pt[g]*sin(fGen.phi[g]);
      if (ae <= fCfg.etaFwd) { fCx[i] += x; fCy[i] += y; } else { fFx[i] += x; fFy[i] += y; }
    }
  }
  for (int i = 0; i < fNSlot; ++i)
    fIdGen = std::max(fIdGen, std::max(fabs(fTx[i] + fCx[i] + fFx[i]), fabs(fTy[i] + fCy[i] + fFy[i])));
  fGenJets.clear();
}

// GenSeed::ReconstructVertices, rule for rule: z per reco vertex, the
// interaction each one owns (unique per vertex, nearest within dzVertexGen,
// among the interactions with a particle in acceptance), the usable vertex
// nearest to every interaction; plus the usable index of the PV (key 0).
void JVA::ReconstructVertices()
{
  // Vertex z from the tracks inside |eta| < 2.5 whatever PUPPI's tracker edge
  // is set to (pupEtaTracker), so that the vertices, their owners and the
  // truth stay the same when the PUPPI knobs are varied.
  puppi::VertexZ(fPF.ztrk, fPF.vref, fPF.charged, fPF.eta, kVtxEta,
                 int(Opt("minTrkPerVertex",3)), fZv, fNtrk);
  const int nv = (int)fZv.size();
  fGenOf.assign(nv, -1); fVtxIdx.assign(nv, -1); fUsable.clear();
  const double dzmax = Opt("dzVertexGen",0.05);
  for (int v = 0; v < nv; ++v) {
    if (fZv[v] > 1e8) continue;
    fVtxIdx[v] = (int)fUsable.size(); fUsable.push_back(v);
    double best = 1e9; int bi = -1;
    for (int i = 0; i < fNSlot; ++i) {
      if (!fNInter[i]) continue;
      const double d = fabs(fZv[v] - fVz[i]);
      if (d < best) { best = d; bi = i; }
    }
    if (bi >= 0 && best < dzmax) fGenOf[v] = bi;
  }
  fNearVtx.assign(fNSlot, -1);
  for (int i = 0; i < fNSlot; ++i) {
    double best = 1e9;
    for (size_t u = 0; u < fUsable.size(); ++u) {
      const int v = fUsable[u]; const double d = fabs(fZv[v] - fVz[i]);
      if (d < best) { best = d; fNearVtx[i] = v; }
    }
  }
  fPVIdx = (nv > 0) ? fVtxIdx[0] : -1;
}

// PUPPI at every usable vertex; the vertex-resolved MET and scalar sum of
// each; wFwd from the PV's weights (the first usable vertex if the PV is not
// usable) and the check that every vertex gives the same.
void JVA::VertexWeights()
{
  const int nu = (int)fUsable.size(), np = (int)fPF.size();
  const double ef = fCfg.etaFwd;
  if ((int)fWv.size() < nu) fWv.resize(nu);
  fVtx.assign(nu, Vertex());
  for (int u = 0; u < nu; ++u) {
    const int key = fUsable[u];
    puppi::Weights(fPF, key, fZv[key], fPupCfg, fWv[u]);
    Vertex &V = fVtx[u];
    const std::vector<float> &w = fWv[u];
    for (int i = 0; i < np; ++i) {
      if (w[i] <= 0 || fabs(fPF.eta[i]) > ef) continue;
      const double p = w[i]*fPF.pt[i];
      V.mx -= p*cos(fPF.phi[i]); V.my -= p*sin(fPF.phi[i]); V.s += p;
    }
  }
  fWFwd.assign(np, 0.f); fFwd.clear();
  if (nu == 0) return;
  const int ref = fPVIdx >= 0 ? fPVIdx : 0;
  for (int i = 0; i < np; ++i) {
    if (fabs(fPF.eta[i]) <= ef) continue;
    fWFwd[i] = fWv[ref][i];
    if (fWFwd[i] > 0) fFwd.push_back(i);
    for (int u = 0; u < nu; ++u) fIdWfwd = std::max(fIdWfwd, (double)fabs(fWv[u][i] - fWFwd[i]));
  }
}

// The vertex-resolved event of every usable vertex (the none event: its
// candidates at w_i(v) and nothing else), clustered once: its jets are reused
// by every method that gives the vertex no forward candidate (ProcessVertex),
// they make the HT, sum pT^2 and MHT at the tuple's thresholds, and in
// recoilMode 1 the recoil jva balances, the MHT above mhtMin with its
// sigma_v terms (jvassoc.h, THE RECOIL).  fPMet keeps the particle MET M_v
// whatever recoilMode; in recoilMode 0 fVtx is that, untouched.
void JVA::VertexRecoil()
{
  const int nu = (int)fUsable.size(), np = (int)fPF.size();
  fPMet = fVtx;
  if ((int)fBase.size() < nu) { fBase.resize(nu); fBaseJets.resize(nu); fBaseCons.resize(nu); }
  fBaseMx.assign(nu, 0.); fBaseMy.assign(nu, 0.);
  fHtT.assign((size_t)nu*kNT, 0.); fHtSqT.assign((size_t)nu*kNT, 0.);
  fMhtXT.assign((size_t)nu*kNT, 0.); fMhtYT.assign((size_t)nu*kNT, 0.);
  static const std::vector<int> kNoFwd;
  Ev ev; std::vector<double> jpt, jphi;
  for (int u = 0; u < nu; ++u) {
    const std::vector<float> &w = fWv[u];
    std::vector<int> &base = fBase[u]; base.clear();
    for (int i = 0; i < np; ++i) if (w[i] > 0 && fabs(fPF.eta[i]) <= fCfg.etaFwd) base.push_back(i);
    ev.build(fPF, base, w, kNoFwd, fWFwd, true);
    fBaseMx[u] = ev.mx; fBaseMy[u] = ev.my;
    fBaseJets[u] = ClusterJets(ev.px, ev.py, ev.pz, ev.E, fCfg.R, fPtMinCluster, fBaseCons[u]);
    ++fNClusterings;
    const std::vector<Jet> &js = fBaseJets[u];
    jpt.resize(js.size()); jphi.resize(js.size());
    for (size_t k = 0; k < js.size(); ++k) { jpt[k] = js[k].pt; jphi[k] = js[k].phi; }
    for (int t = 0; t < kNT; ++t) {                   // the same code as the recoil below
      jva::Vertex x; jva::SetMHT(x, jpt, jphi, kT[t]);
      const size_t n = (size_t)u*kNT + t;
      fHtT[n] = x.ht; fHtSqT[n] = x.sj2; fMhtXT[n] = x.mx; fMhtYT[n] = x.my;
    }
    if (fCfg.recoilMode == 1) jva::SetMHT(fVtx[u], jpt, jphi, fCfg.mhtMin);
  }
  for (int u = 0; u < nu; ++u) {
    fSig2Sum += jva::Sigma2(fVtx[u], fCfg);
    fRecoilS2 += fVtx[u].mx*fVtx[u].mx + fVtx[u].my*fVtx[u].my;
  }
}

// The dominant interaction of every forward cluster and of every vertex-blind
// candidate, from the PF side of the links: linked generated pT wFwd_i f_gi
// pT_g summed by slot (the candidate's own sum without the common wFwd_i),
// share = top/total; the true vertex is nearVtx of the dominant interaction if
// the share is above 0.5 and the vertex within dzOracle of it, else -1 (also
// for an unlinked cluster or candidate).
void JVA::Truth()
{
  const int nk = (int)fCl.size(), np = (int)fPF.size();
  const genlink::Links &L = fLinks;
  fClOf.assign(np, -1);
  for (int k = 0; k < nk; ++k)
    for (size_t j = 0; j < fClMem[k].size(); ++j) {
      const int i = fClMem[k][j];
      if (fClOf[i] >= 0) fIdCover += 1;          // in two clusters
      fClOf[i] = k;
    }
  for (int i = 0; i < np; ++i) {
    const bool blind = fabs(fPF.eta[i]) > fCfg.etaFwd && fWFwd[i] > 0;
    if (blind != (fClOf[i] >= 0)) fIdCover += 1;  // in none, or one that should not be
  }
  fClDom.assign(nk, -1); fClShare.assign(nk, 0.); fClTrue.assign(nk, -1);
  fClTop3.assign((size_t)3*nk, -1); fClTop3Share.assign((size_t)3*nk, 0.);
  fPartDom.assign(np, -1); fPartHost.assign(np, -1);
  auto trueVtx = [&](int inter, double share) {
    if (inter < 0 || !(share > 0.5)) return -1;
    const int v = fNearVtx[inter];
    if (v < 0 || !(fabs(fZv[v] - fVz[inter]) < fDzOracle)) return -1;
    return fVtxIdx[v];
  };
  SlotSum cs, ps; cs.reset(fNSlot); ps.reset(fNSlot);
  for (int k = 0; k < nk; ++k) {
    cs.clear();
    for (size_t j = 0; j < fClMem[k].size(); ++j) {
      const int i = fClMem[k][j];
      ps.clear();
      for (int n = L.pfirst[i]; n < L.pfirst[i] + L.pcount[i]; ++n) {
        const int g = L.pgen[n]; const double x = L.pfrac[n]*fGen.pt[g];
        ps.add(fGen.inter[g], x); cs.add(fGen.inter[g], fWFwd[i]*x);
      }
      double sh; const int d = ps.dominant(sh);
      fPartDom[i] = d; fPartHost[i] = trueVtx(d, sh);
    }
    double sh; const int d = cs.dominant(sh);
    fClDom[k] = d; fClShare[k] = sh; fClTrue[k] = trueVtx(d, sh);
    cs.leading(3, &fClTop3[(size_t)3*k], &fClTop3Share[(size_t)3*k]);
  }
}

// The evidence of every cluster, checked from outside jvassoc.h: C_ku summed
// over u against the raw pT of its charged constituents on usable vertices,
// recounted from the members; T_ku and C_ku of every candidate against the
// evidence lists, and the candidates against the candMode rule applied to
// them (a vertex in either list that the rule takes and cand lacks, or the
// reverse, is a mismatch).
void JVA::CheckEvidence()
{
  const int nu = (int)fUsable.size();
  std::vector<double> t(nu), h(nu);
  for (size_t k = 0; k < fCl.size(); ++k) {
    const Cluster &q = fCl[k];
    double ch = 0, sum = 0;
    for (size_t j = 0; j < fClMem[k].size(); ++j) {
      const int i = fClMem[k][j], key = fPF.vref[i];
      if (fPF.charged[i] && key >= 0 && key < (int)fVtxIdx.size() && fVtxIdx[key] >= 0) ch += fPF.pt[i];
    }
    for (size_t j = 0; j < q.chPt.size(); ++j) sum += q.chPt[j];
    fIdCh = std::max(fIdCh, fabs(ch - sum));
    std::fill(t.begin(), t.end(), 0.); std::fill(h.begin(), h.end(), 0.);
    for (size_t j = 0; j < q.trkVtx.size(); ++j) t[q.trkVtx[j]] = q.trkPt[j];
    for (size_t j = 0; j < q.chVtx.size(); ++j)  h[q.chVtx[j]]  = q.chPt[j];
    std::vector<char> isCand(nu, 0);
    for (size_t j = 0; j < q.cand.size(); ++j) {
      const int u = q.cand[j]; isCand[u] = 1;
      // a candidate's T_ku below trkMinPt is not in trkVtx (only possible in candMode 1-2)
      if ((t[u] > 0 && q.candTrkPt[j] != t[u]) || q.candChPt[j] != h[u]) fIdCand += 1;
    }
    for (int u = 0; u < nu; ++u) {
      const bool bt = t[u] > 0, bh = h[u] >= fCfg.chMinPt && h[u] > 0;
      const bool rule = fCfg.candMode == 1 ? bh : (fCfg.candMode == 2 ? (bt || bh) : bt);
      if (rule != (isCand[u] != 0)) fIdCand += 1;
    }
  }
}

// The oracles once more, from the GEN side of the links (gfirst/gpf instead of
// pfirst/pgen), dense per cluster and per candidate: a second path to the same
// numbers, so that an indexing slip in Truth() cannot pass for a result.
void JVA::CheckTruth()
{
  const int nk = (int)fCl.size(), nf = (int)fFwd.size(), ng = (int)fGen.size(), ns = fNSlot;
  if (!nk) return;
  const genlink::Links &L = fLinks;
  std::vector<int> pos(fPF.size(), -1);
  for (int j = 0; j < nf; ++j) pos[fFwd[j]] = j;
  std::vector<double> cl((size_t)nk*ns, 0.), part((size_t)nf*ns, 0.);
  for (int g = 0; g < ng; ++g)
    for (int n = L.gfirst[g]; n < L.gfirst[g] + L.gcount[g]; ++n) {
      const int i = L.gpf[n];
      if (pos[i] < 0 || fClOf[i] < 0) continue;
      const double x = L.gfrac[n]*fGen.pt[g];
      part[(size_t)pos[i]*ns + fGen.inter[g]] += x;
      cl[(size_t)fClOf[i]*ns + fGen.inter[g]] += fWFwd[i]*x;
    }
  auto trueVtx = [&](int inter, double share) {
    if (inter < 0 || !(share > 0.5)) return -1;
    const int v = fNearVtx[inter];
    if (v < 0 || !(fabs(fZv[v] - fVz[inter]) < fDzOracle)) return -1;
    return fVtxIdx[v];
  };
  for (int k = 0; k < nk; ++k) {
    double sh; const int d = DominantRow(&cl[(size_t)k*ns], ns, sh);
    if (d != fClDom[k] || trueVtx(d, sh) != fClTrue[k] || trueVtx(d, sh) != fHost[kOjet][k]) fIdTruth += 1;
  }
  for (int j = 0; j < nf; ++j) {
    double sh; const int d = DominantRow(&part[(size_t)j*ns], ns, sh);
    if (trueVtx(d, sh) != fPartHost[fFwd[j]]) fIdTruth += 1;
  }
}

// Every method's host per cluster (lv trk jva ojet none) or per candidate
// (opart), and from them the list of vertex-blind candidates each vertex gets.
// With no usable PV (fPVIdx = -1) lv and trk's fallback give to nobody; none
// gives to nobody by definition.
void JVA::Methods()
{
  const int nk = (int)fCl.size(), nu = (int)fUsable.size(), np = (int)fPF.size();
  for (int m = 0; m < kNMeth; ++m) fHost[m].assign(nk, kNull);
  for (int k = 0; k < nk; ++k) {
    fHost[kLv][k] = fPVIdx;
    int t = fPVIdx; double best = -1;             // the most evidence (T_ku in candMode 0), the lowest index on a tie
    for (size_t j = 0; j < fCl[k].cand.size(); ++j) {
      const double e = jva::Evidence(fCl[k], j, fCfg);
      if (e > best) { best = e; t = fCl[k].cand[j]; }
    }
    fHost[kTrk][k] = t;
    fHost[kOjet][k] = fClTrue[k];
  }
  TStopwatch sw; sw.Start();
  fHost[kJva] = jva::Assign(fVtx, fCl, fCfg, &fStats);
  fAssignMs += 1000.*sw.RealTime();

  for (int m = 0; m < kNMeth; ++m) {
    fFwdTo[m].resize(nu);
    for (int u = 0; u < nu; ++u) fFwdTo[m][u].clear();
  }
  for (size_t j = 0; j < fFwd.size(); ++j) {
    const int i = fFwd[j], k = fClOf[i];
    for (int u = 0; u < nu; ++u) fFwdTo[kDup][u].push_back(i);
    if (k < 0) continue;                           // counted in fIdCover
    for (int m = kLv; m < kNMeth; ++m) {
      const int h = (m == kOpart) ? fPartHost[i] : fHost[m][k];
      if (h >= 0 && h < nu) fFwdTo[m][h].push_back(i);
    }
  }
  // the identities: dup gives everything to everybody, the others each
  // candidate to at most one vertex
  for (int u = 0; u < nu; ++u) fIdDup = std::max(fIdDup, fabs(double(fFwdTo[kDup][u].size()) - double(fFwd.size())));
  for (int u = 0; u < nu; ++u) fIdNone += (double)fFwdTo[kNobody][u].size();
  for (int k = 0; k < nk; ++k) if (fHost[kNobody][k] != kNull) fIdNone += 1;
  std::vector<int> cnt(np);
  for (int m = kLv; m < kNMeth; ++m) {
    std::fill(cnt.begin(), cnt.end(), 0);
    for (int u = 0; u < nu; ++u) for (size_t j = 0; j < fFwdTo[m][u].size(); ++j) ++cnt[fFwdTo[m][u][j]];
    for (int i = 0; i < np; ++i) fIdOnce = std::max(fIdOnce, double(cnt[i] - 1));
  }
}

// fwd/: every cluster once per crossing, and its outcome under each
// cluster-level method.  Single = share > 0.5 with a usable vertex within
// dzOracle of its interaction (fClTrue >= 0); anything else, unlinked
// included, is a combination.
void JVA::FillFwd()
{
  for (size_t k = 0; k < fCl.size(); ++k) {
    const Cluster &q = fCl[k];
    const int t = q.cand.empty() ? 1 : 0;
    hclpt[t]->Fill(q.pt);
    hshare->Fill(q.pt, std::min(fClShare[k], kShareMax));
    hcand->Fill((double)q.cand.size());
    const int tv = fClTrue[k];
    for (int c = 0; c < kNCat; ++c) {
      const int h = fHost[CatMeth(c)][k];
      int o;
      if (tv >= 0) o = (h == tv) ? kRight : (h < 0 ? kNulledSingle : kWrong);
      else         o = (h < 0) ? kNulledCombo : kKeptCombo;
      hcat[c][t]->Fill(q.pt, o + 0.5);
    }
    // what the evidence knows about a single-interaction cluster above 5 GeV:
    // is its true vertex among the candidates, the trk choice, among / leading
    // the u with T_ku >= trkMinPt, the u with C_ku > 0 (lowest index on a tie)
    if (tv < 0 || !(q.pt > 5)) continue;
    auto among = [tv](const std::vector<int> &v) { return std::find(v.begin(), v.end(), tv) != v.end(); };
    auto leads = [tv](const std::vector<int> &v, const std::vector<double> &e) {
      int a = -1; double b = -1; for (size_t j = 0; j < v.size(); ++j) if (e[j] > b) { b = e[j]; a = v[j]; } return a == tv; };
    fS5 += 1;
    if (among(q.cand)) fS5In += 1;
    if (!q.cand.empty() && fHost[kTrk][k] == tv) fS5Top += 1;
    if (among(q.trkVtx)) fS5InT += 1;
    if (leads(q.trkVtx, q.trkPt)) fS5TopT += 1;
    if (among(q.chVtx)) fS5InC += 1;
    if (leads(q.chVtx, q.chPt)) fS5TopC += 1;
  }
}

// The pure GenJets of one interaction (anti-kT R from ptMinCluster), on demand.
const std::vector<Jet> &JVA::GenJets(int inter)
{
  std::map<int, std::vector<Jet> >::iterator it = fGenJets.find(inter);
  if (it != fGenJets.end()) return it->second;
  std::vector<double> px, py, pz, E; std::vector< std::vector<int> > cons;
  for (int g = fIBase[inter]; g < fIEnd[inter]; ++g) {
    const double pt = fGen.pt[g], f = fGen.phi[g], m = fGen.mass[g];
    const double x = pt*cos(f), y = pt*sin(f), z = pt*sinh(fGen.eta[g]);
    px.push_back(x); py.push_back(y); pz.push_back(z); E.push_back(sqrt(x*x + y*y + z*z + m*m));
  }
  return fGenJets[inter] = ClusterJets(px, py, pz, E, fCfg.R, fPtMinCluster, cons);
}

////////////////////////////////////////////////////////////////////////////
// The seven events of usable vertex u: the vertex-resolved candidates at
// w_i(u) plus the vertex-blind ones the method gives u, at wFwd; with none
// given, the vertex-resolved event of VertexRecoil, jets and MET as they are.
// Jets with their correction, truth (pT^linkz, owner fraction, dominant
// share), the raw and type-1 MET, and every fill.
void JVA::ProcessVertex(int u)
{
  const int key = fUsable[u], owner = fGenOf[key];
  const bool owned = owner >= 0, isPV = (key == 0);
  const std::vector<float> &w = fWv[u];
  const double sv = fPMet[u].s;
  const genlink::Links &L = fLinks;
  const std::vector<int> &base = fBase[u];

  Ev ev; std::vector<Jet> fresh; std::vector< std::vector<int> > freshCons;
  SlotSum js; js.reset(fNSlot);
  const int cls[2] = {kAll, isPV ? kPV : kPU};
  const double tx = owned ? fTx[owner] : 0., ty = owned ? fTy[owner] : 0.;
  if (owned) for (int c = 0; c < 2; ++c) { fTrueN[cls[c]] += 2; fTrueS2[cls[c]] += tx*tx + ty*ty; }

  for (int m = 0; m < kNMeth; ++m) {
    const std::vector<int> &fw = fFwdTo[m][u];
    const bool cached = fw.empty();
    ev.build(fPF, base, w, fw, fWFwd, !cached);       // without kinematics when cached: idx, wt
    if (!cached) { fresh = ClusterJets(ev.px, ev.py, ev.pz, ev.E, fCfg.R, fPtMinCluster, freshCons); ++fNClusterings; }
    const double mx = cached ? fBaseMx[u] : ev.mx, my = cached ? fBaseMy[u] : ev.my;
    const std::vector<Jet> &jets = cached ? fBaseJets[u] : fresh;
    const std::vector< std::vector<int> > &cons = cached ? fBaseCons[u] : freshCons;
    const std::vector<int> &idx = ev.idx; const std::vector<float> &wt = ev.wt;
    std::vector<RJet> rj(jets.size());
    double t1x = mx, t1y = my;
    for (size_t k = 0; k < jets.size(); ++k) {
      RJet &r = rj[k];
      r.pt = jets[k].pt; r.eta = jets[k].eta; r.phi = jets[k].phi;
      r.ptcorr = r.pt*fJEC(r.pt, r.eta);
      double zx = 0, zy = 0; js.clear();
      const bool hard = r.pt > 5;
      for (size_t c = 0; c < cons[k].size(); ++c) {
        const int i = idx[cons[k][c]]; const double wi = wt[cons[k][c]];
        zx += wi*L.Gpx[i]; zy += wi*L.Gpy[i];
        if (hard)
          for (int n = L.pfirst[i]; n < L.pfirst[i] + L.pcount[i]; ++n)
            js.add(fGen.inter[L.pgen[n]], wi*L.pfrac[n]*fGen.pt[L.pgen[n]]);
      }
      r.linkz = sqrt(zx*zx + zy*zy);
      r.fown = (owned && js.total > 0) ? js.v[owner]/js.total : 0.;
      js.dominant(r.dshare);
      if (r.ptcorr > fT1Min) { t1x -= (r.ptcorr - r.pt)*cos(r.phi); t1y -= (r.ptcorr - r.pt)*sin(r.phi); }
    }

    // ---- MET ------------------------------------------------------------------
    const double amet = sqrt(mx*mx + my*my);
    for (int c = 0; c < 2; ++c) {
      const int cc = cls[c];
      hmetxy[m][cc]->Fill(mx); hmetxy[m][cc]->Fill(my); hmet[m][cc]->Fill(amet);
      fMetN[m][cc] += 2; fMetS2[m][cc] += mx*mx + my*my;
      if (owned) {
        hmetres[m][cc]->Fill(mx - tx); hmetres[m][cc]->Fill(my - ty);
        fResN[m][cc] += 2; fResS2[m][cc] += (mx - tx)*(mx - tx) + (my - ty)*(my - ty);
      }
    }
    hmetxy_npu[m]->Fill(fNPU, mx); hmetxy_npu[m]->Fill(fNPU, my);
    hmetxy_s[m]->Fill(sv, mx); hmetxy_s[m]->Fill(sv, my);
    if (owned) { hmetres_s[m]->Fill(sv, mx - tx); hmetres_s[m]->Fill(sv, my - ty); }
    if (isPV && m <= kLv) { fMetPV[m][0] = mx; fMetPV[m][1] = my; fHavePV[m] = true; fNJetPV[m] = (int)jets.size(); }

    // ---- the forward spike, owned vertices ---------------------------------------
    if (owned)
      for (size_t k = 0; k < rj.size(); ++k) {
        const RJet &r = rj[k]; const double ae = fabs(r.eta); const int rg = Reg(ae);
        for (int x = 0; x < 3; ++x) if (r.pt > kPtX[x]) hjeteta[m][x]->Fill(ae);
        if (rg >= 0) hjetpt[m][rg]->Fill(r.ptcorr);
        if (r.pt > 5) {
          pown[m]->Fill(ae, r.fown);
          if (rg >= 0) hdom[m][rg]->Fill(std::min(r.dshare, kShareMax));
        }
      }

    // ---- dijets ----------------------------------------------------------------
    FillDijet(m, 0, rj, t1x, t1y);
    if (isPV) FillDijet(m, 1, rj, t1x, t1y);
  }

  // the reference: the owner's pure GenJets, once per owned vertex
  if (owned) {
    const std::vector<Jet> &gj = GenJets(owner);
    for (size_t k = 0; k < gj.size(); ++k) {
      const double ae = fabs(gj[k].eta); const int rg = Reg(ae);
      for (int x = 0; x < 3; ++x) if (gj[k].pt > kPtX[x]) hjeteta[kNMeth][x]->Fill(ae);
      if (rg >= 0) hjetpt[kNMeth][rg]->Fill(gj[k].pt);
    }
  }
}

// One vertex event as a dijet: the two leading jets by pT^corr among those
// above ptMinJet within |eta| < 5.191, back to back.  Tag-probe with a barrel
// tag (both barrel: both ways at 0.5); the same-bin pairs for the JER, both
// signs at 0.5 (see JVA.h, THE DIJETS).
void JVA::FillDijet(int m, int c, const std::vector<RJet> &jets, double t1x, double t1y)
{
  TH1D *hs = hsel[m][c];
  hs->Fill(0.5);
  std::vector<const RJet*> s;
  for (size_t k = 0; k < jets.size(); ++k)
    if (jets[k].ptcorr > fPtMinJet && fabs(jets[k].eta) < kEtaJetMax) s.push_back(&jets[k]);
  if (s.size() < 2) return;
  std::stable_sort(s.begin(), s.end(), [](const RJet *a, const RJet *b){ return a->ptcorr > b->ptcorr; });
  hs->Fill(1.5);
  const RJet &j1 = *s[0], &j2 = *s[1];
  if (!(fabs(dphi(j1.phi, j2.phi)) > fDphiMin)) return;
  hs->Fill(2.5);
  const double ptavg = 0.5*(j1.ptcorr + j2.ptcorr);
  const double alpha = s.size() > 2 ? s[2]->ptcorr/ptavg : 0.;
  const bool a0 = alpha < kAlphaMax;

  // ---- tag and probe ------------------------------------------------------------
  const bool b1 = fabs(j1.eta) < kEtaTag, b2 = fabs(j2.eta) < kEtaTag;
  if (b1 || b2) {
    hs->Fill(3.5);
    if (a0) hs->Fill(4.5);
    const RJet *tg[2] = {b1 ? &j1 : &j2, &j1}, *pr[2] = {b1 ? &j2 : &j1, &j2};
    int np = 1; double wgt = 1.;
    if (b1 && b2) { tg[0] = &j1; pr[0] = &j2; tg[1] = &j2; pr[1] = &j1; np = 2; wgt = 0.5; }
    for (int p = 0; p < np; ++p) {
      const RJet &t = *tg[p], &q = *pr[p];
      const int ie = EtaBinDj(fabs(q.eta));
      if (ie < 0) continue;
      const double adb  = (q.ptcorr - t.ptcorr)/(q.ptcorr + t.ptcorr);
      const double ampf = (t1x*cos(t.phi) + t1y*sin(t.phi))/(2.*ptavg);
      const bool truth = q.linkz > 0 && t.linkz > 0;
      const double rt = truth ? (q.ptcorr/q.linkz)/(t.ptcorr/t.linkz) : 0.;
      for (int ia = 0; ia < 2; ++ia) {
        if (ia == 0 && !a0) continue;
        hadb[m][c][ie][ia]->Fill(ptavg, adb, wgt);
        hampf[m][c][ie][ia]->Fill(ptavg, ampf, wgt);
        if (truth) hrtrue[m][c][ie][ia]->Fill(ptavg, rt, wgt);
      }
    }
  }

  // ---- same-bin dijets: the bisector projections ----------------------------------
  const int ie1 = EtaBinDj(fabs(j1.eta)), ie2 = EtaBinDj(fabs(j2.eta));
  if (ie1 < 0 || ie1 != ie2 || !a0) return;
  const double u1x = cos(j1.phi), u1y = sin(j1.phi), u2x = cos(j2.phi), u2y = sin(j2.phi);
  double nx = u1x - u2x, ny = u1y - u2y;
  const double nn = sqrt(nx*nx + ny*ny);
  if (!(nn > 0)) return;
  hs->Fill(5.5);
  nx /= nn; ny /= nn;
  const double qx = -ny, qy = nx;                  // n_perp
  const double sx = j1.ptcorr*u1x + j2.ptcorr*u2x, sy = j1.ptcorr*u1y + j2.ptcorr*u2y;
  const double v[4] = {(sx*nx + sy*ny)/ptavg, (sx*qx + sy*qy)/ptavg, (t1x*nx + t1y*ny)/ptavg, (t1x*qx + t1y*qy)/ptavg};
  TH2F *h[4] = {hdbpar[m][c][ie1], hdbperp[m][c][ie1], hmpfpar[m][c][ie1], hmpfperp[m][c][ie1]};
  for (int k = 0; k < 4; ++k) { h[k]->Fill(ptavg, v[k], 0.5); h[k]->Fill(ptavg, -v[k], 0.5); }
  if (j1.linkz > 0) hrlink[m][c][ie1]->Fill(ptavg, j1.ptcorr/j1.linkz);
  if (j2.linkz > 0) hrlink[m][c][ie1]->Fill(ptavg, j2.ptcorr/j2.linkz);
}

void JVA::FillTuple()
{
  bRun = fReader->Run(); bLumi = fReader->Lumi(); bEvent = fReader->Event();
  bNPU = fNPU; bNVtx = (int)fUsable.size();
  bNComp = fStats.components; bNExact = fStats.exact; bObj = fStats.objective;
  bVtxKey.clear(); bVtxZ.clear(); bVtxOwner.clear(); bVtxMx.clear(); bVtxMy.clear(); bVtxS.clear();
  for (int t = 0; t < kNT; ++t) {
    bVtxHt[t].clear(); bVtxHtSq[t].clear(); bVtxMhtX[t].clear(); bVtxMhtY[t].clear(); bVtxGenMhtX[t].clear(); bVtxGenMhtY[t].clear(); }
  bVtxGenMetX.clear(); bVtxGenMetY.clear(); bVtxGenFwdX.clear(); bVtxGenFwdY.clear(); bVtxGenCenX.clear(); bVtxGenCenY.clear();
  for (size_t u = 0; u < fUsable.size(); ++u) {
    const int v = fUsable[u], own = fGenOf[v];
    bVtxKey.push_back(v); bVtxZ.push_back(fZv[v]); bVtxOwner.push_back(own);
    bVtxMx.push_back(fPMet[u].mx); bVtxMy.push_back(fPMet[u].my); bVtxS.push_back(fPMet[u].s);
    double gm[kNT][2] = {{0}};
    if (own >= 0) {
      const std::vector<Jet> &gj = GenJets(own);      // already made by ProcessVertex
      for (size_t k = 0; k < gj.size(); ++k) {
        if (!(fabs(gj[k].eta) < fCfg.etaFwd)) continue;
        for (int t = 0; t < kNT; ++t) if (gj[k].pt > kT[t]) { gm[t][0] += gj[k].pt*cos(gj[k].phi); gm[t][1] += gj[k].pt*sin(gj[k].phi); }
      }
    }
    for (int t = 0; t < kNT; ++t) {
      const size_t n = u*kNT + t;
      bVtxHt[t].push_back(fHtT[n]); bVtxHtSq[t].push_back(fHtSqT[n]); bVtxMhtX[t].push_back(fMhtXT[n]); bVtxMhtY[t].push_back(fMhtYT[n]);
      bVtxGenMhtX[t].push_back(gm[t][0]); bVtxGenMhtY[t].push_back(gm[t][1]);
    }
    bVtxGenMetX.push_back(own >= 0 ? fTx[own] : 0.); bVtxGenMetY.push_back(own >= 0 ? fTy[own] : 0.);
    bVtxGenFwdX.push_back(own >= 0 ? fFx[own] : 0.); bVtxGenFwdY.push_back(own >= 0 ? fFy[own] : 0.);
    bVtxGenCenX.push_back(own >= 0 ? fCx[own] : 0.); bVtxGenCenY.push_back(own >= 0 ? fCy[own] : 0.);
  }
  bClPt.clear(); bClEta.clear(); bClPhi.clear(); bClPx.clear(); bClPy.clear(); bClShare.clear(); bClNopt.clear();
  bClNcand.clear(); bClOff.clear(); bCandVtx.clear(); bCandTrkPt.clear(); bClDom.clear(); bClTrue.clear();
  bClNch.clear(); bClChOff.clear(); bChVtx.clear(); bChPt.clear(); bClTop3.clear(); bClTop3Share.clear();
  bClLv.clear(); bClTrk.clear(); bClJva.clear(); bClOjet.clear(); bClNone.clear();
  for (size_t k = 0; k < fCl.size(); ++k) {
    const Cluster &q = fCl[k];
    bClPt.push_back(q.pt); bClEta.push_back(q.eta); bClPhi.push_back(q.phi); bClPx.push_back(q.px); bClPy.push_back(q.py);
    bClNopt.push_back((int)q.cand.size());
    bClNcand.push_back((int)q.trkVtx.size()); bClOff.push_back((int)bCandVtx.size());
    for (size_t j = 0; j < q.trkVtx.size(); ++j) { bCandVtx.push_back(q.trkVtx[j]); bCandTrkPt.push_back(q.trkPt[j]); }
    bClNch.push_back((int)q.chVtx.size()); bClChOff.push_back((int)bChVtx.size());
    for (size_t j = 0; j < q.chVtx.size(); ++j) { bChVtx.push_back(q.chVtx[j]); bChPt.push_back(q.chPt[j]); }
    bClDom.push_back(fClDom[k]); bClShare.push_back(fClShare[k]); bClTrue.push_back(fClTrue[k]);
    for (int j = 0; j < 3; ++j) { bClTop3.push_back(fClTop3[3*k + j]); bClTop3Share.push_back(fClTop3Share[3*k + j]); }
    bClLv.push_back(fHost[kLv][k]); bClTrk.push_back(fHost[kTrk][k]); bClJva.push_back(fHost[kJva][k]); bClOjet.push_back(fHost[kOjet][k]);
    bClNone.push_back(fHost[kNobody][k]);
  }
  fTup->Fill();
}

////////////////////////////////////////////////////////////////////////////
void JVA::Loop()
{
  if (!fTree) { printf("JVA: no tree\n"); return; }
  fReader = new nanoaod::Reader(fTree);
  if (!fReader->Ok()) { printf("JVA: the tree lacks required branches, stopping\n"); return; }
  // A table that was asked for and does not load stops the job before any
  // file is written: the dijet part would silently be uncorrected.
  if (!fJECFile.empty() && !fJEC.Load(fJECFile.c_str())) {
    printf("JVA: cannot load the correction %s, stopping\n", fJECFile.c_str()); return; }
  fCfg.R         = Opt("R", 0.4);
  fCfg.etaFwd    = Opt("etaFwd", 2.5);
  fCfg.trkMinPt  = Opt("trkMinPt", 0.5);
  fCfg.sigma0    = Opt("sigma0", 1.0);
  fCfg.sigmaK    = Opt("sigmaK", 1.0);
  fCfg.tauTrk    = Opt("tauTrk", 4.0);
  fCfg.tauAll    = Opt("tauAll", 12.0);
  fCfg.maxCombos = Opt("maxCombos", 200000);
  fCfg.maxSweeps = int(Opt("maxSweeps", 50));
  fCfg.candMode  = int(Opt("candMode", 0));
  fCfg.chMinPt   = Opt("chMinPt", 0.5);
  fCfg.recoilMode = int(Opt("recoilMode", 0));
  fCfg.mhtMin    = Opt("mhtMin", 10.);
  fCfg.sigmaJ    = Opt("sigmaJ", 0.35);
  if (fCfg.candMode < 0 || fCfg.candMode > 2 || fCfg.recoilMode < 0 || fCfg.recoilMode > 1) {
    printf("JVA: candMode %d (0 tracks, 1 in-cluster charged, 2 either) or recoilMode %d (0 particle MET, 1 MHT) unknown, stopping\n",
           fCfg.candMode, fCfg.recoilMode); return; }
  fPtMinCluster  = Opt("ptMinCluster", 1.0);
  fPtMinJet      = Opt("ptMinJet", 3.0);
  fT1Min         = Opt("t1Min", 3.0);
  fDphiMin       = Opt("dphiMin", 2.7);
  fDzOracle      = Opt("dzOracle", 0.2);
  fStoreTuple    = Opt("storeTuple", 1) != 0;
  // PUPPI's own knobs, for the test of a vertex association extended into the
  // forward pixel coverage: past etaTracker PUPPI is vertex-blind, past
  // etaVtxAssoc a charged candidate is not tied to a vertex, and the neutral
  // floor of the region etaBound[0] < |eta| < etaBound[1] (2.5-3.0) is
  // pupFloor1Pt + pupFloor1Slope N_PV (puppi.h defaults 2.5, 2.5, 1.85, 0.070).
  // Setting pupEtaTracker = pupEtaVtxAssoc = etaFwd = 3.0 makes the tracks
  // with 2.5 < |eta| < 3.0 (their vertexRef, ignored by default) and the alpha
  // of the neutrals there per-vertex.
  fPupCfg.etaTracker     = Opt("pupEtaTracker", fPupCfg.etaTracker);
  fPupCfg.etaVtxAssoc    = Opt("pupEtaVtxAssoc", fPupCfg.etaVtxAssoc);
  fPupCfg.neutralPt[1]   = Opt("pupFloor1Pt", fPupCfg.neutralPt[1]);
  fPupCfg.neutralSlope[1]= Opt("pupFloor1Slope", fPupCfg.neutralSlope[1]);
  if (fabs(fCfg.etaFwd - fPupCfg.etaTracker) > 1e-9)
    printf("JVA: WARNING etaFwd %g differs from PUPPI's etaTracker %g: the \"vertex-blind\" candidates are then not all vertex-blind\n",
           fCfg.etaFwd, fPupCfg.etaTracker);
  Book();
  if (!fOut) return;
  if (fStoreTuple) BookTree();

  const Long64_t nall  = fReader->GetEntries();
  const Long64_t first = std::max((Long64_t)0, (Long64_t)Opt("firstEntry", 0));
  const Long64_t last  = std::min(nall, (Long64_t)Opt("lastEntry", nall));
  const Long64_t progressEvery = std::max((Long64_t)1, (Long64_t)Opt("progressEvery", 200));
  const double kv[kNKnob] = {fCfg.R, fPtMinCluster, fCfg.etaFwd, fCfg.trkMinPt, fCfg.sigma0, fCfg.sigmaK,
    fCfg.tauTrk, fCfg.tauAll, fCfg.maxCombos, (double)fCfg.maxSweeps, fPtMinJet, fT1Min, fDphiMin,
    Opt("dzVertexGen",0.05), fDzOracle, Opt("minTrkPerVertex",3), fStoreTuple ? 1. : 0.,
    fJEC.IsLoaded() ? 1. : 0., (double)first, (double)last, (double)progressEvery,
    strcmp(ClusterBackend(),"FastJet") == 0 ? 1. : 0.,
    (double)fCfg.candMode, fCfg.chMinPt, (double)fCfg.recoilMode, fCfg.mhtMin, fCfg.sigmaJ,
    fPupCfg.etaTracker, fPupCfg.etaVtxAssoc, fPupCfg.neutralPt[1], fPupCfg.neutralSlope[1]};
  for (int i = 0; i < kNKnob; ++i) hopts->SetBinContent(i+1, kv[i]);
  hcount->Fill(5.5);                                          // jobs

  printf("JVA: entries [%lld, %lld) of %lld, %d interaction slots, %s clustering, R %g from %g GeV, etaFwd %g, "
         "trkMinPt %g, sigma^2 = %g^2 + %g^2 S, tauTrk %g, tauAll %g, maxCombos %g, maxSweeps %d, correction %s, puppi %s\n",
         first, last, nall, fReader->NSlots(), ClusterBackend(), fCfg.R, fPtMinCluster, fCfg.etaFwd, fCfg.trkMinPt,
         fCfg.sigma0, fCfg.sigmaK, fCfg.tauTrk, fCfg.tauAll, fCfg.maxCombos, fCfg.maxSweeps,
         fJEC.IsLoaded() ? fJECFile.c_str() : "none (the dijet part is uncorrected)", puppi::kVersion);
  printf("  candidates (candMode %d): %s; recoil (recoilMode %d): %s\n", fCfg.candMode,
         fCfg.candMode == 1 ? Form("in-cluster charged C_ku >= %g GeV", fCfg.chMinPt)
         : fCfg.candMode == 2 ? Form("in-cone tracks T_ku >= %g GeV or in-cluster charged C_ku >= %g GeV", fCfg.trkMinPt, fCfg.chMinPt)
         : Form("in-cone tracks T_ku >= %g GeV", fCfg.trkMinPt), fCfg.recoilMode,
         fCfg.recoilMode == 1 ? Form("MHT of the vertex-resolved jets above %g GeV raw, sigma^2 = %g^2 + %g^2 sum pT_j^2 + %g^2 (S - HT)",
                                     fCfg.mhtMin, fCfg.sigma0, fCfg.sigmaJ, fCfg.sigmaK)
                              : "particle MET of the vertex-resolved part");
  printf("  dijets: pT^corr > %g, |eta| < %g, dphi > %g, tag |eta| < %g, alpha < %g; type-1 from %g GeV; oracle dz %g cm\n",
         fPtMinJet, kEtaJetMax, fDphiMin, kEtaTag, kAlphaMax, fT1Min, fDzOracle);

  // An empty range writes a valid empty file (see GenSeed::Loop).
  if (last <= first || !fReader->SetEntriesRange(first, last)) {
    printf("JVA: empty range [%lld, %lld), nothing to do\n", first, last); Write(); return; }
  TStopwatch sw; sw.Start();
  Long64_t nproc = 0;
  while (fReader->Next()) {
    const Long64_t jentry = fReader->Entry();
    if (jentry >= last) break;
    hcount->Fill(4.5);
    const int half = int(jentry % 2);                         // on the GLOBAL entry
    hcount->Fill(half + 0.5);
    hcount->AddBinContent(3 + half, fReader->NPU());
    ++nproc;

    ReadGen();
    fReader->FillPF(fPF);
    puppi::Prepare(fPF, fPupCfg);
    ReconstructVertices();
    genlink::Link(fPF, fGen, fReader->PVz(), fLinkCfg, fLinks);
    VertexWeights();
    VertexRecoil();
    if (fUsable.empty()) { fCl.clear(); fClMem.clear(); ++fNNoPV; }
    else fCl = jva::BuildClusters(fPF, fWFwd, fVtxIdx, fCfg, &fClMem);
    if (fPVIdx < 0 && !fUsable.empty()) ++fNNoPV;
    Truth();
    CheckEvidence();
    Methods();
    CheckTruth();
    FillFwd();

    fHavePV[0] = fHavePV[1] = false;
    int nown = 0;
    for (size_t u = 0; u < fUsable.size(); ++u) { ProcessVertex((int)u); if (fGenOf[fUsable[u]] >= 0) ++nown; }
    hcount->AddBinContent(7, (double)fUsable.size());
    hcount->AddBinContent(8, (double)nown);
    if (fHavePV[0] && fHavePV[1])
      fIdPV = std::max(fIdPV, std::max(std::max(fabs(fMetPV[0][0] - fMetPV[1][0]), fabs(fMetPV[0][1] - fMetPV[1][1])),
                                       fabs(double(fNJetPV[0] - fNJetPV[1]))));
    fNVtx += fUsable.size(); fNOwned += nown; fNCl += fCl.size(); fNFwdCand += fFwd.size();
    for (size_t k = 0; k < fCl.size(); ++k) { if (!fCl[k].cand.empty()) ++fNClTrk; fNCand += fCl[k].cand.size(); }
    fNComp += fStats.components; fNExact += fStats.exact; fNSweeps += fStats.sweeps; fLargest += fStats.largest;
    fObjSum += fStats.objective;
    if (fStoreTuple) FillTuple();

    if (nproc % progressEvery == 0) {
      const double t = sw.RealTime(); sw.Continue();
      ProcInfo_t pi; gSystem->GetProcInfo(&pi);
      printf("  %lld / %lld  (%.1f ms/crossing, Assign %.2f ms, %.1f forward clusters, %.0f s to go, entry %lld, %.0f MB resident)\n",
             jentry - first + 1, last - first, 1000.*t/nproc, fAssignMs/nproc, double(fNCl)/nproc,
             t*(last - jentry - 1)/nproc, jentry, pi.fMemResident/1024.);
      fflush(stdout);
    }
  }
  const double t = sw.RealTime();
  fNEvents = nproc;
  PrintSummary(t);
  CheckIdentities();
  Write();
}

// The run in numbers: the association, the MET per method and class, the
// cluster outcomes, the dijet selection.
void JVA::PrintSummary(double t)
{
  const double n = fNEvents > 0 ? double(fNEvents) : 1.;
  printf("JVA: done, %lld crossings in %.1f s: %.1f ms/crossing (Assign %.2f ms); %.1f usable vertices (%.1f owned) per crossing%s\n",
         fNEvents, t, 1000.*t/n, fAssignMs/n, fNVtx/n, fNOwned/n,
         fNNoPV ? Form(", %lld crossings without a usable PV", fNNoPV) : "");
  printf("  forward: %.1f vertex-blind candidates with wFwd > 0 and %.2f clusters per crossing, %.3f of the clusters with candidates "
         "(candMode %d), %.2f candidates per cluster;\n"
         "  jva: %.2f components per crossing, %.3f of them solved exactly, largest %.1f clusters, %.2f descent sweeps + pair passes, objective %.1f per crossing\n",
         fNFwdCand/n, fNCl/n, fNCl ? double(fNClTrk)/fNCl : 0., fCfg.candMode, fNCl ? double(fNCand)/fNCl : 0.,
         fNComp/n, fNComp ? double(fNExact)/fNComp : 0., fLargest/n, fNSweeps/n, fObjSum/n);
  printf("  recoil (recoilMode %d): rms %.2f GeV per component, <sigma_v^2> %.1f GeV^2 over the usable vertices; %.1f clusterings per crossing\n",
         fCfg.recoilMode, fNVtx ? sqrt(fRecoilS2/(2.*fNVtx)) : 0., fNVtx ? fSig2Sum/fNVtx : 0., fNClusterings/n);
  if (fS5 > 0)
    printf("  single-interaction clusters above 5 GeV (%.0f): true vertex among the candidates %.3f, the trk choice %.3f;\n"
           "    among / leading the u with T_ku >= %g GeV %.3f / %.3f, with C_ku > 0 %.3f / %.3f\n",
           fS5, fS5In/fS5, fS5Top/fS5, fCfg.trkMinPt, fS5InT/fS5, fS5TopT/fS5, fS5InC/fS5, fS5TopC/fS5);
  printf("  MET rms per component [GeV] (all = every usable vertex, pv = key 0, pu = the others); res = MET - MET_true(owner), owned vertices\n");
  printf("  %-6s %8s %8s %8s   %8s %8s %8s\n", "", "all", "pv", "pu", "res all", "res pv", "res pu");
  for (int m = 0; m < kNMeth; ++m) {
    printf("  %-6s", MethTag(m));
    for (int c = 0; c < kNClass; ++c) printf(" %8.2f", fMetN[m][c] > 0 ? sqrt(fMetS2[m][c]/fMetN[m][c]) : 0.);
    printf("  ");
    for (int c = 0; c < kNClass; ++c) printf(" %8.2f", fResN[m][c] > 0 ? sqrt(fResS2[m][c]/fResN[m][c]) : 0.);
    printf("\n");
  }
  printf("  %-6s %8s %8s %8s  ", "true", "", "", "");
  for (int c = 0; c < kNClass; ++c) printf(" %8.2f", fTrueN[c] > 0 ? sqrt(fTrueS2[c]/fTrueN[c]) : 0.);
  printf("   (rms of MET_true(owner) itself)\n");
  // outcomes, clusters above 5 GeV
  const int b5 = hclpt[0]->GetXaxis()->FindFixBin(5.0001);
  printf("  forward clusters above 5 GeV: outcome fractions  (right wrong nulled_single nulled_combo kept_combo)\n");
  for (int t = 0; t < 2; ++t) {
    const double tot = hclpt[t]->Integral(b5, kNPtExt + 1);
    printf("  %-5s tracks, %6.0f clusters:\n", t == 0 ? "with" : "no", tot);
    for (int c = 0; c < kNCat; ++c) {
      printf("    %-5s", MethTag(CatMeth(c)));
      for (int o = 0; o < kNOut; ++o)
        printf(" %6.3f", tot > 0 ? hcat[c][t]->Integral(b5, kNPtExt + 1, o + 1, o + 1)/tot : 0.);
      printf("\n");
    }
  }
  printf("  dijet selection, all vertices (vertices, >=2 jets, dphi, tag, alpha<0.3, samebin):\n");
  for (int m = 0; m < kNMeth; ++m) {
    printf("    %-5s", MethTag(m));
    for (int s = 1; s <= kNSel; ++s) printf(" %9.0f", hsel[m][0]->GetBinContent(s));
    printf("\n");
  }
}

// The identities, as the worst violation, which must be zero (see JVA.h).
void JVA::CheckIdentities()
{
  double wcat = 0, woj = 0, wnone = 0;
  for (int c = 0; c < kNCat; ++c)
    for (int t = 0; t < 2; ++t)
      for (int i = 0; i <= kNPtExt + 1; ++i) {
        double s = 0;
        for (int o = 1; o <= kNOut; ++o) s += hcat[c][t]->GetBinContent(i, o);
        wcat = std::max(wcat, fabs(s - hclpt[t]->GetBinContent(i)));
        if (CatMeth(c) == kOjet)
          woj = std::max(woj, hcat[c][t]->GetBinContent(i, kWrong + 1) + hcat[c][t]->GetBinContent(i, kNulledSingle + 1)
                              + hcat[c][t]->GetBinContent(i, kKeptCombo + 1));
        if (CatMeth(c) == kNobody)
          wnone = std::max(wnone, hcat[c][t]->GetBinContent(i, kRight + 1) + hcat[c][t]->GetBinContent(i, kWrong + 1)
                                  + hcat[c][t]->GetBinContent(i, kKeptCombo + 1));
      }
  auto line = [](const char *what, double v) { printf("    %-72s %9.3g  %s\n", what, v, v > 1e-9 ? "VIOLATED" : "ok"); };
  printf("  identities (worst violation, must be 0):\n");
  line("dup: every vertex gets every vertex-blind candidate", fIdDup);
  line("all but dup: each vertex-blind candidate to at most one vertex", fIdOnce);
  line("every vertex-blind candidate with wFwd > 0 in exactly one cluster", fIdCover);
  line("wFwd the same at every vertex (max |w_i(v) - wFwd_i|)", fIdWfwd);
  line("dup == lv at the PV (max |dMET| [GeV], |d njets|)", fIdPV);
  line("ojet and opart == truth recomputed from the gen side (mismatches)", fIdTruth);
  line("ojet never wrong, nulled_single or kept_combo (clusters)", woj);
  line("none: no vertex-blind candidate at any vertex, every cluster null", fIdNone + wnone);
  line("sum over outcomes of hcat = hclpt, bin by bin", wcat);
  line("sum_u C_ku = recounted charged constituent pT (max |d| [GeV])", fIdCh);
  line("candidates = the candMode rule on the evidence (mismatches)", fIdCand);
  line("MET_true = -(gencen + genfwd), every interaction (max |d| [GeV])", fIdGen);
}

void JVA::Write()
{
  const double pa = hcount->GetBinContent(3), pb = hcount->GetBinContent(4);
  printf("  events %.0f + %.0f, interactions %.0f + %.0f, entries %.0f, jobs %.0f, vertices %.0f (owned %.0f), tuple %lld entries\n",
         hcount->GetBinContent(1), hcount->GetBinContent(2), pa, pb, hcount->GetBinContent(5), hcount->GetBinContent(6),
         hcount->GetBinContent(7), hcount->GetBinContent(8), fTup ? fTup->GetEntries() : 0);
  fOut->Write();
  fOut->Close();
  printf("  wrote %s\n", fOutName.c_str());
}
