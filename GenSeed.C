// genseed - Copyright 2026 Mikko Voutilainen - Apache License 2.0
// GenSeed.C - see GenSeed.h for what this produces and why.
//
//   root -l -b -q 'runGenSeed.C("files.txt","_p1")'                       // pass 1
//   root -l -b -q 'drawGenSeed.C+("rootfiles/GenSeed_p1.root","_v1")'      // -> text/jec_v1.txt
//   root -l -b -q 'runGenSeed.C("files.txt","_v1",0,1,0,"text/jec_v1.txt")' // pass 2
//
// The per-event order is the one the dependencies dictate: the generated
// list and the pure jets (vertex independent), the PF event and puppi::Prepare
// (vertex independent), the vertices and their ownership, the linker (once,
// vertex independent), then every usable vertex in turn - PUPPI weights, the
// clustering with the own unlinked ghosts, the seeds and the reco-side
// topology of its jets, and, for the interactions whose nearest vertex it is,
// the gen side of their pure jets and their RecoJetSeeds - and at the end the
// dedup across hypotheses and every fill, once per kept jet, plus the tuple.
// Nothing is filled inside ProcessVertex except hnghost: every histogram sees
// the event once, after the dedup, from the stored RJet and PJet records, so
// the identities between the histograms hold by construction.
#include "GenSeed.h"

#include <TStopwatch.h>
#include <TSystem.h>
#include <TDirectory.h>
#include <TNamed.h>
#include <algorithm>
#include <cstdio>
#include <cstring>

using namespace genseed;

namespace {
  inline double dphi(double a, double b) {
    double d = a - b;
    if (d > M_PI) d -= 2*M_PI; else if (d < -M_PI) d += 2*M_PI;
    return d;
  }
  const char *kHalf[2] = {"a","b"};
  inline int Reg(double aeta) { return aeta < 2.5 ? 0 : (aeta < 3.0 ? 1 : 2); }
  const char *kRegTag[3] = {"trk","tran","hf"};
  // the knobs written to hist/hopts, in bin order
  const char *kKnob[] = {"R","ptMinCluster","ptMinGen","ptStoreReco","ptStoreGen","etaMaxPart",
    "shareFrac","dzVertexGen","dzRecoSeed","minTrkPerVertex","dupDR","dupPtFrac","kernelScale",
    "useDz","useCharge","parallaxRef","bendSign","jecLoaded","maxVertices","firstEntry","lastEntry",
    "progressEvery","useFastJet"};
  const int kNKnob = sizeof(kKnob)/sizeof(kKnob[0]);
  struct Acc { double sx, sy, sz, se, zx, zy, own, far, unl; };   // link 4-vector, linkz transverse, decomposition
  template <class T> void Vec(TTree *t, const char *n, std::vector<T> &v) { t->Branch(n, &v); }
}

double GenSeed::DeltaR2(double e1, double p1, double e2, double p2) {
  const double de = e1-e2, df = dphi(p1,p2); return de*de + df*df;
}

GenSeed::GenSeed(TTree *tree, const char *outname)
  : fOutName(outname), fJECFile(""), fTree(tree), fOut(0), fTup(0), fReader(0),
    fR(0.4), fShare(1./3), fPtMinCluster(1.), fPtMinGen(1.), fPtStoreReco(3.), fPtStoreGen(3.),
    fDzRecoSeed(0.2), fNPU(0), fNSlot(0), fLinkMs(0), fNEvents(0) {}
GenSeed::~GenSeed() { delete fReader; }

////////////////////////////////////////////////////////////////////////////
void GenSeed::Book()
{
  fOut = new TFile(fOutName.c_str(), "RECREATE");
  if (fOut->IsZombie()) { printf("GenSeed: cannot open %s\n", fOutName.c_str()); fOut = 0; return; }
  const int nb = kNPtExt; const double *xb = kPtBinsExt;

  TDirectory *dh = fOut->mkdir("hist"); dh->cd();
  hcount = new TH1D("hcount","counters",6,0,6);
  const char *cl[6] = {"events_a","events_b","npu_a","npu_b","entries","jobs"};
  for (int i = 0; i < 6; ++i) hcount->GetXaxis()->SetBinLabel(i+1, cl[i]);
  hopts = new TH1D("hopts","knobs, value x jobs",kNKnob,0,kNKnob);
  for (int i = 0; i < kNKnob; ++i) hopts->GetXaxis()->SetBinLabel(i+1, kKnob[i]);
  hnghost  = new TH1D("hnghost",";own unlinked ghosts per vertex clustering;clusterings",200,0,2000);
  hlinkms  = new TH1D("hlinkms",";Link() time per event [ms];events",200,0,100);
  hncopies = new TH1D("hncopies",";vertex copies per kept reco jet;jets",60,0.5,60.5);
  // What the file was made with that is not a number, as text next to hopts.
  // TNamed has no Merge, so hadd keeps every job's copy as a cycle of its own
  // (Get returns the last); harmless, the jobs of one merge share both (the
  // table is checked to have loaded, see Loop).
  TNamed("puppiVersion", puppi::kVersion).Write();
  TNamed("jecFile", fJEC.IsLoaded() ? fJECFile.c_str() : "none").Write();

  TDirectory *ds = fOut->mkdir("spec"); ds->cd();
  for (int iy = 0; iy < kNY; ++iy) for (int h = 0; h < 2; ++h) {
    const std::string yt = YTag(iy), yl = YLabel(iy);
    hspec_gen[iy][h]  = new TH1D(Form("hgen_%s_%s",yt.c_str(),kHalf[h]),  Form("%s;p_{T}^{pure} [GeV];pure gen jets",yl.c_str()),nb,xb);
    hspec_reco[iy][h] = new TH1D(Form("hreco_%s_%s",yt.c_str(),kHalf[h]), Form("%s;p_{T}^{corr} [GeV];kept reco jets",yl.c_str()),nb,xb);
    hspec_seed[iy][h] = new TH1D(Form("hseed_%s_%s",yt.c_str(),kHalf[h]), Form("%s;p_{T}^{link} [GeV];kept reco jets",yl.c_str()),nb,xb);
    hspec_rs[iy][h]   = new TH1D(Form("hrseed_%s_%s",yt.c_str(),kHalf[h]),Form("%s;p_{T}^{rs} [GeV];RecoJetSeeds",yl.c_str()),nb,xb);
  }

  TDirectory *de = fOut->mkdir("eff"); de->cd();
  for (int ie = 0; ie < jec::kNEta; ++ie) {
    const std::string et = jec::EtaTag(ie), el = jec::EtaLabel(ie);
    hgall[ie]      = new TH1D(Form("hgall_%s",et.c_str()),     Form("%s;p_{T}^{pure} [GeV];pure gen jets with a vertex",el.c_str()),nb,xb);
    hgmat_dr[ie]   = new TH1D(Form("hgmat_dr_%s",et.c_str()),  Form("%s, #DeltaR < 0.2 partner;p_{T}^{pure} [GeV];pure gen jets",el.c_str()),nb,xb);
    hgmat_dom[ie]  = new TH1D(Form("hgmat_dom_%s",et.c_str()), Form("%s, mutual dominance;p_{T}^{pure} [GeV];pure gen jets",el.c_str()),nb,xb);
    hgmat_rs[ie]   = new TH1D(Form("hgmat_rs_%s",et.c_str()),  Form("%s, RecoJetSeed p_{T} > 0;p_{T}^{pure} [GeV];pure gen jets",el.c_str()),nb,xb);
    hgmat_link[ie] = new TH1D(Form("hgmat_link_%s",et.c_str()),Form("%s, link topology not lost;p_{T}^{pure} [GeV];pure gen jets",el.c_str()),nb,xb);
  }

  TDirectory *df = fOut->mkdir("fake"); df->cd();
  for (int ie = 0; ie < jec::kNEta; ++ie) {
    const std::string et = jec::EtaTag(ie), el = jec::EtaLabel(ie);
    hrall[ie]      = new TH1D(Form("hrall_%s",et.c_str()),     Form("%s;p_{T}^{corr} [GeV];kept reco jets",el.c_str()),nb,xb);
    hrnone[ie]     = new TH1D(Form("hrnone_%s",et.c_str()),    Form("%s, no pure feeder;p_{T}^{corr} [GeV];kept reco jets",el.c_str()),nb,xb);
    hrunpaired[ie] = new TH1D(Form("hrunpaired_%s",et.c_str()),Form("%s, no mutual partner;p_{T}^{corr} [GeV];kept reco jets",el.c_str()),nb,xb);
  }

  TDirectory *dr = fOut->mkdir("resp"); dr->cd();
  const char *pu[3] = {"N_{PU} < 35","35 #leq N_{PU} #leq 55","N_{PU} > 55"};
  for (int ie = 0; ie < jec::kNEta; ++ie) {
    const std::string et = jec::EtaTag(ie), el = jec::EtaLabel(ie);
    hresp_link[ie]  = new TH2D(Form("hresp_link_%s",et.c_str()), Form("%s;p_{T}^{link} [GeV];p_{T}^{raw} / p_{T}^{link}",el.c_str()),nb,xb,kNR,0,kRMax);
    hresp_linkz[ie] = new TH2D(Form("hresp_linkz_%s",et.c_str()),Form("%s;p_{T}^{linkz} [GeV];p_{T}^{raw} / p_{T}^{linkz}",el.c_str()),nb,xb,kNR,0,kRMax);
    hresp_dr[ie]    = new TH2D(Form("hresp_dr_%s",et.c_str()),   Form("%s, #DeltaR < 0.2 pairs;p_{T}^{pure} [GeV];p_{T}^{raw} / p_{T}^{pure}",el.c_str()),nb,xb,kNR,0,kRMax);
    hresp_rs[ie]    = new TH2D(Form("hresp_rs_%s",et.c_str()),   Form("%s;p_{T}^{pure} [GeV];p_{T}^{rs} / p_{T}^{pure}",el.c_str()),nb,xb,kNR,0,kRMax);
    hresp_rs0[ie]   = new TH2D(Form("hresp_rs0_%s",et.c_str()),  Form("%s, no PUPPI weights;p_{T}^{pure} [GeV];p_{T}^{rs0} / p_{T}^{pure}",el.c_str()),nb,xb,kNR,0,kRMax);
    hresp_reco_rs[ie] = new TH2D(Form("hresp_reco_rs_%s",et.c_str()),Form("%s, mutual partner;p_{T}^{pure} [GeV];p_{T}^{raw} / p_{T}^{rs}",el.c_str()),nb,xb,kNR,0,kRMax);
    for (int t = 0; t < kNTopo; ++t)
      hrtopo[t][ie] = new TH2D(Form("hrtopo_link_%s_%s",RTopoTag(t),et.c_str()),
        Form("%s, %s;p_{T}^{link} [GeV];p_{T}^{raw} / p_{T}^{link}",TopoLabel(t),el.c_str()),nb,xb,kNR,0,kRMax);
    for (int s = 0; s < 3; ++s)
      hresppu[s][ie] = new TH2D(Form("hresppu_link_pu%d_%s",s,et.c_str()),
        Form("%s, %s;p_{T}^{link} [GeV];p_{T}^{raw} / p_{T}^{link}",pu[s],el.c_str()),nb,xb,kNR,0,kRMax);
    hfar[ie]  = new TProfile(Form("hfar_%s",et.c_str()), Form("%s;p_{T}^{link} [GeV];far fraction of the seed",el.c_str()),nb,xb);
    hunl[ie]  = new TProfile(Form("hunl_%s",et.c_str()), Form("%s;p_{T}^{link} [GeV];unlinked fraction of the seed",el.c_str()),nb,xb);
    horph[ie] = new TProfile(Form("horph_%s",et.c_str()),Form("%s;p_{T}^{link} [GeV];orphan fraction of the reco p_{T}",el.c_str()),nb,xb);
  }

  TDirectory *dt = fOut->mkdir("topo"); dt->cd();
  for (int iy = 0; iy < kNY; ++iy) {
    const std::string yt = YTag(iy), yl = YLabel(iy);
    htall[iy] = new TH1D(Form("hall_%s",yt.c_str()),Form("%s;p_{T}^{pure} [GeV];pure gen jets with a vertex",yl.c_str()),nb,xb);
    for (int t = 0; t < kNTopo; ++t) {
      htopo[t][iy]   = new TH1D(Form("h%s_%s",TopoTag(t),yt.c_str()),  Form("%s, %s;p_{T}^{pure} [GeV];pure gen jets",TopoLabel(t),yl.c_str()),nb,xb);
      hrtopo1[t][iy] = new TH1D(Form("hr%s_%s",RTopoTag(t),yt.c_str()),Form("%s, %s;p_{T}^{corr} [GeV];kept reco jets",TopoLabel(t),yl.c_str()),nb,xb);
    }
  }

  TDirectory *db = fOut->mkdir("bridge"); db->cd();
  for (int iy = 0; iy < kNY; ++iy) for (int h = 0; h < 2; ++h) {
    const std::string yt = YTag(iy), yl = YLabel(iy);
    hbridge[iy][h] = new TH2D(Form("hbridge_%s_%s",yt.c_str(),kHalf[h]),Form("%s, mutual pairs;p_{T}^{link} [GeV];p_{T}^{pure} [GeV]",yl.c_str()),nb,xb,nb,xb);
    hbfake[iy][h]  = new TH1D(Form("hbfake_%s_%s",yt.c_str(),kHalf[h]), Form("%s, unpaired seeds;p_{T}^{link} [GeV];kept reco jets",yl.c_str()),nb,xb);
    hbmiss[iy][h]  = new TH1D(Form("hbmiss_%s_%s",yt.c_str(),kHalf[h]), Form("%s, unpaired pure jets;p_{T}^{pure} [GeV];pure gen jets",yl.c_str()),nb,xb);
  }

  TDirectory *dm = fOut->mkdir("map"); dm->cd();
  hjet_reco = new TH2D("hjet_reco","kept reco jets, p_{T}^{corr} > 5 GeV;#eta;#phi",100,-5,5,72,-M_PI,M_PI);
  hjet_gen  = new TH2D("hjet_gen","pure gen jets, p_{T} > 5 GeV;#eta;#phi",100,-5,5,72,-M_PI,M_PI);
  hjet_rs   = new TH2D("hjet_rs","RecoJetSeeds, p_{T} > 5 GeV;#eta;#phi",100,-5,5,72,-M_PI,M_PI);
  presp_link = new TProfile2D("presp_link","p_{T}^{raw} / p_{T}^{link}, p_{T}^{link} > 5 GeV;#eta;#phi",100,-5,5,72,-M_PI,M_PI);
  presp_rs   = new TProfile2D("presp_rs","p_{T}^{rs} / p_{T}^{pure}, p_{T}^{pure} > 5 GeV;#eta;#phi",100,-5,5,72,-M_PI,M_PI);
  porph  = new TProfile2D("porph","orphan fraction of the reco p_{T}, kept jets;#eta;#phi",100,-5,5,72,-M_PI,M_PI);
  pflink = new TProfile2D("pflink","linked fraction of the generated p_{T};#eta;#phi",100,-5,5,72,-M_PI,M_PI);
  pfar   = new TProfile2D("pfar","far fraction of the seed, kept jets;#eta;#phi",100,-5,5,72,-M_PI,M_PI);

  TDirectory *dl = fOut->mkdir("link"); dl->cd();
  for (int c = 0; c < genlink::kNDiagClass; ++c) {
    const char *ct = genlink::DiagTag(c);
    for (int r = 0; r < 3; ++r)
      hlinkfrac[c][r] = new TProfile(Form("hlinkfrac_%s_%s",ct,kRegTag[r]),Form("%s, %s;p_{T}^{gen} [GeV];linked p_{T} fraction",ct,kRegTag[r]),60,0,30);
  }
  for (int r = 0; r < 3; ++r)
    hlinkdz[r] = new TH1D(Form("hlinkdz_%s",kRegTag[r]),Form("%s;z_{track} - z_{gen} [cm];pass-1 links",kRegTag[r]),200,-0.5,0.5);
  hparallax = new TH2D("hparallax","single-source #gamma links, |#eta| < 1.4;(z_{g} - PV_{z}) / (r_{ECAL} cosh#eta_{g});#eta_{PF} - #eta_{gen}",100,-0.1,0.1,150,-0.15,0.15);
  hbend = new TH2D("hbend","untracked charged #rightarrow calo links, barrel;1 / p_{T}^{gen} [GeV^{-1}];q (#phi_{PF} - #phi_{gen})",100,0,2,150,-1.5,1.5);
  for (int c = 0; c < genlink::kNPFClass; ++c) {
    const char *ct = genlink::PFClassTag(c);
    hpfratio[c] = new TH1D(Form("hpfratio_%s",ct),Form("PF %s, single gen partner;p_{T}^{PF} / G_{T};candidates",ct),300,0,3);
    for (int r = 0; r < 3; ++r)
      horphan[c][r] = new TProfile(Form("horphan_%s_%s",ct,kRegTag[r]),Form("PF %s, %s;N_{PU};orphan p_{T} fraction",ct,kRegTag[r]),24,0,120);
  }
  fOut->cd();
}

// The tuple: one entry per event, indices into the stored arrays of the
// same event, -1 when the partner is below its storage threshold.
void GenSeed::BookTree()
{
  fOut->cd();
  fTup = new TTree("genseed","genseed: pure gen jets, RecoJetSeeds, kept reco jets with GenJetSeeds");
  fTup->Branch("run",&bRun,"run/i"); fTup->Branch("lumi",&bLumi,"lumi/i");
  fTup->Branch("event",&bEvent,"event/l"); fTup->Branch("nPU",&bNPU,"nPU/I");
  fTup->Branch("nVtxUsable",&bNVtx,"nVtxUsable/I"); fTup->Branch("pvz",&bPVz,"pvz/F");
  fTup->Branch("lumiWeight",&bLumiW,"lumiWeight/F");
  Vec(fTup,"vtx_z",bVtxZ); Vec(fTup,"vtx_ntrk",bVtxNtrk); Vec(fTup,"vtx_gen",bVtxGen);
  Vec(fTup,"gj_pt",bGjPt); Vec(fTup,"gj_eta",bGjEta); Vec(fTup,"gj_phi",bGjPhi); Vec(fTup,"gj_mass",bGjMass);
  Vec(fTup,"gj_inter",bGjInter); Vec(fTup,"gj_ncon",bGjNcon); Vec(fTup,"gj_vtx",bGjVtx);
  Vec(fTup,"gj_dr02reco",bGjDr02); Vec(fTup,"gj_topo",bGjTopo); Vec(fTup,"gj_kmain",bGjKmain); Vec(fTup,"gj_share",bGjShare);
  Vec(fTup,"gj_rs_pt",bRsPt); Vec(fTup,"gj_rs0_pt",bRs0Pt); Vec(fTup,"gj_rs_eta",bRsEta); Vec(fTup,"gj_rs_phi",bRsPhi); Vec(fTup,"gj_rs_mass",bRsMass);
  Vec(fTup,"gj_rs_npf",bRsNpf); Vec(fTup,"gj_rs_flead",bRsFlead); Vec(fTup,"gj_rs_flinked",bRsFlinked); Vec(fTup,"gj_rs_fpuppi",bRsFpuppi);
  Vec(fTup,"rj_pt",bRjPt); Vec(fTup,"rj_ptcorr",bRjPtcorr); Vec(fTup,"rj_eta",bRjEta); Vec(fTup,"rj_phi",bRjPhi); Vec(fTup,"rj_mass",bRjMass);
  Vec(fTup,"rj_vtx",bRjVtx); Vec(fTup,"rj_gen",bRjGen);
  Vec(fTup,"rj_seed_link",bRjLink); Vec(fTup,"rj_seed_linkz",bRjLinkz);
  Vec(fTup,"rj_fown",bRjFown); Vec(fTup,"rj_ffar",bRjFfar); Vec(fTup,"rj_funl",bRjFunl); Vec(fTup,"rj_forph",bRjForph);
  Vec(fTup,"rj_topo",bRjTopo); Vec(fTup,"rj_pstar",bRjPstar); Vec(fTup,"rj_pshare",bRjPshare);
  Vec(fTup,"rj_bpair",bRjBpair); Vec(fTup,"rj_ncopies",bRjNcopies);
}

////////////////////////////////////////////////////////////////////////////
// The generated list of every interaction (nanoreader keeps each one
// contiguous), the four-vectors, and the pure jets: each interaction
// clustered on its own, from ptMinGen, every particle remembering its jet.
void GenSeed::ReadGen()
{
  fNPU = fReader->NPU();
  fNSlot = fReader->FillGen(fGen, fVz, Opt("etaMaxPart",5.5), true);
  const int ng = (int)fGen.size();
  fGPx.resize(ng); fGPy.resize(ng); fGPz.resize(ng); fGE.resize(ng);
  fGPure.assign(ng, -1); fNInter.assign(fNSlot, 0);
  fPure.clear(); fPureCons.clear();
  for (int g = 0; g < ng; ++g) {
    const double pt = fGen.pt[g], f = fGen.phi[g], m = fGen.mass[g];
    fGPx[g] = pt*cos(f); fGPy[g] = pt*sin(f); fGPz[g] = pt*sinh(fGen.eta[g]);
    fGE[g] = sqrt(fGPx[g]*fGPx[g] + fGPy[g]*fGPy[g] + fGPz[g]*fGPz[g] + m*m);
    ++fNInter[fGen.inter[g]];
  }
  std::vector<double> px, py, pz, E; std::vector< std::vector<int> > cons;
  for (int base = 0; base < ng; ) {
    const int inter = fGen.inter[base];
    int end = base; while (end < ng && fGen.inter[end] == inter) ++end;
    px.assign(fGPx.begin()+base, fGPx.begin()+end); py.assign(fGPy.begin()+base, fGPy.begin()+end);
    pz.assign(fGPz.begin()+base, fGPz.begin()+end); E.assign(fGE.begin()+base, fGE.begin()+end);
    std::vector<Jet> js = Cluster(px, py, pz, E, fR, fPtMinGen, cons);
    for (size_t k = 0; k < js.size(); ++k) {
      PJet p; p.pt = js[k].pt; p.y = js[k].y; p.eta = js[k].eta; p.phi = js[k].phi; p.m = js[k].m;
      p.ncon = js[k].ncon; p.inter = inter; p.vtx = -1; p.dr02 = p.kmain = p.dom = -1;
      p.topo = kLost; p.share = 0;
      p.rs.pt = p.rs.eta = p.rs.phi = p.rs.m = -1; p.rs.npf = -1;   // -1: no vertex, see RSeed
      p.rs.flead = p.rs.flinked = p.rs.fpuppi = -1; p.rs.vtx = -1; p.rs.pt0 = -1;
      const int ip = (int)fPure.size();
      fPure.push_back(p);
      std::vector<int> gi(cons[k].size());
      for (size_t m = 0; m < cons[k].size(); ++m) { gi[m] = base + cons[k][m]; fGPure[gi[m]] = ip; }
      fPureCons.push_back(gi);
    }
    base = end;
  }
}

// z per reco vertex, the interaction each one owns (unique, nearest within
// dzVertexGen, among the interactions with a particle in acceptance), and
// the usable vertex nearest to every interaction.
void GenSeed::ReconstructVertices()
{
  puppi::VertexZ(fPF.ztrk, fPF.vref, fPF.charged, fPF.eta, fPupCfg.etaTracker,
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
}

// What the linker did this event, before any vertex is looked at.
void GenSeed::LinkDiagnostics()
{
  using namespace genlink;
  const int ng = (int)fGen.size(), np = (int)fPF.size();
  const double pvz = fReader->PVz();
  for (int g = 0; g < ng; ++g) {
    const double ptg = fGen.pt[g], eg = fGen.eta[g];
    const int cls = DiagClassOf(fGen.pdg[g]), r = Reg(fabs(eg));
    const double lk = fLinks.linked[g] ? 1. : 0.;
    hlinkfrac[cls][r]->Fill(ptg, lk, ptg);
    pflink->Fill(eg, fGen.phi[g], lk, ptg);
  }
  for (size_t k = 0; k < fLinks.recs.size(); ++k) {
    const LinkRec &l = fLinks.recs[k];
    const double eg = fGen.eta[l.g], aeta = fabs(eg);
    const int r = Reg(aeta);
    if (l.pass == 1) hlinkdz[r]->Fill(l.dz);
    const int gc = GenClassOf(fGen.pdg[l.g], fGen.q[l.g]);
    const int pc = PFClassOf(fPF.pdg[l.i], fPF.q[l.i]);
    if (l.f > 0.9 && gc == kGGam && pc == kPGam && aeta < 1.4)
      hparallax->Fill(0.01*(fGen.z[l.g] - pvz)/(fLinkCfg.rEcal*cosh(eg)), fPF.eta[l.i] - eg);
    if (l.pass == 2 && gc == kGTrk && aeta < 1.3 && (pc == kPGam || pc == kPNH))
      hbend->Fill(1./fGen.pt[l.g], fGen.q[l.g]*dphi(fPF.phi[l.i], fGen.phi[l.g]));
  }
  for (int i = 0; i < np; ++i) {
    const int pc = PFClassOf(fPF.pdg[i], fPF.q[i]), r = Reg(fabs(fPF.eta[i]));
    horphan[pc][r]->Fill(fNPU, fLinks.orphan[i] ? 1. : 0., fPF.pt[i]);
    if (fLinks.pcount[i] == 1 && fLinks.pfrac[fLinks.pfirst[i]] > 0.9 && fLinks.Gpt[i] > 0)
      hpfratio[pc]->Fill(fPF.pt[i]/fLinks.Gpt[i]);
  }
  hlinkms->Fill(fLinkMs);
}

////////////////////////////////////////////////////////////////////////////
// One vertex hypothesis: PUPPI weights, the clustering of the weighted
// candidates with ghosts for the own interaction's unlinked particles, the
// seeds and reco-side topology of every jet; then, for every interaction
// whose nearest usable vertex this is (within dzRecoSeed), the gen side of
// its pure jets - dR partner, topology, mutual dominance - and their
// RecoJetSeeds with the weights of this vertex.
void GenSeed::ProcessVertex(int v)
{
  const int g = fGenOf[v];
  puppi::Weights(fPF, v, fZv[v], fPupCfg, fW);
  const int np = (int)fPF.size(), ng = (int)fGen.size();
  const genlink::Links &L = fLinks;

  std::vector<double> px, py, pz, E; std::vector<int> pfidx, ghost;
  px.reserve(np + 512); py.reserve(np + 512); pz.reserve(np + 512); E.reserve(np + 512);
  for (int i = 0; i < np; ++i) {
    if (fW[i] <= 0) continue;
    const double s = fW[i]*fPF.pt[i], m = fW[i]*fPF.mass[i];
    const double x = s*cos(fPF.phi[i]), y = s*sin(fPF.phi[i]), z = s*sinh(fPF.eta[i]);
    px.push_back(x); py.push_back(y); pz.push_back(z); E.push_back(sqrt(x*x + y*y + z*z + m*m));
    pfidx.push_back(i);
  }
  const int nreal = (int)px.size();
  if (g >= 0)
    for (int gi = 0; gi < ng; ++gi) {
      if (fGen.inter[gi] != g || L.linked[gi]) continue;
      px.push_back(kGhost*fGPx[gi]); py.push_back(kGhost*fGPy[gi]);
      pz.push_back(kGhost*fGPz[gi]); E.push_back(kGhost*fGE[gi]);
      ghost.push_back(gi);
    }
  hnghost->Fill((double)ghost.size());

  std::vector< std::vector<int> > cons;
  std::vector<Jet> js = Cluster(px, py, pz, E, fR, fPtMinCluster, cons);
  const int nj = (int)js.size();

  // ---- the seeds ---------------------------------------------------------
  const Acc zero = {0,0,0,0,0,0,0,0,0};
  std::vector<Acc> acc(nj, zero);
  std::vector< std::map<int,double> > fromPure(nj);
  std::vector<double> rattr(nj,0.), rorph(nj,0.);
  auto addGen = [&](int k, int gi, double w, bool ghostly) {
    if (w <= 0) return;
    Acc &a = acc[k];
    a.sx += w*fGPx[gi]; a.sy += w*fGPy[gi]; a.sz += w*fGPz[gi]; a.se += w*fGE[gi];
    if (!ghostly) { a.zx += w*fGPx[gi]; a.zy += w*fGPy[gi]; }
    const double p = w*fGen.pt[gi];
    if (g >= 0 && fGen.inter[gi] == g) a.own += p; else a.far += p;
    if (ghostly) a.unl += p;
    const int ip = fGPure[gi];
    if (ip >= 0) fromPure[k][ip] += p;
  };
  for (int k = 0; k < nj; ++k)
    for (size_t m = 0; m < cons[k].size(); ++m) {
      const int ic = cons[k][m];
      if (ic < nreal) {                          // a PF candidate routes its gen partners
        const int i = pfidx[ic]; const double w = fW[i];
        if (L.Gpt[i] > 0) rattr[k] += w*fPF.pt[i]; else rorph[k] += w*fPF.pt[i];
        for (int n = L.pfirst[i]; n < L.pfirst[i] + L.pcount[i]; ++n)
          addGen(k, L.pgen[n], w*L.pfrac[n], false);
      } else addGen(k, ghost[ic - nreal], 1., true);   // an unlinked own ghost
    }

  // ---- reco jets and their reco-side topology (a partition) --------------
  const int base = (int)fReco.size();
  std::map<int, std::map<int,double> > gaveTo;        // pure -> (jet -> pT)
  for (int k = 0; k < nj; ++k) {
    RJet r;
    r.pt = js[k].pt; r.y = js[k].y; r.eta = js[k].eta; r.phi = js[k].phi; r.m = js[k].m;
    r.ptcorr = r.pt*fJEC(r.pt, r.eta);
    r.vtx = v; r.gen = g; r.dompure = -1; r.ncopies = 1; r.rep = base + k;
    const Acc &a = acc[k];
    r.link = sqrt(a.sx*a.sx + a.sy*a.sy); r.linkz = sqrt(a.zx*a.zx + a.zy*a.zy);
    const double ssum = a.own + a.far, rsum = rattr[k] + rorph[k];
    r.fown = ssum > 0 ? a.own/ssum : 0.; r.ffar = ssum > 0 ? a.far/ssum : 0.;
    r.funl = ssum > 0 ? a.unl/ssum : 0.; r.forph = rsum > 0 ? rorph[k]/rsum : 0.;
    r.topo = kLost; r.pstar = -1; r.pshare = 0; r.bpair = 0;
    for (std::map<int,double>::const_iterator it = fromPure[k].begin(); it != fromPure[k].end(); ++it) {
      gaveTo[it->first][k] = it->second;
      if (it->second > r.pshare) { r.pshare = it->second; r.pstar = it->first; }
    }
    fReco.push_back(r);
  }
  std::vector<int> nin(nj, 0);
  for (int k = 0; k < nj; ++k) {
    RJet &r = fReco[base + k];
    for (std::map<int,double>::const_iterator it = fromPure[k].begin(); it != fromPure[k].end(); ++it)
      if (r.link > 0 && it->second > fShare*r.link) ++nin[k];
    if (nin[k] == 0) continue;                       // none
    int nout = 0, kmain = -1; double kmost = 0;
    const double ppt = fPure[r.pstar].pt;
    const std::map<int,double> &gt = gaveTo[r.pstar];
    for (std::map<int,double>::const_iterator it = gt.begin(); it != gt.end(); ++it) {
      if (it->second > fShare*ppt) ++nout;
      if (it->second > kmost) { kmost = it->second; kmain = it->first; }
    }
    if (nin[k] == 1) r.topo = (nout == 1) ? kOneOne : kSplit;
    else             r.topo = (nout == 1) ? kMerge  : kTangle;
    r.bpair = (kmain == k && r.pshare > fShare*ppt) ? 1 : 0;
  }

  // ---- gen side: the interactions whose nearest usable vertex this is -----
  const double dr2max = kDRMatch*kDRMatch;
  for (size_t p = 0; p < fPure.size(); ++p) {
    PJet &J = fPure[p];
    const int i = J.inter;
    if (fNearVtx[i] != v || fabs(fZv[v] - fVz[i]) >= fDzRecoSeed) continue;
    J.vtx = v;
    const double gpt = J.pt;
    double best = 1e9; int bk = -1;                  // classic: nearest jet inside dR < 0.2
    for (int k = 0; k < nj; ++k) {
      const double d2 = DeltaR2(J.eta, J.phi, js[k].eta, js[k].phi);
      if (d2 < best) { best = d2; bk = k; }
    }
    if (bk >= 0 && best < dr2max) J.dr02 = base + bk;
    int nr = 0, kmain = -1; double most = 0;         // who took it
    std::map<int, std::map<int,double> >::const_iterator git = gaveTo.find((int)p);
    if (git != gaveTo.end())
      for (std::map<int,double>::const_iterator it = git->second.begin(); it != git->second.end(); ++it) {
        if (it->second > fShare*gpt) ++nr;
        if (it->second > most) { most = it->second; kmain = it->first; }
      }
    const int npure = kmain >= 0 ? nin[kmain] : 0;
    if      (nr == 0) J.topo = kLost;
    else if (nr >= 2) J.topo = (npure >= 2) ? kTangle : kSplit;
    else              J.topo = (npure >= 2) ? kMerge  : kOneOne;
    if (kmain >= 0) { J.kmain = base + kmain; J.share = most; }
    // Mutual dominance is ONE relation seen from both ends, so dom is bpair:
    // J is k's p* and k is a bpair (k is J's main jet and took more than
    // shareFrac of J, AND J gave more than shareFrac of k's pT^link, nin > 0).
    // Without the last condition a seed that is mostly far pileup would be
    // J's partner here and a fake in fake/ and bridge/.
    if (kmain >= 0 && fReco[base + kmain].bpair && fReco[base + kmain].pstar == (int)p) {
      J.dom = base + kmain; fReco[base + kmain].dompure = (int)p;
    }
    RecoJetSeed((int)p, v);
  }
}

// The RecoJetSeed of pure jet J with the PUPPI weights of vertex v (fW), and
// rs0: the same candidates at weight s_iJ alone, clustered the same way -
// what the detector made of J before PUPPI decided anything.
void GenSeed::RecoJetSeed(int p, int v)
{
  PJet &J = fPure[p];
  const genlink::Links &L = fLinks;
  const std::vector<int> &gi = fPureCons[p];
  std::map<int,double> share;                      // PF -> sum_{g in J} f_gi pT_g
  double ptsum = 0, ptlinked = 0;
  for (size_t m = 0; m < gi.size(); ++m) {
    const int g = gi[m]; const double ptg = fGen.pt[g];
    ptsum += ptg;
    if (!L.linked[g]) continue;
    ptlinked += ptg;
    for (int n = L.gfirst[g]; n < L.gfirst[g] + L.gcount[g]; ++n) share[L.gpf[n]] += L.gfrac[n]*ptg;
  }
  RSeed &rs = J.rs;
  // A vertex was found: from here on "nothing" is 0, not -1 (see RSeed).
  rs.vtx = v; rs.pt = rs.eta = rs.phi = rs.m = 0; rs.npf = 0; rs.flead = rs.fpuppi = 0; rs.pt0 = 0;
  rs.flinked = ptsum > 0 ? ptlinked/ptsum : 0;
  std::vector<double> px, py, pz, E, qx, qy, qz, qE; std::vector< std::vector<int> > cons;
  auto push = [&](int i, double f, std::vector<double> &x, std::vector<double> &y,
                  std::vector<double> &z, std::vector<double> &e) {
    const double pt = f*fPF.pt[i], m = f*fPF.mass[i];
    const double a = pt*cos(fPF.phi[i]), b = pt*sin(fPF.phi[i]), c = pt*sinh(fPF.eta[i]);
    x.push_back(a); y.push_back(b); z.push_back(c); e.push_back(sqrt(a*a + b*b + c*c + m*m));
  };
  double sumSPt = 0, sumSWPt = 0;
  for (std::map<int,double>::const_iterator it = share.begin(); it != share.end(); ++it) {
    const int i = it->first;
    if (fSpt[i] <= 0) continue;
    // over the SCALAR attributed pT S_i: in [0, 1] and summing to at most 1
    // over all pure jets, so no candidate is credited twice
    const double s = it->second/fSpt[i], sw = s*fW[i];
    sumSPt += s*fPF.pt[i]; sumSWPt += sw*fPF.pt[i];
    if (s  > 0) push(i, s,  qx, qy, qz, qE);
    if (sw > 0) push(i, sw, px, py, pz, E);
  }
  rs.fpuppi = sumSPt > 0 ? sumSWPt/sumSPt : 0;
  if (!qx.empty()) {
    std::vector<Jet> j0 = Cluster(qx, qy, qz, qE, fR, 0., cons);
    if (!j0.empty()) rs.pt0 = j0[0].pt;
  }
  if (px.empty()) return;
  std::vector<Jet> js = Cluster(px, py, pz, E, fR, 0., cons);
  if (js.empty()) return;
  rs.pt = js[0].pt; rs.eta = js[0].eta; rs.phi = js[0].phi; rs.m = js[0].m; rs.npf = js[0].ncon;
  rs.flead = sumSWPt > 0 ? rs.pt/sumSWPt : 0;
}

////////////////////////////////////////////////////////////////////////////
// The dedup across vertex hypotheses, then every fill on the kept list and
// the pure list, and the tuple.
void GenSeed::FillEvent(int half)
{
  const int nr = (int)fReco.size();
  // ---- dedup: group the copies, keep the one with the largest own share ---
  std::vector<int> ord(nr); for (int k = 0; k < nr; ++k) ord[k] = k;
  std::sort(ord.begin(), ord.end(), [&](int a, int b){ return fReco[a].pt > fReco[b].pt; });
  const double drdup = Opt("dupDR",0.15), fdup = Opt("dupPtFrac",0.10);
  std::vector< std::vector<int> > grp;
  for (int o = 0; o < nr; ++o) {
    const int k = ord[o]; int m = -1;
    for (size_t q = 0; q < grp.size() && m < 0; ++q) {
      const RJet &s = fReco[grp[q][0]];
      if (DeltaR2(fReco[k].eta,fReco[k].phi,s.eta,s.phi) > drdup*drdup) continue;
      if (fabs(fReco[k].pt - s.pt) > fdup*std::max(fReco[k].pt, s.pt)) continue;
      m = (int)q;
    }
    if (m < 0) grp.push_back(std::vector<int>(1,k)); else grp[m].push_back(k);
  }
  std::vector<int> sel; sel.reserve(grp.size());
  for (size_t q = 0; q < grp.size(); ++q) {
    int keep = grp[q][0];                            // pT order breaks the ties
    for (size_t m = 1; m < grp[q].size(); ++m) if (fReco[grp[q][m]].fown > fReco[keep].fown) keep = grp[q][m];
    for (size_t m = 0; m < grp[q].size(); ++m) { fReco[grp[q][m]].rep = keep; fReco[grp[q][m]].ncopies = (int)grp[q].size(); }
    hncopies->Fill((double)grp[q].size());
    sel.push_back(keep);
  }
  const int ns = (int)sel.size();
  const int ipu = fNPU < 35 ? 0 : (fNPU <= 55 ? 1 : 2);

  // ---- the kept reco jets: spectra, fakes, responses, maps ----------------
  std::vector<char> paired(ns, 0);
  for (int s = 0; s < ns; ++s) {
    const RJet &r = fReco[sel[s]];
    const int iy = YBin(r.y), ie = jec::EtaBin(r.eta);
    const int ic = PtBinExt(r.ptcorr), il = PtBinExt(r.link), iz = PtBinExt(r.linkz);
    if (iy >= 0) {
      if (ic >= 0) { hspec_reco[iy][half]->Fill(r.ptcorr); hrtopo1[r.topo][iy]->Fill(r.ptcorr); }
      if (il >= 0) hspec_seed[iy][half]->Fill(r.link);
    }
    if (ie >= 0) {
      if (ic >= 0) {
        hrall[ie]->Fill(r.ptcorr);
        if (r.topo == kLost) hrnone[ie]->Fill(r.ptcorr);
        if (!r.bpair) hrunpaired[ie]->Fill(r.ptcorr);
      }
      if (il >= 0) {
        const double R = r.pt/r.link;
        hresp_link[ie]->Fill(r.link, R); hrtopo[r.topo][ie]->Fill(r.link, R); hresppu[ipu][ie]->Fill(r.link, R);
        hfar[ie]->Fill(r.link, r.ffar); hunl[ie]->Fill(r.link, r.funl); horph[ie]->Fill(r.link, r.forph);
      }
      if (iz >= 0) hresp_linkz[ie]->Fill(r.linkz, r.pt/r.linkz);
    }
    if (r.ptcorr > 5) hjet_reco->Fill(r.eta, r.phi);
    if (r.link > 5) presp_link->Fill(r.eta, r.phi, r.pt/r.link);
    porph->Fill(r.eta, r.phi, r.forph); pfar->Fill(r.eta, r.phi, r.ffar);
  }

  // ---- the pure jets: spectra, efficiencies, topology, RecoJetSeed --------
  const int npu = (int)fPure.size();
  for (int p = 0; p < npu; ++p) {
    const PJet &J = fPure[p];
    const int iy = YBin(J.y), ie = jec::EtaBin(J.eta), ib = PtBinExt(J.pt);
    if (ib >= 0 && iy >= 0) {
      hspec_gen[iy][half]->Fill(J.pt);
      // the gen-side topology is a statement about one vertex: a jet whose
      // interaction has none is not "lost", it was never looked at
      if (J.vtx >= 0) { htall[iy]->Fill(J.pt); htopo[J.topo][iy]->Fill(J.pt); }
    }
    if (J.pt > 5) hjet_gen->Fill(J.eta, J.phi);
    // The RecoJetSeed spectrum and map, binned in the pure jet's |y| but with
    // no cut on its pT or |eta|: a 1-2 GeV pure jet whose pT^rs lands at 3 GeV
    // is in the spectrum like any other, as in hseed and hreco.
    if (J.vtx >= 0 && J.rs.pt > 0) {
      if (iy >= 0 && PtBinExt(J.rs.pt) >= 0) hspec_rs[iy][half]->Fill(J.rs.pt);
      if (J.rs.pt > 5) hjet_rs->Fill(J.rs.eta, J.rs.phi);
    }
    if (ib < 0 || ie < 0 || J.vtx < 0) continue;
    hgall[ie]->Fill(J.pt);
    if (J.dr02 >= 0) { hgmat_dr[ie]->Fill(J.pt); hresp_dr[ie]->Fill(J.pt, fReco[J.dr02].pt/J.pt); }
    if (J.dom >= 0) hgmat_dom[ie]->Fill(J.pt);
    if (J.topo != kLost) hgmat_link[ie]->Fill(J.pt);
    // Before PUPPI: filled whenever anything of J was reconstructed, also
    // where PUPPI then kept none of it (rs.pt = 0).
    if (J.rs.pt0 > 0) hresp_rs0[ie]->Fill(J.pt, J.rs.pt0/J.pt);
    // The responses are CONDITIONAL on a RecoJetSeed (pT^rs > 0): the fraction
    // that has one is hgmat_rs/hgall, and the two are read together.
    if (J.rs.pt > 0) {
      const double R = J.rs.pt/J.pt;
      hgmat_rs[ie]->Fill(J.pt);
      hresp_rs[ie]->Fill(J.pt, R);
      if (J.pt > 5) presp_rs->Fill(J.eta, J.phi, R);
      // reco/pure = (rs/pure) x (reco/rs): the second factor is what pileup
      // and the clustering of the whole event add, taken on the kept copy of
      // J's mutual partner like every other fill here
      if (J.dom >= 0) hresp_reco_rs[ie]->Fill(J.pt, fReco[fReco[J.dom].rep].pt/J.rs.pt);
    }
  }

  // ---- the bridge seed -> pure over the mutual pairs ---------------------
  // A pure jet claimed by kept jets of two vertices (rare after the dedup)
  // goes to the one it gave most to, so that every pure jet is in the matrix
  // or in hbmiss exactly once; a pair across |y| bins is a miss and a fake,
  // so that the identities hold per |y| bin.
  {
    std::map<int,int> claim;
    for (int s = 0; s < ns; ++s) {
      const RJet &r = fReco[sel[s]];
      if (!r.bpair || r.pstar < 0) continue;
      std::map<int,int>::iterator it = claim.find(r.pstar);
      if (it == claim.end() || r.pshare > fReco[sel[it->second]].pshare) claim[r.pstar] = s;
    }
    std::vector<char> pairedP(npu, 0);
    for (std::map<int,int>::const_iterator it = claim.begin(); it != claim.end(); ++it) {
      const int p = it->first, s = it->second; const RJet &r = fReco[sel[s]];
      const int is = PtBinExt(r.link), ip = PtBinExt(fPure[p].pt), iys = YBin(r.y), iyp = YBin(fPure[p].y);
      if (is < 0 || ip < 0 || iys < 0 || iys != iyp) continue;
      hbridge[iys][half]->Fill(r.link, fPure[p].pt);
      paired[s] = 1; pairedP[p] = 1;
    }
    for (int s = 0; s < ns; ++s) {
      const RJet &r = fReco[sel[s]];
      const int is = PtBinExt(r.link), iy = YBin(r.y);
      if (is >= 0 && iy >= 0 && !paired[s]) hbfake[iy][half]->Fill(r.link);
    }
    for (int p = 0; p < npu; ++p) {
      const int iy = YBin(fPure[p].y), ib = PtBinExt(fPure[p].pt);
      if (iy >= 0 && ib >= 0 && !pairedP[p]) hbmiss[iy][half]->Fill(fPure[p].pt);
    }
  }

  // ---- the tuple ----------------------------------------------------------
  bRun = fReader->Run(); bLumi = fReader->Lumi(); bEvent = fReader->Event();
  bNPU = fNPU; bNVtx = (int)fUsable.size(); bPVz = fReader->PVz(); bLumiW = LumiPbInv(fNPU);
  bVtxZ.clear(); bVtxNtrk.clear(); bVtxGen.clear();
  for (size_t u = 0; u < fUsable.size(); ++u) {
    const int v = fUsable[u];
    bVtxZ.push_back(fZv[v]); bVtxNtrk.push_back(fNtrk[v]); bVtxGen.push_back(fGenOf[v]);
  }
  std::vector<int> gjIdx(npu, -1), rjIdx(nr, -1);      // stored indices, -1 below threshold
  int ngj = 0, nrj = 0;
  for (int p = 0; p < npu; ++p) if (fPure[p].pt > fPtStoreGen) gjIdx[p] = ngj++;
  for (int s = 0; s < ns; ++s) if (fReco[sel[s]].pt > fPtStoreReco) rjIdx[sel[s]] = nrj++;
  auto rj = [&](int k) { return k < 0 ? -1 : rjIdx[fReco[k].rep]; };   // through the kept copy
  auto vt = [&](int v) { return v < 0 ? -1 : fVtxIdx[v]; };
  bGjPt.clear(); bGjEta.clear(); bGjPhi.clear(); bGjMass.clear(); bGjShare.clear();
  bGjInter.clear(); bGjNcon.clear(); bGjVtx.clear(); bGjDr02.clear(); bGjTopo.clear(); bGjKmain.clear();
  bRsPt.clear(); bRs0Pt.clear(); bRsEta.clear(); bRsPhi.clear(); bRsMass.clear(); bRsFlead.clear(); bRsFlinked.clear(); bRsFpuppi.clear(); bRsNpf.clear();
  for (int p = 0; p < npu; ++p) {
    if (gjIdx[p] < 0) continue;
    const PJet &J = fPure[p];
    bGjPt.push_back(J.pt); bGjEta.push_back(J.eta); bGjPhi.push_back(J.phi); bGjMass.push_back(J.m);
    bGjInter.push_back(J.inter); bGjNcon.push_back(J.ncon); bGjVtx.push_back(vt(J.vtx));
    bGjDr02.push_back(rj(J.dr02)); bGjTopo.push_back(J.topo); bGjKmain.push_back(rj(J.kmain)); bGjShare.push_back(J.share);
    bRsPt.push_back(J.rs.pt); bRs0Pt.push_back(J.rs.pt0); bRsEta.push_back(J.rs.eta); bRsPhi.push_back(J.rs.phi); bRsMass.push_back(J.rs.m);
    bRsNpf.push_back(J.rs.npf); bRsFlead.push_back(J.rs.flead); bRsFlinked.push_back(J.rs.flinked); bRsFpuppi.push_back(J.rs.fpuppi);
  }
  bRjPt.clear(); bRjPtcorr.clear(); bRjEta.clear(); bRjPhi.clear(); bRjMass.clear(); bRjLink.clear(); bRjLinkz.clear();
  bRjFown.clear(); bRjFfar.clear(); bRjFunl.clear(); bRjForph.clear(); bRjPshare.clear();
  bRjVtx.clear(); bRjGen.clear(); bRjTopo.clear(); bRjPstar.clear(); bRjBpair.clear(); bRjNcopies.clear();
  for (int s = 0; s < ns; ++s) {
    const RJet &r = fReco[sel[s]];
    if (rjIdx[sel[s]] < 0) continue;
    bRjPt.push_back(r.pt); bRjPtcorr.push_back(r.ptcorr); bRjEta.push_back(r.eta); bRjPhi.push_back(r.phi); bRjMass.push_back(r.m);
    bRjVtx.push_back(vt(r.vtx)); bRjGen.push_back(r.gen);
    bRjLink.push_back(r.link); bRjLinkz.push_back(r.linkz);
    bRjFown.push_back(r.fown); bRjFfar.push_back(r.ffar); bRjFunl.push_back(r.funl); bRjForph.push_back(r.forph);
    bRjTopo.push_back(r.topo); bRjPstar.push_back(r.pstar < 0 ? -1 : gjIdx[r.pstar]); bRjPshare.push_back(r.pshare);
    bRjBpair.push_back(r.bpair); bRjNcopies.push_back(r.ncopies);
  }
  fTup->Fill();
}

////////////////////////////////////////////////////////////////////////////
void GenSeed::Loop()
{
  if (!fTree) { printf("GenSeed: no tree\n"); return; }
  fReader = new nanoaod::Reader(fTree);
  if (!fReader->Ok()) { printf("GenSeed: the tree lacks required branches, stopping\n"); return; }
  // A table that was asked for and does not load stops the job before any
  // file is written: carrying on would make a pass-2 file with ptcorr = raw.
  if (!fJECFile.empty() && !fJEC.Load(fJECFile.c_str())) {
    printf("GenSeed: cannot load the correction %s, stopping\n", fJECFile.c_str()); return; }
  fR            = Opt("R", 0.4);
  fPtMinCluster = Opt("ptMinCluster", 1.0);
  fPtMinGen     = Opt("ptMinGen", 1.0);
  fPtStoreReco  = Opt("ptStoreReco", 3.0);
  fPtStoreGen   = Opt("ptStoreGen", 3.0);
  fShare        = Opt("shareFrac", 1./3);
  fDzRecoSeed   = Opt("dzRecoSeed", 0.2);
  fLinkCfg.kernelScale = Opt("kernelScale", fLinkCfg.kernelScale);
  fLinkCfg.useDz       = int(Opt("useDz", fLinkCfg.useDz));
  fLinkCfg.useCharge   = int(Opt("useCharge", fLinkCfg.useCharge));
  fLinkCfg.parallaxRef = int(Opt("parallaxRef", fLinkCfg.parallaxRef));
  fLinkCfg.bendSign    = int(Opt("bendSign", fLinkCfg.bendSign));
  Book();
  if (!fOut) return;
  BookTree();

  const Long64_t nall  = fReader->GetEntries();
  const Long64_t first = std::max((Long64_t)0, (Long64_t)Opt("firstEntry", 0));
  const Long64_t last  = std::min(nall, (Long64_t)Opt("lastEntry", nall));
  const int maxv = int(Opt("maxVertices", -1));
  const Long64_t progressEvery = std::max((Long64_t)1, (Long64_t)Opt("progressEvery", 200));
  const double kv[kNKnob] = {fR, fPtMinCluster, fPtMinGen, fPtStoreReco, fPtStoreGen, Opt("etaMaxPart",5.5),
    fShare, Opt("dzVertexGen",0.05), fDzRecoSeed, Opt("minTrkPerVertex",3), Opt("dupDR",0.15), Opt("dupPtFrac",0.10),
    fLinkCfg.kernelScale, (double)fLinkCfg.useDz, (double)fLinkCfg.useCharge, (double)fLinkCfg.parallaxRef,
    (double)fLinkCfg.bendSign, fJEC.IsLoaded() ? 1. : 0., (double)maxv, (double)first, (double)last,
    (double)progressEvery, strcmp(ClusterBackend(),"FastJet") == 0 ? 1. : 0.};
  for (int i = 0; i < kNKnob; ++i) hopts->SetBinContent(i+1, kv[i]);
  hcount->Fill(5.5);                                          // jobs

  printf("GenSeed: entries [%lld, %lld) of %lld, %d interaction slots, %s clustering, R %g, ghosts x %.0e, "
         "reco from %g GeV, pure gen from %g GeV, share > %.3f, dzRecoSeed %g cm, correction %s\n",
         first, last, nall, fReader->NSlots(), ClusterBackend(), fR, kGhost, fPtMinCluster, fPtMinGen,
         fShare, fDzRecoSeed, fJEC.IsLoaded() ? fJECFile.c_str() : "none (pass 1)");
  printf("  linker: kernelScale %g useDz %d useCharge %d parallaxRef %d bendSign %d; puppi %s\n",
         fLinkCfg.kernelScale, fLinkCfg.useDz, fLinkCfg.useCharge, fLinkCfg.parallaxRef, fLinkCfg.bendSign,
         puppi::kVersion);

  // An empty range (more jobs than event pairs) writes a valid empty file with
  // its jobs counter, so that hadd and run.sh see the job as done.  It must
  // not reach TTreeReader, which reads [k, k) as "from k to the end".
  if (last <= first || !fReader->SetEntriesRange(first, last)) {
    printf("GenSeed: empty range [%lld, %lld), nothing to do\n", first, last); Write(); return; }
  TStopwatch sw; sw.Start(); TStopwatch swl;
  double linkTot = 0; Long64_t nproc = 0, nclu = 0, npure = 0, nkept = 0;
  while (fReader->Next()) {
    const Long64_t jentry = fReader->Entry();
    if (jentry >= last) break;                                // the range is ours, whatever the reader does
    hcount->Fill(4.5);
    const int half = int(jentry % 2);                         // on the GLOBAL entry
    hcount->Fill(half + 0.5);
    hcount->AddBinContent(3 + half, fReader->NPU());
    ++nproc;

    ReadGen();
    fReader->FillPF(fPF);
    puppi::Prepare(fPF, fPupCfg);
    ReconstructVertices();
    swl.Start();
    genlink::Link(fPF, fGen, fReader->PVz(), fLinkCfg, fLinks);
    fLinkMs = 1000.*swl.RealTime(); linkTot += fLinkMs;
    // S_i = sum_g f_gi pT_g, the scalar generated pT attributed to candidate
    // i: the denominator of the RecoJetSeed shares.  Links::Gpt is the
    // magnitude of the VECTOR sum, smaller for partners in different
    // directions, and would credit a shared candidate more than once.
    fSpt.assign(fPF.size(), 0.);
    for (int i = 0; i < (int)fPF.size(); ++i)
      for (int n = fLinks.pfirst[i]; n < fLinks.pfirst[i] + fLinks.pcount[i]; ++n)
        fSpt[i] += fLinks.pfrac[n]*fGen.pt[fLinks.pgen[n]];
    LinkDiagnostics();

    fReco.clear();
    for (size_t u = 0; u < fUsable.size() && (maxv < 0 || (int)u < maxv); ++u) { ProcessVertex(fUsable[u]); ++nclu; }
    FillEvent(half);
    npure += fPure.size(); nkept += bRjNcopies.size();

    if (nproc % progressEvery == 0) {
      const double t = sw.RealTime(); sw.Continue();
      ProcInfo_t pi; gSystem->GetProcInfo(&pi);
      printf("  %lld / %lld  (%.1f ms/event, link %.1f ms, %.0f ghosts/clustering, %.0f s to go, entry %lld, %.0f MB resident)\n",
             jentry - first + 1, last - first, 1000.*t/nproc, linkTot/nproc, hnghost->GetMean(),
             t*(last - jentry - 1)/nproc, jentry, pi.fMemResident/1024.);
      fflush(stdout);
    }
  }
  const double t = sw.RealTime();
  fNEvents = nproc;
  printf("GenSeed: done, %lld events in %.1f s: %.1f ms/event, of which Link %.1f ms; %.1f clusterings per event, "
         "%.0f own unlinked ghosts per clustering, %.1f pure gen jets and %.1f stored reco jets per event\n",
         nproc, t, nproc ? 1000.*t/nproc : 0., nproc ? linkTot/nproc : 0., nproc ? double(nclu)/nproc : 0.,
         hnghost->GetMean(), nproc ? double(npure)/nproc : 0., nproc ? double(nkept)/nproc : 0.);
  CheckIdentities();
  Write();
}

// The identities the histogram file promises, printed as the worst absolute
// violation, which must be zero:  bridge  hseed = hbfake + proj_x(hbridge),
// hgen = hbmiss + proj_y(hbridge) per |y| and half;  topology  sum of the
// gen-side classes = hall (both over the pure jets with a vertex) and of the
// reco-side classes = hreco_a + hreco_b per |y|;  efficiency  every
// numerator <= hgall per |eta| bin.
void GenSeed::CheckIdentities()
{
  const int nb = kNPtExt;
  double wb = 0, wt = 0, we = 0;
  printf("  identities (worst |lhs - rhs| over bins):\n  %-8s", "");
  for (int iy = 0; iy < kNY; ++iy) printf(" %8s", YTag(iy).c_str());
  printf("   worst\n  %-8s", "bridge");
  for (int iy = 0; iy < kNY; ++iy) {
    double wy = 0;
    for (int h = 0; h < 2; ++h)
      for (int i = 1; i <= nb; ++i) {
        double cs = 0, cp = 0;
        for (int k = 1; k <= nb; ++k) { cs += hbridge[iy][h]->GetBinContent(i,k); cp += hbridge[iy][h]->GetBinContent(k,i); }
        wy = std::max(wy, fabs(hspec_seed[iy][h]->GetBinContent(i) - hbfake[iy][h]->GetBinContent(i) - cs));
        wy = std::max(wy, fabs(hspec_gen[iy][h]->GetBinContent(i) - hbmiss[iy][h]->GetBinContent(i) - cp));
      }
    printf(" %8.1e", wy); wb = std::max(wb, wy);
  }
  printf("   %.1e %s\n  %-8s", wb, wb > 1e-6 ? "VIOLATED" : "ok", "topo");
  for (int iy = 0; iy < kNY; ++iy) {
    double wy = 0;
    for (int i = 1; i <= nb; ++i) {
      double sg = 0, sr = 0;
      for (int t = 0; t < kNTopo; ++t) { sg += htopo[t][iy]->GetBinContent(i); sr += hrtopo1[t][iy]->GetBinContent(i); }
      wy = std::max(wy, fabs(sg - htall[iy]->GetBinContent(i)));
      wy = std::max(wy, fabs(sr - hspec_reco[iy][0]->GetBinContent(i) - hspec_reco[iy][1]->GetBinContent(i)));
    }
    printf(" %8.1e", wy); wt = std::max(wt, wy);
  }
  printf("   %.1e %s\n  %-8s", wt, wt > 1e-6 ? "VIOLATED" : "ok", "eff");
  for (int ie = 0; ie < jec::kNEta; ++ie) {
    double wy = 0;
    TH1D *num[4] = {hgmat_dr[ie], hgmat_dom[ie], hgmat_rs[ie], hgmat_link[ie]};
    for (int n = 0; n < 4; ++n) for (int i = 1; i <= nb; ++i)
      wy = std::max(wy, std::max(0., num[n]->GetBinContent(i) - hgall[ie]->GetBinContent(i)));
    printf(" %8.1e", wy); we = std::max(we, wy);
  }
  printf("   %.1e %s\n", we, we > 1e-6 ? "VIOLATED" : "ok");
}

void GenSeed::Write()
{
  const double pa = hcount->GetBinContent(3), pb = hcount->GetBinContent(4);
  printf("  events %.0f + %.0f, interactions %.0f + %.0f, L = %.3g + %.3g pb^-1, entries %.0f, jobs %.0f, tuple %lld entries\n",
         hcount->GetBinContent(1), hcount->GetBinContent(2), pa, pb, LumiPbInv(pa), LumiPbInv(pb),
         hcount->GetBinContent(5), hcount->GetBinContent(6), fTup ? fTup->GetEntries() : 0);
  fOut->Write();
  fOut->Close();
  printf("  wrote %s\n", fOutName.c_str());
}
