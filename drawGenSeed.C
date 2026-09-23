// drawGenSeed.C - every plot, table and text file of the genseed package,
// from the histogram file that GenSeed.C writes.
//
//   root -l -b -q 'drawGenSeed.C+("rootfiles/GenSeed_v1.root","_v1")'
//
// Writes plots/*.pdf (cropped), doc/plots.tex, doc/tables.tex and
// text/jec<tag>.txt.  Everything is also printed as tables.
//
// WHY THIS EXISTS.  GenSeed.C produces two new collections and this macro is
// the test of both.  The routed LINK SEED gives every reconstructed jet a
// generated partner at generated scale, so that R = pT_raw / pT_link is a pure
// detector response; v4 showed that the seed follows PUPPI and R loses the
// right tail the v3 ghost seed had.  The RECOJETSEED does the opposite: for
// every pure generated jet it collects the PF candidates linked to its
// particles, weighted by their linked share and their PUPPI weight at the own
// vertex, and clusters them.  If that works, pT_rs / pT_pure is a response
// with the efficiency of the linker (a RecoJetSeed exists whenever PUPPI kept
// anything of the jet's reconstructed candidates, rs0 whenever anything was
// linked) and the Gaussianity of the classic dR < 0.2 pair (nothing from
// another interaction is in it).  Section 4 puts the four
// definitions - reco/link, reco/linkz, reco/pure by dR, RecoJetSeed/pure - on
// the same axes, and section 5 shows the distributions themselves.  Two more
// ride along where the file has them: rs0/pure, the RecoJetSeed before PUPPI
// (the same candidates at their share alone: what PUPPI took away is the gap
// to rs/pure), and reco/rs, the raw pT of the mutual-dominance reco partner
// over the RecoJetSeed (the pileup PUPPI kept and what the clustering added),
// so that reco/pure = rs/pure x reco/rs factorises the pairing.
//
// WHAT IS MEASURED, AND HOW.  A single core-fit protocol is used everywhere,
// ported from drawSeed4.C: the median first, then a Gaussian fit iterated to
// convergence in the ASYMMETRIC window [mu - 2 sigma, mu + 1 sigma] - v4's
// window, made for a response whose tail is on the right (pileup and merged
// neighbours add pT), and kept for EVERY definition so that the numbers
// compare with v4 on the same events - alongside the quantile widths
// (q50 - q16) and (q84 - q50) as the fit-free estimators.  A bin with fewer
// than 200 entries, a fit that does not converge, a core wider than 1.1 rms
// or a core mean outside [0.5, 1.3] x median FAILS and is dropped; it is
// never replaced by the rms, which a tail would have inflated.  The
// RecoJetSeed responses (rs, rs0) are the exception to the premise: they are
// built from a subset of the jet's own candidates, so R <= ~1 with a flat low
// side, their tail is on the LEFT, the window sits on it, and their core is
// undefined in most bins.  For them the quantile widths are the measure, and
// the generated text says so; a mirrored window would be a separately
// labelled metric, never a per-distribution switch (reco/link has sigma_L >
// sigma_R too, and flipping it would break the v4 comparison).  Every plot
// carries the thresholds line read from hist/hopts, so that no picture is
// looked at without knowing what was clustered and what was stored.
//
// RAW OR CORRECTED.  The reco-side axes (spectra, fakes, reco-side topology,
// the occupancy map) are in p_T^corr = p_T^raw x JEC.  In pass 1 no JEC is
// loaded (hopts jecLoaded = 0) and p_T^corr IS p_T^raw, so every label and
// caption then says p_T^raw / uncorrected; hist/jecFile names the table when
// there is one.  The responses themselves are always raw.
//
// THE PILEUP PROOF (section 6): a seed that knows what PUPPI kept gives the
// same R at any N_PU.  The link response is sliced in N_PU and the medians
// and core widths of the slices should lie on top of each other; the far,
// unlinked and orphan fractions of the seed say how much of it is not the
// own interaction.
//
// THE BRIDGE (section 8): whatever PUPPI decided is no longer in R but in the
// map from the seed to the pure generated jet, C(pT) = N_pure / N_seed per
// bin, measured per half a/b and decomposed into the fake fraction of the
// seeds, the miss fraction of the pure jets and the migration of the paired
// ones, with the identities hseed = hbfake + proj_x and hgen = hbmiss +
// proj_y checked bin by bin.
//
// THE LINKER DIAGNOSTICS (section 10) come last in the document but they are
// what everything above depends on: the linked fraction per particle class
// and detector region, the parallax slope (|slope| = 1 if PF points from the
// PV and the gen particle from its own vertex; v4 measured +0.98), the bend
// slope (-0.735 rad GeV at the ECAL for B = 3.8 T; v4 measured -0.97, the
// deposit sits deeper than the ECAL front face), the dz of the track links,
// the single-particle response per PF class and the orphan fraction against
// N_PU (1-2% in v4).
//
// THE DOCUMENT.  doc/plots.tex gets one figure environment per plot group,
// in the order of the sections, with a caption that says what is shown and
// quotes the key numbers; doc/tables.tex gets the efficiency, the response
// ladder, the topology fractions and the bridge as LaTeX tables.  Both are
// \input by doc/genseed.tex; the figure path is the macro \plotdir, which
// defaults to ../plots/ and can be redefined before the \input.  A number
// that does not exist - a failed core fit, an empty bin, a missing
// histogram - is written as "--" or "not defined", never as 0, and a
// caption that names rows ("first row: ...") gets those rows from TexFigure.
//
// AXES.  Every derived TH1 clones the axis of the histogram it is derived
// from (Like()); nothing here assumes a bin count or a bin edge.  Names are
// tried as given and under hist/, so a file with everything inside one
// directory works too.
//
// genseed - Copyright 2026 Mikko Voutilainen - Apache License 2.0
#include "tdrstyle_mod22.C"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TSystem.h>
#include <TError.h>
#include <TMath.h>
#include <fstream>
#include <memory>
#include <vector>
#include <set>
#include <string>
#include <cstring>
#include <initializer_list>

namespace {
  // -------------------------------------------------------------------
  // binning tags: the contract with GenSeed.C
  const int    kNEta = 10; const double kEtaW = 0.5;
  TString EtaTag(int i)   { return Form("e%02d", i); }
  TString EtaLabel(int i) { return Form("%.1f < |#eta| < %.1f", i*kEtaW, (i+1)*kEtaW); }
  TString EtaTex(int i)   { return Form("$%.1f < |\\eta| < %.1f$", i*kEtaW, (i+1)*kEtaW); }
  const int    kNY = 10; const double kYW = 0.5;
  TString YTag(int i)   { return Form("y%02d", i); }
  TString YLabel(int i) { return Form("%.1f < |y| < %.1f", i*kYW, (i+1)*kYW); }
  TString YTex(int i)   { return Form("$%.1f < |y| < %.1f$", i*kYW, (i+1)*kYW); }
  const char* HalfTag(int h) { return h == 0 ? "a" : "b"; }

  // the response definitions of sections 4 and 5.  The first kNRMain are in
  // every file.  kRS0 is the RecoJetSeed before PUPPI (the same candidates at
  // their share s_iJ alone); kRecoRS is the raw pT of the mutual-dominance
  // reco partner over pT^rs, against pT^pure: the pileup PUPPI kept plus what
  // the clustering added or took away, so reco/pure = rs/pure x reco/rs.  The
  // two are skipped quietly in a file that predates them (which has
  // resp/hresp_rs_corr instead: RS x JEC, never a closure, and no longer read)
  enum RDef { kLink = 0, kLinkZ, kDR, kRS, kRS0, kRecoRS, kNRDef };
  const int kNRMain = 4;
  const char* RTag(int d) {
    static const char *s[kNRDef] = {"link","linkz","dr","rs","rs0","recors"}; return s[d]; }
  const char* RHist(int d) {
    static const char *s[kNRDef] = {"resp/hresp_link_%s","resp/hresp_linkz_%s","resp/hresp_dr_%s",
                                    "resp/hresp_rs_%s","resp/hresp_rs0_%s","resp/hresp_reco_rs_%s"}; return s[d]; }
  const char* RLabel(int d) {
    static const char *s[kNRDef] = {"p_{T}^{raw} / p_{T}^{link} (routed seed)",
      "p_{T}^{raw} / p_{T}^{linkz} (linked only)", "p_{T}^{raw} / p_{T}^{pure} (#DeltaR < 0.2)",
      "p_{T}^{rs} / p_{T}^{pure} (RecoJetSeed)", "p_{T}^{rs0} / p_{T}^{pure} (RecoJetSeed, no PUPPI)",
      "p_{T}^{raw} / p_{T}^{rs} (dominance partner)"}; return s[d]; }
  const char* RShort(int d) {
    static const char *s[kNRDef] = {"reco/link","reco/linkz","reco/pure (dR)","RecoJetSeed/pure",
                                    "rs0/pure (no PUPPI)","reco/rs (partner)"}; return s[d]; }
  const char* RTex(int d) {
    static const char *s[kNRDef] = {"$\\pt^{\\mathrm{raw}}/\\pt^{\\mathrm{link}}$", "$\\pt^{\\mathrm{raw}}/\\pt^{\\mathrm{linkz}}$",
      "$\\pt^{\\mathrm{raw}}/\\pt^{\\mathrm{pure}}$ ($\\Delta R<0.2$)", "$\\pt^{\\mathrm{rs}}/\\pt^{\\mathrm{pure}}$",
      "$\\pt^{\\mathrm{rs0}}/\\pt^{\\mathrm{pure}}$", "$\\pt^{\\mathrm{raw}}/\\pt^{\\mathrm{rs}}$"}; return s[d]; }
  const char* RPt(int d) {
    static const char *s[kNRDef] = {"p_{T}^{link}","p_{T}^{linkz}","p_{T}^{pure}","p_{T}^{pure}","p_{T}^{pure}","p_{T}^{pure}"}; return s[d]; }
  // colour, marker and (for the distributions) line style; the two optional
  // ones are open markers and broken lines, so they read as the add-ons
  const int RColour[kNRDef] = {kBlack, kGreen+2, kBlue+1, kRed+1, kOrange+2, kMagenta+1};
  const int RMarker[kNRDef] = {kFullCircle, kFullCross, kFullSquare, kFullTriangleUp, kOpenTriangleUp, kOpenDiamond};
  const int RLine[kNRDef]   = {kSolid, kSolid, kSolid, kSolid, kDashed, kDashDotted};
  const char* RColTex(int d) {
    static const char *s[kNRDef] = {"black","green","blue","red","open orange","open magenta"}; return s[d]; }

  // efficiency definitions (section 2): eff/hgmat_<tag>_<e> over eff/hgall_<e>
  const int kNEff = 4;
  const char* EffTag(int k)   { static const char *s[kNEff] = {"dr","dom","link","rs"}; return s[k]; }
  const char* EffLabel(int k) { static const char *s[kNEff] = {"#DeltaR < 0.2 partner",
    "mutual-dominance partner", "linked (gen topology not lost)", "has a RecoJetSeed"}; return s[k]; }
  const char* EffTex(int k)   { static const char *s[kNEff] = {"$\\Delta R<0.2$","dominance","linked","RecoJetSeed"}; return s[k]; }
  const int EffColour[kNEff] = {kBlue+1, kMagenta+2, kBlack, kRed+1};
  const int EffMarker[kNEff] = {kFullSquare, kFullDiamond, kFullCircle, kFullTriangleUp};

  // raw or corrected: hopts jecLoaded, set once the knobs are read.  Every
  // reco-side axis is filled with pT^corr, which is pT^raw when no JEC was
  // loaded (pass 1), and must then say so
  bool gJecOn = false;
  const char* PtReco()     { return gJecOn ? "p_{T}^{corr}" : "p_{T}^{raw}"; }
  const char* PtRecoTex()  { return gJecOn ? "$\\pt^{\\mathrm{corr}}$" : "$\\pt^{\\mathrm{raw}}$"; }
  const char* RecoState()  { return gJecOn ? "corrected" : "uncorrected"; }

  // spectra (section 1): spec/h<tag>_<y>_<half>
  const int kNSpec = 4;
  const char* SpecTag(int k)   { static const char *s[kNSpec] = {"gen","reco","seed","rseed"}; return s[k]; }
  const char* SpecLabel(int k) { static const char *s[kNSpec] = {"pure gen jets, p_{T}^{pure}",
    "reco jets, p_{T}^{corr}", "link seeds, p_{T}^{link}", "RecoJetSeeds, p_{T}^{rs}"};
    return (k == 1 && !gJecOn) ? "reco jets, p_{T}^{raw} (uncorrected)" : s[k]; }
  const char* SpecTex(int k)   { static const char *s[kNSpec] = {"gen","reco","seed","RecoJetSeed"}; return s[k]; }
  const int SpecColour[kNSpec] = {kBlack, kRed+1, kBlue+1, kGreen+2};
  const int SpecMarker[kNSpec] = {kFullCircle, kFullSquare, kFullTriangleUp, kFullDiamond};

  // topology (section 7): a partition on each side
  const int kNTopo = 5;
  const char* GTopoTag(int t) { static const char *s[kNTopo] = {"1to1","split","merge","tangle","lost"}; return s[t]; }
  const char* RTopoTag(int t) { static const char *s[kNTopo] = {"1to1","split","merge","tangle","none"}; return s[t]; }
  const char* GTopoLabel(int t) { static const char *s[kNTopo] = {"1 #leftrightarrow 1",
    "split (1 pure #rightarrow 2+ reco)", "merge (2+ pure #rightarrow 1 reco)", "tangle", "lost (no reco jet)"}; return s[t]; }
  const char* RTopoLabel(int t) { static const char *s[kNTopo] = {"1 #leftrightarrow 1",
    "split (1 pure #rightarrow 2+ reco)", "merge (2+ pure #rightarrow 1 reco)", "tangle", "no pure contributor"}; return s[t]; }
  const int TopoColour[kNTopo] = {kBlack, kRed+1, kBlue+1, kGreen+2, kOrange+2};
  const int TopoMarker[kNTopo] = {kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullDiamond};

  // linker classes and regions (genlink.h DiagClassOf, GenSeed.C regions)
  const int kNCls = 6;
  const char* ClsTag(int c)   { static const char *s[kNCls] = {"pich","kp","ph","nh","v0","lep"}; return s[c]; }
  const char* ClsLabel(int c) { static const char *s[kNCls] = {"#pi^{#pm}", "K^{#pm}, p", "#gamma",
                                                                "n, K_{L}^{0}", "V^{0}", "e, #mu"}; return s[c]; }
  const char* ClsTex(int c)   { static const char *s[kNCls] = {"$\\pi^\\pm$", "$K^\\pm$, p", "$\\gamma$",
                                                                "n, $K^0_L$", "$V^0$", "e, $\\mu$"}; return s[c]; }
  const int kNReg = 3;
  const char* RegTag(int r)   { static const char *s[kNReg] = {"trk","tran","hf"}; return s[r]; }
  const char* RegLabel(int r) { static const char *s[kNReg] = {"|#eta| < 2.5", "2.5 < |#eta| < 3.0", "|#eta| > 3.0"}; return s[r]; }
  const char* RegTex(int r)   { static const char *s[kNReg] = {"$|\\eta|<2.5$", "$2.5<|\\eta|<3.0$", "$|\\eta|>3.0$"}; return s[r]; }
  const int RegEta[kNReg][2] = {{0,4},{5,5},{6,9}};
  const int kNPF = 6;
  const char* PFTag(int c)   { static const char *s[kNPF] = {"trk","gam","nh","hfh","hfe","other"}; return s[c]; }
  const char* PFLabel(int c) { static const char *s[kNPF] = {"charged hadron", "photon", "neutral hadron",
                                                              "HF hadron", "HF em", "other"}; return s[c]; }
  const int ClsColour[6] = {kBlack, kBlue+1, kRed+1, kGreen+2, kMagenta+2, kOrange+2};
  const int ClsMarker[6] = {kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullDiamond, kFullCross};

  const double kMinN = 200.;      // entries a distribution needs to be described
  const double kMinDen = 20.;     // denominator a fraction needs to be drawn
  const double kPtDrawHi = 100.;  // upper edge of the pT axes of the metric plots
  const double kPtSpecHi = 300.;  // and of the spectra
  const double kSigmaInel_mb = 80.0;  // 13.6 TeV minimum bias, the pileup library's value

  // -------------------------------------------------------------------
  // getters: a name is tried as given and under hist/; a missing histogram
  // is reported by name and counted
  TFile *gF = 0;
  int    gMissing = 0;
  TH1* GetH(const char *name, bool quiet = false) {
    TH1 *h = gF ? (TH1*)gF->Get(name) : 0;
    if (!h && gF) h = (TH1*)gF->Get(Form("hist/%s", name));
    if (!h) { if (!quiet) { printf("drawGenSeed: missing %s\n", name); ++gMissing; } return 0; }
    h->SetDirectory(0);
    return h;
  }
  TH1D*       Get1(const char *n) { return dynamic_cast<TH1D*>(GetH(n)); }
  TH2D*       Get2(const char *n) { return dynamic_cast<TH2D*>(GetH(n)); }
  TProfile*   GetP(const char *n) { return dynamic_cast<TProfile*>(GetH(n)); }
  TProfile2D* GetP2(const char *n){ return dynamic_cast<TProfile2D*>(GetH(n)); }
  // the text GenSeed.C keeps next to hopts (hist/puppiVersion, hist/jecFile);
  // "" when the file predates it
  TString GetText(const char *name) {
    TNamed *o = gF ? dynamic_cast<TNamed*>(gF->Get(name)) : 0;
    if (!o && gF) o = dynamic_cast<TNamed*>(gF->Get(Form("hist/%s", name)));
    return o ? TString(o->GetTitle()) : TString("");
  }
  // sum of Form(pat, EtaTag(ie)) over ie in [lo, hi]
  template<class T> T* SumEta(const char *pat, int lo, int hi, const char *name) {
    T *sum = 0;
    for (int ie = lo; ie <= hi; ++ie) {
      T *h = dynamic_cast<T*>(GetH(Form(pat, EtaTag(ie).Data())));
      if (!h) continue;
      if (!sum) { sum = (T*)h->Clone(name); sum->SetDirectory(0); } else sum->Add(h);
    }
    return sum;
  }
  template<class T> T* SumY(const char *pat, int lo, int hi, const char *name) {
    T *sum = 0;
    for (int iy = lo; iy <= hi; ++iy) {
      T *h = dynamic_cast<T*>(GetH(Form(pat, YTag(iy).Data())));
      if (!h) continue;
      if (!sum) { sum = (T*)h->Clone(name); sum->SetDirectory(0); } else sum->Add(h);
    }
    return sum;
  }
  // a + b of a per-half histogram
  template<class T> T* SumAB(const char *pat, const char *name) {
    T *a = dynamic_cast<T*>(GetH(Form(pat, "a"))), *b = dynamic_cast<T*>(GetH(Form(pat, "b")));
    if (!a || !b) return 0;
    T *s = (T*)a->Clone(name); s->SetDirectory(0); s->Add(b);
    return s;
  }

  // empty TH1D with the x axis of src: the ONLY way a pT axis is made here
  TH1D* Like(const TH1 *src, const char *name) {
    const TAxis *a = src->GetXaxis();
    TH1D *h = a->GetXbins()->GetSize() > 0
      ? new TH1D(name, "", a->GetNbins(), a->GetXbins()->GetArray())
      : new TH1D(name, "", a->GetNbins(), a->GetXmin(), a->GetXmax());
    h->SetDirectory(0); h->Sumw2();
    return h;
  }
  // where a falling spectrum puts the mean of a bin, near enough
  double GeoCentre(const TAxis *a, int b) { return sqrt(a->GetBinLowEdge(b)*a->GetBinUpEdge(b)); }
  TString BinLabel(const TAxis *a, int b) { return Form("%g-%g", a->GetBinLowEdge(b), a->GetBinUpEdge(b)); }
  TString BinTex(const TAxis *a, int b)   { return Form("%g--%g", a->GetBinLowEdge(b), a->GetBinUpEdge(b)); }
  // the bin of a pT axis that starts at pt (a hair inside, so an edge hits its own bin)
  int BinAt(const TH1 *h, double pt) { return h->GetXaxis()->FindFixBin(pt + 1e-3); }

  // A fraction of exactly 0 or 1 has a binomial error of 0, and ROOT's
  // painter skips a point with y = 0 and no error, so a well-measured "none"
  // would read as "no data".  The error is floored at 1/d (the step of one
  // entry), which also makes error > 0 the test of a defined bin (Val below).
  double BinomErr(double p, double d) { return std::max(sqrt(std::max(p*(1-p),0.)/d), 1./d); }
  TH1D* Fraction(const TH1 *num, const TH1 *den, const char *name, double minDen = kMinDen) {
    TH1D *h = Like(num, name);
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
      const double a = num->GetBinContent(i), d = den->GetBinContent(i);
      if (d < minDen) continue;
      const double p = a/d;
      h->SetBinContent(i, p);
      h->SetBinError(i, BinomErr(p, d));
    }
    return h;
  }
  // A number for a caption or a table cell, or "--" where it does not exist.
  // Every derived histogram here gives a filled bin a nonzero error
  // (Fraction and P(R>1.5) floor it, Ratio and the fits cannot give 0), so
  // error > 0 is "defined"; a profile is asked for its entries instead,
  // because a bin whose entries all agree (a linked fraction of exactly 1)
  // has error 0.  minEff: effective entries a profile cell needs (a weighted
  // profile's entries are not particles).
  bool Has(const TH1 *h, int b) { return h && b >= 1 && b <= h->GetNbinsX() && h->GetBinError(b) > 0; }
  TString Val(const TH1 *h, int b, const char *fmt = "%.2f", const char *miss = "--") {
    return Has(h, b) ? TString(Form(fmt, h->GetBinContent(b))) : TString(miss); }
  bool HasP(const TProfile *p, int b, double minEff = 0) {
    return p && b >= 1 && b <= p->GetNbinsX() && p->GetBinEntries(b) > 0 && p->GetBinEffectiveEntries(b) >= minEff; }
  TString ValP(const TProfile *p, int b, const char *fmt = "%.3f", double minEff = 0, const char *miss = "--") {
    return HasP(p, b, minEff) ? TString(Form(fmt, p->GetBinContent(b))) : TString(miss); }
  // x projection with the R under/overflow included: every jet counts
  TH1D* ProjX(TH2 *h, const char *name) { TH1D *p = h->ProjectionX(name, 0, -1); p->SetDirectory(0); return p; }
  TH1D* ProjY(TH2 *h, const char *name) { TH1D *p = h->ProjectionY(name, 0, -1); p->SetDirectory(0); return p; }
  TH1D* ProfX(TProfile *p, const char *name) { TH1D *h = p->ProjectionX(name); h->SetDirectory(0); return h; }
  // normalised copy for shape comparisons
  TH1D* Unit(TH1D *h, const char *name) {
    TH1D *u = (TH1D*)h->Clone(name); u->SetDirectory(0);
    const double n = u->Integral(0, u->GetNbinsX()+1);
    if (n > 0) u->Scale(1./n);
    return u;
  }
  // dN/dpT of a count histogram
  TH1D* Density(const TH1D *h, const char *name) {
    TH1D *d = Like(h, name);
    for (int b = 1; b <= d->GetNbinsX(); ++b) {
      const double n = h->GetBinContent(b), w = h->GetBinWidth(b);
      if (n <= 0) continue;
      d->SetBinContent(b, n/w); d->SetBinError(b, sqrt(n)/w);
    }
    return d;
  }
  // num/den with Poisson errors, drawn where both are populated
  TH1D* Ratio(const TH1D *num, const TH1D *den, const char *name) {
    TH1D *r = Like(num, name);
    for (int b = 1; b <= r->GetNbinsX(); ++b) {
      const double n = num->GetBinContent(b), d = den->GetBinContent(b);
      if (n <= 0 || d < kMinDen) continue;
      r->SetBinContent(b, n/d); r->SetBinError(b, n/d*sqrt(1./n + 1./d));
    }
    return r;
  }
  // open marker of the same shape, for the b half
  int OpenOf(int m) {
    switch (m) {
    case kFullCircle: return kOpenCircle;  case kFullSquare: return kOpenSquare;
    case kFullTriangleUp: return kOpenTriangleUp;  case kFullTriangleDown: return kOpenTriangleDown;
    case kFullDiamond: return kOpenDiamond;  case kFullCross: return kOpenCross;
    default: return kOpenCircle; }
  }

  // -------------------------------------------------------------------
  // knobs and the thresholds line every plot carries
  TH1   *gOpts = 0;
  double gJobs = 1;
  double Knob(const char *k, double def) {
    if (!gOpts) return def;
    const int b = gOpts->GetXaxis()->FindFixBin(k);
    return (b > 0 && b <= gOpts->GetNbinsX()) ? gOpts->GetBinContent(b)/gJobs : def;
  }
  // two short lines hanging from the top right corner inside the frame (align
  // 33: right, top), or stacked upwards from a bottom corner when up = true
  TString gThr1, gThr2;
  void Thresholds(double x = 0.94, double y = 0.925, int align = 33, double size = 0.019, bool up = false) {
    TLatex t; t.SetNDC(); t.SetTextSize(size); t.SetTextAlign(align);
    t.DrawLatex(x, up ? y + 1.3*size : y, gThr1);
    t.DrawLatex(x, up ? y : y - 1.3*size, gThr2);
  }
  void Note(double x, double y, const char *s, double size = 0.026, int align = 11, int colour = kBlack) {
    TLatex t; t.SetNDC(); t.SetTextSize(size); t.SetTextAlign(align); t.SetTextColor(colour); t.DrawLatex(x, y, s); }
  void UnitLine(double x1, double x2, double y = 1.) {
    TLine l; l.SetLineStyle(kDashed); l.SetLineColor(kGray+1); l.DrawLine(x1, y, x2, y); }
  TString LumiText(double pbinv) {
    if (pbinv >= 1)    return Form("%.3g pb^{-1}", pbinv);
    if (pbinv >= 1e-3) return Form("%.3g nb^{-1}", pbinv*1e3);
    return Form("%.3g #mub^{-1}", pbinv*1e6);
  }

  // -------------------------------------------------------------------
  // the PDF writer: plots/<name>.pdf, cropped to the ink with ghostscript
  // (croppdf.sh logic: ROOT writes an A4 MediaBox whatever the canvas size)
  TString gGs; int gNPdf = 0, gNCrop = 0;
  std::set<std::string> gSaved;
  void FindGs() {
    for (const char *p : {"/usr/local/bin/gs", "/opt/homebrew/bin/gs", "/usr/bin/gs"})
      if (!gSystem->AccessPathName(p, kExecutePermission)) { gGs = p; break; }
    if (gGs.IsNull()) {
      TString w = gSystem->GetFromPipe("command -v gs 2>/dev/null || true");
      w = w.Strip(TString::kBoth, '\n'); w = w.Strip(TString::kBoth);
      if (!w.IsNull()) gGs = w;
    }
    if (gGs.IsNull()) printf("drawGenSeed: ghostscript (gs) not found, the PDFs are left uncropped\n");
    else printf("drawGenSeed: cropping with %s\n", gGs.Data());
  }
  void Crop(const char *f) {
    if (gGs.IsNull()) return;
    const TString bb = gSystem->GetFromPipe(Form("%s -q -dBATCH -dNOPAUSE -sDEVICE=bbox '%s' 2>&1 | grep '^%%%%BoundingBox:' || true",
                                                 gGs.Data(), f));
    int x0, y0, x1, y1;
    // the line is "%%BoundingBox: x0 y0 x1 y1": two literal percent signs
    if (sscanf(bb.Data(), "%%%%BoundingBox: %d %d %d %d", &x0, &y0, &x1, &y1) != 4) { printf("  crop: no bbox for %s\n", f); return; }
    x0 -= 1; y0 -= 1; x1 += 1; y1 += 1;      // one point so anti-aliased edges are not shaved
    const int w = x1 - x0, h = y1 - y0;
    if (w <= 0 || h <= 0) { printf("  crop: empty bbox for %s\n", f); return; }
    const TString tmp = Form("%s.crop", f);
    const int st = gSystem->Exec(Form("%s -q -o '%s' -sDEVICE=pdfwrite -dCompatibilityLevel=1.5 -dDEVICEWIDTHPOINTS=%d -dDEVICEHEIGHTPOINTS=%d -dFIXEDMEDIA "
                                      "-c '<</PageOffset [%d %d]>> setpagedevice' -f '%s' 2>/dev/null",
                                      gGs.Data(), tmp.Data(), w, h, -x0, -y0, f));
    if (st == 0 && !gSystem->AccessPathName(tmp)) { gSystem->Rename(tmp, f); ++gNCrop; }
    else { gSystem->Unlink(tmp); printf("  crop FAILED: %s\n", f); }
  }
  const char* SavePdf(TCanvas *c, const char *name) {
    const TString f = Form("plots/%s.pdf", name);
    c->SaveAs(f);
    Crop(f);
    ++gNPdf; gSaved.insert(name);
    return name;
  }

  // a canvas for an eta-phi map: room for the palette on the right, the CMS
  // header placed by CMS_lumi from the margins it finds
  TCanvas* MapCanvas(const char *name, TH1 *h, const char *ztitle, double zlo = 0, double zhi = 0) {
    setTDRStyle();
    TCanvas *c = new TCanvas(name, name, 50, 50, 700, 600);
    c->SetLeftMargin(0.12); c->SetRightMargin(0.18); c->SetTopMargin(0.07); c->SetBottomMargin(0.13);
    h->GetXaxis()->SetTitle("#eta"); h->GetYaxis()->SetTitle("#phi");
    h->GetYaxis()->SetTitleOffset(0.95);
    h->GetZaxis()->SetTitle(ztitle); h->GetZaxis()->SetTitleOffset(1.35);
    h->GetZaxis()->SetTitleSize(0.045); h->GetZaxis()->SetLabelSize(0.04);
    if (zhi > zlo) { h->SetMinimum(zlo); h->SetMaximum(zhi); }
    h->Draw("COLZ");
    CMS_lumi(c, 8, 11);
    c->Update();
    return c;
  }

  // -------------------------------------------------------------------
  // THE CORE-FIT PROTOCOL (drawSeed4.C).  Everything that can be said about
  // one response distribution, computed the same way everywhere.
  struct Shape {
    double n, ovf;                    // entries (all, incl. R overflow), overflow fraction
    double med, mean, rms;            // location and width, moments and median
    double q16, q25, q75, q84;        // quantiles
    double sigL, sigR, bowley;        // (q50-q16)/q50, (q84-q50)/q50, (q75+q25-2q50)/(q75-q25)
    double pHi, pLo;                  // P(R > 1.5), P(R < 0.5)
    double muC, sigC, muCerr, sigCerr;// Gaussian core
    double fcore, tailR, tailL, ks;   // core fraction, beyond +-2 sigma_core, KS to the core Gaussian
    int    npass; bool conv, ok, core;// ok: n >= kMinN; core: the core fit passed
  };
  // fraction of the distribution above x, bin-linear at the boundary
  double FracAbove(const TH1D *p, double x) {
    const int nb = p->GetNbinsX(), b = p->FindFixBin(x);
    if (b > nb) return p->GetBinContent(nb+1);
    if (b < 1)  return p->Integral(0, nb+1);
    const double w = p->GetBinWidth(b);
    return p->Integral(b+1, nb+1) + p->GetBinContent(b)*(p->GetBinLowEdge(b)+w - x)/w;
  }
  int gNDesc = 0, gNConv = 0, gNWide = 0, gNMu = 0;   // protocol bookkeeping
  Shape Describe(TH1D *p)
  {
    Shape s; memset(&s, 0, sizeof(s));
    if (!p) return s;
    const int nb = p->GetNbinsX();
    s.n = p->Integral(0, nb+1);
    if (s.n < kMinN) return s;
    s.ovf = p->GetBinContent(nb+1)/s.n;
    double pr[5] = {0.16, 0.25, 0.5, 0.75, 0.84}, q[5];
    p->GetQuantiles(5, q, pr);
    s.q16 = q[0]; s.q25 = q[1]; s.med = q[2]; s.q75 = q[3]; s.q84 = q[4];
    s.mean = p->GetMean(); s.rms = p->GetRMS();
    s.sigL = s.med > 0 ? (s.med - s.q16)/s.med : 0;
    s.sigR = s.med > 0 ? (s.q84 - s.med)/s.med : 0;
    s.bowley = (s.q75 - s.q25) > 0 ? (s.q75 + s.q25 - 2*s.med)/(s.q75 - s.q25) : 0;
    s.pHi = FracAbove(p, 1.5)/s.n;
    s.pLo = 1. - FracAbove(p, 0.5)/s.n;
    s.ok = true;
    // the core: iterate the asymmetric window to convergence
    double mu = s.med, sig = 0.5*(s.q84 - s.q16);
    if (sig <= 0) sig = 0.1*std::max(s.med, 1e-3);
    TF1 fg(Form("fg_%s", p->GetName()), "gaus", 0, p->GetXaxis()->GetXmax());
    for (s.npass = 0; s.npass < 20; ++s.npass) {
      fg.SetRange(mu - 2*sig, mu + sig);
      fg.SetParameters(p->GetMaximum(), mu, sig);
      if (p->Fit(&fg, "QRN0L") != 0) break;
      const double mu1 = fg.GetParameter(1), sig1 = fabs(fg.GetParameter(2));
      if (sig1 <= 0) break;
      const bool c = fabs(sig1 - sig)/sig < 0.01;
      // move the window half way: an edge sitting on a bin boundary would
      // otherwise flip the same bin in and out of the fit for ever
      mu = c ? mu1 : 0.5*(mu + mu1); sig = c ? sig1 : 0.5*(sig + sig1);
      if (c) { s.conv = true; break; }
    }
    s.muC = mu; s.sigC = sig;
    s.muCerr = fg.GetParError(1); s.sigCerr = fg.GetParError(2);
    // A core wider than the whole distribution is a fit that ran away, not a
    // core.  The bound is 1.1 rms: once the response IS Gaussian,
    // sigma_core = rms is the right answer and must not be called a failure.
    s.core = s.conv && sig <= 1.1*s.rms && mu >= 0.5*s.med && mu <= 1.3*s.med;
    ++gNDesc; if (s.conv) ++gNConv;
    if (s.conv && sig > 1.1*s.rms) ++gNWide;
    if (s.conv && (mu < 0.5*s.med || mu > 1.3*s.med)) ++gNMu;
    if (!s.core) return s;
    s.fcore = fg.GetParameter(0)*sig*sqrt(2*M_PI)/p->GetBinWidth(1)/s.n;
    s.tailR = FracAbove(p, mu + 2*sig)/s.n;
    s.tailL = 1. - FracAbove(p, mu - 2*sig)/s.n;
    double cum = p->GetBinContent(0)/s.n, ks = 0;
    for (int i = 1; i <= nb; ++i) {
      cum += p->GetBinContent(i)/s.n;
      const double g = 0.5*(1 + TMath::Erf((p->GetXaxis()->GetBinUpEdge(i) - mu)/(sig*M_SQRT2)));
      ks = std::max(ks, fabs(cum - g));
    }
    s.ks = ks;
    return s;
  }
  Shape DescribeBin(TH2 *h, int b) {
    Shape s; memset(&s, 0, sizeof(s));
    if (!h || b < 1 || b > h->GetNbinsX()) return s;
    std::unique_ptr<TH1D> p(h->ProjectionY(Form("_py_%s_%d", h->GetName(), b), b, b));
    p->SetDirectory(0);
    return Describe(p.get());
  }
  const char* MetricHeader() {
    return "      n   ovf    med rms/md sigC/mu muC/md mean/md  sR/sL bowley  P>1.5  P<0.5  fcore  tailR  tailL     KS"; }
  TString MetricRow(const Shape &s) {
    if (!s.ok) return Form("%7.0f      -", s.n);
    TString r = Form("%7.0f %5.3f %6.3f %6.3f", s.n, s.ovf, s.med, s.med > 0 ? s.rms/s.med : 0);
    if (s.core) r += Form(" %7.3f %6.3f", s.sigC/s.muC, s.muC/s.med);
    else        r += Form(" %7s %6s", s.conv ? "FAIL" : "noconv", "-");
    r += Form(" %7.3f %6.3f %6.3f %6.3f %6.3f", s.med > 0 ? s.mean/s.med : 0,
              s.sigL > 0 ? s.sigR/s.sigL : 0, s.bowley, s.pHi, s.pLo);
    if (s.core) r += Form(" %6.3f %6.3f %6.3f %6.3f", s.fcore, s.tailR, s.tailL, s.ks);
    else        r += "      -      -      -      -";
    return r;
  }
  // the same row for LaTeX: n, median, sigma_core/mu, sigma_R/sigma_L, P(R>1.5), f_core
  TString MetricTex(const Shape &s) {
    if (!s.ok) return Form("%.0f & -- & -- & -- & -- & --", s.n);
    const TString rl = s.sigL > 0 ? Form("%.2f", s.sigR/s.sigL) : "--";
    const TString sc = s.core ? Form("%.3f", s.sigC/s.muC) : "fail", fc = s.core ? Form("%.2f", s.fcore) : "--";
    return Form("%.0f & %.3f & %s & %s & %.3f & %s", s.n, s.med, sc.Data(), rl.Data(), s.pHi, fc.Data());
  }

  // one metric of a (pT, R) histogram against pT, on the histogram's own axis
  enum Metric { kMMedian = 0, kMSigCore, kMSigRL, kMPHi, kMFCore, kNMetric };
  const char* MetricLabel(int m) {
    static const char *s[kNMetric] = {"median of R", "#sigma_{core} / #mu_{core}",
      "#sigma_{R} / #sigma_{L} (quantile widths)", "P(R > 1.5)", "f_{core}"}; return s[m]; }
  const char* MetricTag(int m) {
    static const char *s[kNMetric] = {"median","sigcore","sigrl","phi","fcore"}; return s[m]; }
  const char* MetricTex(int m) {
    static const char *s[kNMetric] = {"median of $R$", "$\\sigma_{\\mathrm{core}}/\\mu_{\\mathrm{core}}$",
      "$\\sigma_R/\\sigma_L$", "$P(R>1.5)$", "$f_{\\mathrm{core}}$"}; return s[m]; }
  TH1D* MetricVsPt(TH2 *h, const char *name, int m)
  {
    if (!h) return 0;
    TH1D *o = Like(h, name);
    for (int b = 1; b <= o->GetNbinsX(); ++b) {
      const Shape s = DescribeBin(h, b);
      if (!s.ok || s.med <= 0) continue;
      double v = 0, e = 0;
      const double emed = 1.2533*0.5*(s.q84 - s.q16)/sqrt(s.n);
      switch (m) {
      case kMMedian: v = s.med; e = emed; break;
      case kMSigRL:  if (s.sigL <= 0) continue;
                     v = s.sigR/s.sigL; e = v*1.8/sqrt(s.n); break;   // quantile SE, Gaussian approx.
      case kMPHi:    v = s.pHi; e = BinomErr(v, s.n); break;             // an exact 0 stays drawn
      case kMSigCore: if (!s.core) continue;
                     v = s.sigC/s.muC; e = sqrt(pow(s.sigCerr/s.muC,2) + pow(v*s.muCerr/s.muC,2)); break;
      case kMFCore:  if (!s.core) continue;
                     v = s.fcore; e = v*sqrt(pow(s.sigCerr/s.sigC,2) + 1./s.n); break;
      default: continue;
      }
      o->SetBinContent(b, v); o->SetBinError(b, e);
    }
    return o;
  }
  bool HasPoints(const TH1 *h) {
    if (!h) return false;
    for (int b = 1; b <= h->GetNbinsX(); ++b) if (h->GetBinError(b) > 0 || h->GetBinContent(b) != 0) return true;
    return false;
  }
  double AxLo(const TH1 *h) { return h->GetXaxis()->GetXmin(); }
  // the mean of y per x bin of a TH2, for a fit: bins with fewer than minN
  // entries are left empty, and the error never falls below what the y bin
  // width can resolve (a bin whose entries all sit in one y bin has a spread
  // of rounding noise, which would give it an infinite weight)
  TH1D* MeanY(TH2 *h2, const char *name, double minN = 10) {
    std::unique_ptr<TProfile> pr(h2->ProfileX(Form("_pr_%s", name))); pr->SetDirectory(0);
    TH1D *h = Like(h2, name);
    const double wy = h2->GetYaxis()->GetBinWidth(1);
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
      const double n = pr->GetBinEntries(b);
      if (n < minN) continue;
      h->SetBinContent(b, pr->GetBinContent(b));
      h->SetBinError(b, std::max(pr->GetBinError(b), wy/sqrt(12.*n)));
    }
    return h;
  }

  // -------------------------------------------------------------------
  // the LaTeX writers.  A figure lists only PDFs that were actually saved.
  std::ofstream gTexP, gTexT;
  TString TexEscape(const char *s) { TString t(s); t.ReplaceAll("_", "\\_"); return t; }
  void TexSection(const char *title, const char *text) {
    gTexP << "\n\\subsection{" << title << "}\n" << text << "\n";
  }
  // rows: optional sizes of the row groups over the INPUT list (e.g. {3,2,3}),
  // for a caption that says "first row: ..."; a new row starts at every group
  // boundary and a group wider than a line still wraps.  The caller counts
  // the groups from what it actually saved (or lists the missing names too:
  // they are dropped here), so the rows stay what the caption says when a
  // histogram is missing.  Without it the rows are whatever fits the width,
  // and the caption must read "left to right, top to bottom".  The glue
  // between panels is \hfil, the order of \centering's own, so that a short
  // row is centred; \hfill would outrank it and push the panels to the edges.
  void TexFigure(const std::vector<TString> &pdfs, const TString &caption, const char *label, double width = 0,
                 const std::vector<int> &rows = {}) {
    std::vector<TString> ok; std::vector<int> grp;
    for (size_t k = 0; k < pdfs.size(); ++k) {
      if (!gSaved.count(pdfs[k].Data())) continue;
      int g = 0, edge = 0;
      for (int r : rows) { edge += r; if ((int)k >= edge) ++g; }
      ok.push_back(pdfs[k]); grp.push_back(g);
    }
    if (ok.empty()) return;
    const size_t n = ok.size();
    if (width <= 0) width = n == 1 ? 0.55 : (n == 2 || n == 4) ? 0.48 : n <= 6 ? 0.32 : 0.24;
    const size_t perRow = std::max<size_t>(1, size_t(1./width + 1e-6));
    gTexP << "\\begin{figure}[htbp]\\centering\n";
    size_t inRow = 0;
    for (size_t k = 0; k < n; ++k) {
      ++inRow;
      const bool brk = k + 1 < n && (grp[k+1] != grp[k] || inRow == perRow);
      if (brk) inRow = 0;
      gTexP << Form("\\includegraphics[width=%.2f\\textwidth]{\\plotdir %s.pdf}%s\n", width, ok[k].Data(),
                    k + 1 == n ? "" : brk ? "\\\\" : "\\hfil");
    }
    gTexP << "\\caption{" << caption << "}\n\\label{fig:gs_" << label << "}\n\\end{figure}\n\n";
  }
  void TexTableBegin(const char *caption, const char *label, const char *cols, const char *header) {
    gTexT << "\\begin{table}[htbp]\\centering\\small\n\\caption{" << caption << "}\n\\label{tab:gs_" << label << "}\n"
          << "\\begin{tabular}{" << cols << "}\\hline\n" << header << " \\\\ \\hline\n";
  }
  void TexTableEnd() { gTexT << "\\hline\n\\end{tabular}\n\\end{table}\n\n"; }
}

void drawGenSeed(const char *fname = "rootfiles/GenSeed_v1.root", const char *tag = "_v1")
{
  gF = TFile::Open(fname);
  if (!gF || gF->IsZombie()) { printf("cannot open %s\n", fname); return; }
  setTDRStyle();
  gErrorIgnoreLevel = kWarning;   // no "pdf file created" between the tables
  gSystem->mkdir("plots", kTRUE); gSystem->mkdir("doc", kTRUE); gSystem->mkdir("text", kTRUE);
  FindGs();

  // ---- counters, header, knobs ------------------------------------------
  double events = 0, npu = 0;
  if (TH1 *hc = GetH("hcount")) {
    auto at = [hc](const char *l) { const int b = hc->GetXaxis()->FindFixBin(l); return b > 0 ? hc->GetBinContent(b) : 0.; };
    if (at("jobs") > 0) gJobs = at("jobs");
    events = at("events_a") + at("events_b"); npu = at("npu_a") + at("npu_b");
    printf("=== hcount: %.0f + %.0f events, %.0f entries, %g jobs, sum N_PU %.0f, <PU> = %.1f, L = %.3g ub-1\n",
           at("events_a"), at("events_b"), at("entries"), gJobs, npu, events > 0 ? npu/events : 0,
           npu/(kSigmaInel_mb*1e9)*1e6);
  }
  extraText = "Simulation"; extraText2 = "";
  lumi_136TeV = Form("%s, #LTPU#GT = %.0f", LumiText(npu/(kSigmaInel_mb*1e9)).Data(), events > 0 ? npu/events : 0.);
  gOpts = GetH("hopts");
  if (gOpts) {
    printf("=== knobs (hopts / jobs):");
    for (int b = 1; b <= gOpts->GetNbinsX(); ++b)
      printf(" %s=%g", gOpts->GetXaxis()->GetBinLabel(b), gOpts->GetBinContent(b)/gJobs);
    printf("\n");
  }
  const double R = Knob("R", 0.4), ptMinCluster = Knob("ptMinCluster", 1.), ptMinGen = Knob("ptMinGen", 1.);
  const double ptStoreReco = Knob("ptStoreReco", 3.), ptStoreGen = Knob("ptStoreGen", 3.), shareFrac = Knob("shareFrac", 1./3.);
  const double dzRecoSeed = Knob("dzRecoSeed", 0.2), jecLoaded = Knob("jecLoaded", 0.);
  gJecOn = jecLoaded > 0.5;
  if (jecLoaded > 1e-6 && jecLoaded < 1 - 1e-6)
    printf("drawGenSeed: WARNING jecLoaded = %.3f - jobs with and without a JEC were merged; labelled as %s\n",
           jecLoaded, RecoState());
  // what the file says about itself beyond the numbers (GenSeed.C writes both)
  const TString puppiVersion = GetText("puppiVersion"), jecFile = GetText("jecFile");
  printf("=== reco pT is %s (hopts jecLoaded = %g)%s%s%s%s\n", gJecOn ? "p_T^corr" : "p_T^raw, uncorrected", jecLoaded,
         puppiVersion.IsNull() ? "" : "; PUPPI ", puppiVersion.Data(), jecFile.IsNull() ? "" : "; jecFile ", jecFile.Data());
  // the share in words, for the generated text: it is a knob
  const bool shareThird = fabs(shareFrac - 1./3.) < 1e-3;
  const TString shareTex = shareThird ? TString("a third") : TString(Form("%.2f", shareFrac));
  gThr1 = Form("anti-k_{T} R = %g, reco jets p_{T}^{raw} > %g GeV (stored > %g GeV)", R, ptMinCluster, ptStoreReco);
  gThr2 = Form("pure gen jets p_{T} > %g GeV (stored > %g GeV), share > %s", ptMinGen, ptStoreGen,
               shareThird ? "1/3" : Form("%.2f", shareFrac));

  gTexP.open("doc/plots.tex"); gTexT.open("doc/tables.tex");
  gTexP << "% doc/plots.tex - generated by drawGenSeed.C from " << fname << "; do not edit\n"
        << "\\providecommand{\\plotdir}{../plots/}\n\\providecommand{\\pt}{p_{\\mathrm{T}}}\n";
  gTexT << "% doc/tables.tex - generated by drawGenSeed.C from " << fname << "; do not edit\n"
        << "\\providecommand{\\pt}{p_{\\mathrm{T}}}\n";
  TString thrTex = Form("Anti-$k_{\\mathrm{T}}$ $R=%g$; reconstructed jets are clustered from $\\pt^{\\mathrm{raw}}>%g$\\,GeV "
                        "and stored above %g\\,GeV, pure generated jets from %g\\,GeV and stored above %g\\,GeV; "
                        "the sample has %.0f events at $\\langle N_{\\mathrm{PU}}\\rangle=%.0f$.",
                        R, ptMinCluster, ptStoreReco, ptMinGen, ptStoreGen, events, events > 0 ? npu/events : 0.);
  thrTex += gJecOn ? Form(" The reconstructed jets are corrected with \\texttt{%s}.", TexEscape(jecFile.IsNull() ? "?" : jecFile.Data()).Data())
                   : TString(" No jet energy correction was loaded for this file, so the reconstructed-jet $\\pt$ on every "
                             "axis below is the uncorrected $\\pt^{\\mathrm{raw}}$.");
  if (!puppiVersion.IsNull()) thrTex += Form(" PUPPI weights from \\texttt{puppi.h} version \\texttt{%s}.", TexEscape(puppiVersion).Data());

  // ======================================================================
  // 1. spectra per |y|: gen, reco, seed, RecoJetSeed, and their ratio to gen
  // ======================================================================
  printf("\n=== 1. spectra: jets above 5 GeV per collection and |y| bin, and the ratio to gen at 5-6 / 10-12 / 20-24 GeV ===\n");
  printf("%-16s %9s %9s %9s %9s   %-23s %-23s %-23s\n", "|y| bin", "gen", "reco", "seed", "rseed", "reco/gen", "seed/gen", "rseed/gen");
  TexSection("Spectra", Form("Figure~\\ref{fig:gs_spectra} shows $\\mathrm{d}N/\\mathrm{d}\\pt$ of the pure generated jets, the "
             "reconstructed jets (%s, in %s, de-duplicated across vertices), the routed link seeds and the RecoJetSeeds, "
             "per $|y|$ bin, with the ratio of each collection to the generated one. The seeds count reconstructed jets, so "
             "seed/gen is the bridge $1/C$ of figure~\\ref{fig:gs_bridge}; RecoJetSeed/gen is the efficiency of the linker "
             "convolved with the response. %s", RecoState(), PtRecoTex(), thrTex.Data()));
  {
    std::vector<TString> figs; TString cap;
    for (int iy = 0; iy < kNY; ++iy) {
      TH1D *h[kNSpec] = {0};
      for (int k = 0; k < kNSpec; ++k) h[k] = SumAB<TH1D>(Form("spec/h%s_%s_%%s", SpecTag(k), YTag(iy).Data()), Form("sp%d%d", iy, k));
      if (!h[0]) continue;
      TH1D *dn[kNSpec] = {0}, *ratio[kNSpec] = {0};
      double ymax = 0;
      for (int k = 0; k < kNSpec; ++k) if (h[k]) {
        dn[k] = Density(h[k], Form("dn%d%d", iy, k)); ymax = std::max(ymax, dn[k]->GetMaximum());
        if (k) ratio[k] = Ratio(h[k], h[0], Form("rs%d%d", iy, k));
      }
      if (ymax <= 0) continue;
      const double xlo = AxLo(h[0]);
      TH1D *up = tdrHist(Form("fsu%d", iy), "dN / dp_{T} [GeV^{-1}]", ymax*1e-7, ymax*30, "p_{T} [GeV]", xlo, kPtSpecHi);
      TH1D *dw = tdrHist(Form("fsd%d", iy), "Ratio to gen", 0., 2.4, Form("p_{T}^{pure}, %s, p_{T}^{link}, p_{T}^{rs} [GeV]", PtReco()), xlo, kPtSpecHi);
      up->GetXaxis()->SetLabelSize(0.);
      std::unique_ptr<TCanvas> c(tdrDiCanvas(Form("c_sp%d", iy), up, dw, 8, 11));
      c->cd(1); gPad->SetLogx(); gPad->SetLogy();
      TLegend *leg = tdrLeg(0.50, 0.55, 0.94, 0.82);
      leg->SetTextSize(0.033); leg->SetHeader(YLabel(iy));
      for (int k = 0; k < kNSpec; ++k) if (dn[k]) {
        tdrDraw(dn[k], "Pz", SpecMarker[k], SpecColour[k], kSolid, -1, kNone, 0, 0.9);
        leg->AddEntry(dn[k], SpecLabel(k), "PL");
      }
      Thresholds(0.94, 0.925, 33, 0.024);
      c->cd(2); gPad->SetLogx();
      UnitLine(xlo, kPtSpecHi);
      for (int k = 1; k < kNSpec; ++k) if (ratio[k]) tdrDraw(ratio[k], "Pz", SpecMarker[k], SpecColour[k], kSolid, -1, kNone, 0, 0.9);
      fixOverlay();
      figs.push_back(SavePdf(c.get(), Form("spectrum_%s", YTag(iy).Data())));
      // the table row
      const int b5 = BinAt(h[0], 5);
      printf("%-16s", YLabel(iy).Data());
      for (int k = 0; k < kNSpec; ++k) printf(" %9s", h[k] ? Form("%.0f", h[k]->Integral(b5, h[k]->GetNbinsX()+1)) : "-");
      printf("  ");
      for (int k = 1; k < kNSpec; ++k) {
        TString cell;
        for (double pt : {5., 10., 20.}) {
          const int b = BinAt(h[0], pt);
          cell += ratio[k] && ratio[k]->GetBinError(b) > 0 ? Form(" %6.3f", ratio[k]->GetBinContent(b)) : "      -";
        }
        printf(" %-23s", cell.Data());
      }
      printf("\n");
      if (iy == 0) {
        TString c1;
        for (int k = 1; k < kNSpec; ++k) { const int b = BinAt(h[0], 5);
          c1 += Form("%s%s/gen %s", k > 1 ? ", " : "", SpecTex(k), Val(ratio[k], b).Data()); }
        cap = Form("At $|y|<0.5$ and 5--6\\,GeV: %s; %.0f generated jets above 5\\,GeV in that bin of $|y|$.", c1.Data(),
                   h[0]->Integral(b5, h[0]->GetNbinsX()+1));
      }
    }
    TexFigure(figs, Form("Spectra $\\mathrm{d}N/\\mathrm{d}\\pt$ of the pure generated jets (black), the %s reconstructed jets "
              "(red, in %s), the link seeds (blue) and the RecoJetSeeds (green), one panel per $|y|$ bin from $|y|<0.5$ to "
              "$4.5<|y|<5.0$ (left to right, top to bottom), with the ratio of each collection to the generated spectrum below. %s",
              RecoState(), PtRecoTex(), cap.Data()), "spectra");
  }

  // ======================================================================
  // 2. efficiencies vs pT in three regions
  // ======================================================================
  printf("\n=== 2. pure gen jets with a partner, per definition and region ===\n");
  TexSection("Pairing efficiency", Form("The fraction of pure generated jets that have a reconstructed partner, per pairing rule: "
             "the classic $\\Delta R<0.2$ match to the nearest reconstructed jet at the own vertex, the mutual-dominance "
             "partner (the pure jet gave the reconstructed jet more than %s of its seed $\\pt^{\\mathrm{link}}$, the "
             "reconstructed jet received more than %s of $\\pt^{\\mathrm{pure}}$, and each is the other's largest), the linker (the gen-side "
             "topology is not \\emph{lost}, i.e.\\ one reconstructed jet at its vertex received more than %s of its $\\pt$) "
             "and the existence of a RecoJetSeed with $\\pt>0$. The last two are the efficiencies of the two new collections. "
             "The denominator is every pure generated jet whose interaction has a usable vertex within "
             "\\texttt{dzRecoSeed}~$=%g$\\,cm.", shareTex.Data(), shareTex.Data(), shareTex.Data(), dzRecoSeed));
  {
    std::vector<TString> figs; TString cap;
    for (int r = 0; r < kNReg; ++r) {
      TH1D *all = SumEta<TH1D>("eff/hgall_%s", RegEta[r][0], RegEta[r][1], Form("gall%d", r));
      if (!all) continue;
      TH1D *eff[kNEff] = {0};
      for (int k = 0; k < kNEff; ++k) {
        TH1D *m = SumEta<TH1D>(Form("eff/hgmat_%s_%%s", EffTag(k)), RegEta[r][0], RegEta[r][1], Form("gmat%d%d", k, r));
        if (m) eff[k] = Fraction(m, all, Form("eff%d%d", k, r));
      }
      const double xlo = AxLo(all);
      TH1D *frame = tdrHist(Form("feff%d", r), "Pure gen jets with a reco partner", 0., 1.39, "p_{T}^{pure} [GeV]", xlo, kPtDrawHi);
      std::unique_ptr<TCanvas> c(tdrCanvas(Form("c_eff%d", r), frame, 8, 11, kSquare));
      c->SetLogx();
      UnitLine(xlo, kPtDrawHi);
      TLegend *leg = tdrLeg(0.45, 0.17, 0.94, 0.42);
      leg->SetTextSize(0.028); leg->SetHeader(RegLabel(r));
      for (int k = 0; k < kNEff; ++k) if (eff[k]) {
        tdrDraw(eff[k], "Pz", EffMarker[k], EffColour[k], kSolid, -1, kNone, 0, 0.9);
        leg->AddEntry(eff[k], EffLabel(k), "PL");
      }
      Thresholds();
      fixOverlay();
      figs.push_back(SavePdf(c.get(), Form("efficiency_%s", RegTag(r))));
      printf("  %s\n  %-10s %9s", RegLabel(r), "pT pure", "n");
      for (int k = 0; k < kNEff; ++k) printf("%8s", EffTag(k));
      printf("\n");
      if (r == 0) TexTableBegin(("Pairing efficiency of the pure generated jets at " + TString(RegTex(0)) +
                                 ", per $\\pt^{\\mathrm{pure}}$ bin: fraction with a $\\Delta R<0.2$ partner, a mutual-dominance "
                                 "partner, a linked reconstructed jet (topology not lost) and a RecoJetSeed.").Data(),
                                "efficiency", "lrrrrr", "$\\pt^{\\mathrm{pure}}$ [GeV] & $N$ & $\\Delta R<0.2$ & dominance & linked & RecoJetSeed");
      for (int b = 1; b <= all->GetNbinsX(); ++b) {
        if (all->GetXaxis()->GetBinLowEdge(b) > 60 || all->GetBinContent(b) < kMinDen) continue;
        printf("  %-10s %9.0f", BinLabel(all->GetXaxis(), b).Data(), all->GetBinContent(b));
        for (int k = 0; k < kNEff; ++k) printf("%8s", eff[k] ? Form("%.3f", eff[k]->GetBinContent(b)) : "-");
        printf("\n");
        if (r == 0) {
          gTexT << BinTex(all->GetXaxis(), b) << Form(" & %.0f", all->GetBinContent(b));
          for (int k = 0; k < kNEff; ++k) gTexT << (eff[k] ? Form(" & %.3f", eff[k]->GetBinContent(b)) : " & --");
          gTexT << " \\\\\n";
        }
      }
      if (r == 0) {
        TexTableEnd();
        const int b = BinAt(all, 5);
        cap = "At " + TString(RegTex(0)) + " and 5--6\\,GeV: ";
        for (int k = 0; k < kNEff; ++k) cap += Form("%s%s %s", k ? ", " : "", EffTex(k), Val(eff[k], b).Data());
        cap += ".";
      }
    }
    TexFigure(figs, "Fraction of pure generated jets with a reconstructed partner against $\\pt^{\\mathrm{pure}}$, for the "
              "$\\Delta R<0.2$ rule, mutual dominance, the linker (gen-side topology not lost) and the RecoJetSeed, at "
              "$|\\eta|<2.5$ (left), $2.5<|\\eta|<3.0$ (middle) and $|\\eta|>3.0$ (right). " + cap, "efficiency");
  }

  // ======================================================================
  // 3. fake rates vs pT: topology none, no bpair
  // ======================================================================
  printf("\n=== 3. kept reco jets without a pure contributor (none) and without a mutual-dominance partner (unpaired), vs %s ===\n",
         gJecOn ? "pT^corr" : "pT^raw (uncorrected)");
  TexSection("Fake rates", Form("The fraction of kept reconstructed jets, against their %s $\\pt$ (%s), to which no single "
             "pure generated jet, of any interaction, gives more than %s of the seed $\\pt^{\\mathrm{link}}$ (topology "
             "\\emph{none}: pileup that PUPPI kept, soft own particles outside every pure jet, or orphan PF), and the fraction "
             "without a mutual-dominance partner (\\emph{unpaired}, which adds the jets whose largest pure contributor has "
             "its main jet elsewhere or gave them less than %s of its own $\\pt$).",
             RecoState(), PtRecoTex(), shareTex.Data(), shareTex.Data()));
  {
    std::vector<TString> figs; TString cap;
    const char *fpat[2] = {"fake/hrnone_%s", "fake/hrunpaired_%s"};
    const char *flab[2] = {"no pure contributor (topology none)", "no mutual-dominance partner"};
    for (int r = 0; r < kNReg; ++r) {
      TH1D *all = SumEta<TH1D>("fake/hrall_%s", RegEta[r][0], RegEta[r][1], Form("rall%d", r));
      if (!all) continue;
      TH1D *fr[2] = {0};
      for (int k = 0; k < 2; ++k) {
        TH1D *m = SumEta<TH1D>(fpat[k], RegEta[r][0], RegEta[r][1], Form("rfk%d%d", k, r));
        if (m) fr[k] = Fraction(m, all, Form("fake%d%d", k, r));
      }
      const double xlo = AxLo(all);
      TH1D *frame = tdrHist(Form("ffk%d", r), "Fraction of kept reco jets", 0., 1.39, Form("%s [GeV]", PtReco()), xlo, kPtDrawHi);
      std::unique_ptr<TCanvas> c(tdrCanvas(Form("c_fk%d", r), frame, 8, 11, kSquare));
      c->SetLogx();
      TLegend *leg = tdrLeg(0.40, 0.66, 0.94, 0.84);
      leg->SetTextSize(0.028); leg->SetHeader(RegLabel(r));
      for (int k = 0; k < 2; ++k) if (fr[k]) {
        tdrDraw(fr[k], "Pz", k == 0 ? kFullCircle : kFullSquare, k == 0 ? kBlack : kRed+1, kSolid, -1, kNone, 0, 0.9);
        leg->AddEntry(fr[k], flab[k], "PL");
      }
      Thresholds();
      fixOverlay();
      figs.push_back(SavePdf(c.get(), Form("fakes_%s", RegTag(r))));
      printf("  %s\n  %-10s %9s %8s %8s\n", RegLabel(r), gJecOn ? "pT corr" : "pT raw", "n", "none", "unpaired");
      for (int b = 1; b <= all->GetNbinsX(); ++b) {
        if (all->GetXaxis()->GetBinLowEdge(b) > 60 || all->GetBinContent(b) < kMinDen) continue;
        printf("  %-10s %9.0f %8s %8s\n", BinLabel(all->GetXaxis(), b).Data(), all->GetBinContent(b),
               fr[0] ? Form("%.3f", fr[0]->GetBinContent(b)) : "-", fr[1] ? Form("%.3f", fr[1]->GetBinContent(b)) : "-");
      }
      if (r == 0) {
        const int b5 = BinAt(all, 5), b10 = BinAt(all, 10);
        cap = Form("At %s the fraction with no pure contributor is %s at 5--6\\,GeV and %s at 10--12\\,GeV; unpaired %s and %s.",
                   RegTex(0), Val(fr[0], b5, "%.3f").Data(), Val(fr[0], b10, "%.3f").Data(),
                   Val(fr[1], b5, "%.3f").Data(), Val(fr[1], b10, "%.3f").Data());
      }
    }
    TexFigure(figs, Form("Fake rates of the kept reconstructed jets against %s (%s): the fraction with no pure "
              "generated contributor (black) and the fraction without a mutual-dominance partner (red), at $|\\eta|<2.5$, "
              "$2.5<|\\eta|<3.0$ and $|\\eta|>3.0$. %s", PtRecoTex(), RecoState(), cap.Data()), "fakes");
  }

  // ======================================================================
  // 4. response summaries vs pT per |eta| bin, the ladder table, the JEC
  // ======================================================================
  TexSection("Response", "$R$ is the reconstructed over the generated $\\pt$ in bins of the generated one, for four "
             "definitions of the pair: the reconstructed jet over its routed link seed (everything linked into the jet, at "
             "generated scale, weighted by PUPPI, plus the unlinked own ghosts), over the linked-only seed, the classic "
             "$\\Delta R<0.2$ pair over the pure generated jet, and the RecoJetSeed over the pure generated jet. Two more "
             "ride along where the file has them: the RecoJetSeed before PUPPI, $\\pt^{\\mathrm{rs0}}/\\pt^{\\mathrm{pure}}$ "
             "(the same candidates at their share alone, filled whenever anything of the jet was linked), and the factor "
             "from the RecoJetSeed to its mutual-dominance reconstructed partner, $\\pt^{\\mathrm{raw}}/\\pt^{\\mathrm{rs}}$ "
             "against $\\pt^{\\mathrm{pure}}$ (the pileup PUPPI kept and what the clustering added), so that "
             "reco/pure $=$ rs/pure $\\times$ reco/rs. Each metric "
             "is one point per $\\pt$ bin: the median; $\\sigma_{\\mathrm{core}}/\\mu_{\\mathrm{core}}$ from a Gaussian fitted "
             "iteratively in $[\\mu-2\\sigma,\\mu+\\sigma]$, dropped when the bin has fewer than 200 entries, the fit does not "
             "converge, $\\sigma_{\\mathrm{core}}>1.1$\\,rms or $\\mu_{\\mathrm{core}}$ leaves $[0.5,1.3]\\times$median; the "
             "ratio of the quantile widths $(q_{84}-q_{50})/(q_{50}-q_{16})$, which is 1 for a symmetric distribution; the "
             "outer fraction $P(R>1.5)$; and the core fraction $f_{\\mathrm{core}}$, the integral of the fitted Gaussian over "
             "the entries, 1 for a Gaussian. The window is the v4 one and is kept for every definition, so that the numbers "
             "compare with v4 on the same events. The RecoJetSeed responses ($\\pt^{\\mathrm{rs}}$ and "
             "$\\pt^{\\mathrm{rs0}}$ over $\\pt^{\\mathrm{pure}}$) have their tail on the \\emph{left}: they are built from a "
             "subset of the jet's own candidates, so $R$ stays below about 1 with a flat low side, the window sits on that tail, and "
             "their core is often undefined under this protocol. A missing red or orange point in the core metrics is therefore "
             "the protocol, not the collection; for the RecoJetSeed the quantile widths and $P(R>1.5)$ are the relevant measures.");
  const int metrics[kNMetric] = {kMMedian, kMSigCore, kMSigRL, kMPHi, kMFCore};
  // the smallest y range of each metric; the frame opens up from it where
  // the points need that (below)
  const double mlo[kNMetric] = {0.5, 0., 0.5, 0., 0.}, mhi[kNMetric] = {1.5, 1.09, 5.4, 0.85, 2.2};
  const double ladderPt[4] = {5, 8, 12, 18};
  TString capResp;
  for (int ie = 0; ie < kNEta; ++ie) {
    const TString e = EtaTag(ie);
    TH2D *hr[kNRDef];
    int nany = 0;
    for (int d = 0; d < kNRDef; ++d) { hr[d] = d >= kNRMain ? dynamic_cast<TH2D*>(GetH(Form(RHist(d), e.Data()), true)) : Get2(Form(RHist(d), e.Data()));
                                       if (hr[d]) ++nany; }
    if (!nany) continue;
    // the ladder table at every eta bin; LaTeX for e00, every definition the
    // file has (at most 4 x 6 rows: one page)
    printf("\n=== 4. response ladder, %s ===\n%-10s %-22s %s\n", EtaLabel(ie).Data(), "pT bin", "definition", MetricHeader());
    if (ie == 0) TexTableBegin("Response ladder at $|\\eta|<0.5$: entries, median, core width, right/left quantile width ratio, "
                               "outer fraction and core fraction of each response definition in four $\\pt$ bins of the "
                               "denominator ($\\pt^{\\mathrm{link}}$, $\\pt^{\\mathrm{linkz}}$ or $\\pt^{\\mathrm{pure}}$; the "
                               "factor $\\pt^{\\mathrm{raw}}/\\pt^{\\mathrm{rs}}$ of the mutual-dominance partner, where present, "
                               "is binned in $\\pt^{\\mathrm{pure}}$ too). \\emph{fail}: the core fit fails the protocol, which "
                               "for the left-tailed RecoJetSeed responses is the rule rather than the exception.",
                               "ladder", "llrrrrrr", "$\\pt$ [GeV] & $R$ & $N$ & median & $\\sigma_{\\mathrm{core}}/\\mu$ & $\\sigma_R/\\sigma_L$ & $P(R>1.5)$ & $f_{\\mathrm{core}}$");
    for (double pt : ladderPt) {
      bool first = true;
      for (int d = 0; d < kNRDef; ++d) {
        if (!hr[d]) continue;
        const int b = BinAt(hr[d], pt);
        const Shape s = DescribeBin(hr[d], b);
        printf("%-10s %-22s %s\n", first ? BinLabel(hr[d]->GetXaxis(), b).Data() : "", RShort(d), MetricRow(s).Data());
        if (ie == 0) gTexT << (first ? BinTex(hr[d]->GetXaxis(), b) : TString("")) << " & " << RTex(d) << " & " << MetricTex(s) << " \\\\\n";
        first = false;
      }
      if (ie == 0 && pt != ladderPt[3]) gTexT << "\\hline\n";
    }
    if (ie == 0) TexTableEnd();
    // the metrics against pT, one canvas each
    std::vector<TString> figs;
    for (int m = 0; m < kNMetric; ++m) {
      TH1D *h[kNRDef] = {0}; int np = 0;
      for (int d = 0; d < kNRDef; ++d) {
        // reco/rs is a factor, not a response: where a reco jet is exactly its
        // RecoJetSeed (nothing from pileup survived) it has a spike at R = 1,
        // and a "core" fitted to that spike (sigma/mu ~ 0.04) is not a width
        if (d == kRecoRS && (metrics[m] == kMSigCore || metrics[m] == kMFCore)) continue;
        h[d] = MetricVsPt(hr[d], Form("rm%d%d%d", ie, m, d), metrics[m]); if (HasPoints(h[d])) ++np;
      }
      if (!np) continue;
      TH2D *any = 0; for (int d = 0; d < kNRDef; ++d) if (hr[d]) { any = hr[d]; break; }
      const double xlo = AxLo(any);
      // up to four: one column with the full labels; more: two columns of
      // the short ones, so the legend stays in the top quarter
      const bool twoCol = np > 4;
      // The frame: the default range, opened up where the points need it -
      // downwards for a low median (the RecoJetSeed in HF, ~0.25), upwards so
      // that the highest point stays under the legend (reco/rs reaches a
      // median of 2 and P(R>1.5) of 0.8 at 2.5 < |eta| < 3).  below: the
      // fraction of the frame height under the legend, with a little air.
      double ylo = mlo[m], yhi = mhi[m], dmin = 1e9, dmax = -1e9;
      for (int d = 0; d < kNRDef; ++d)
        for (int b = 1; h[d] && b <= h[d]->GetNbinsX(); ++b) {
          if (!Has(h[d], b) || h[d]->GetXaxis()->GetBinLowEdge(b) >= kPtDrawHi) continue;
          dmin = std::min(dmin, h[d]->GetBinContent(b) - h[d]->GetBinError(b));
          dmax = std::max(dmax, h[d]->GetBinContent(b) + h[d]->GetBinError(b));
        }
      const double below = twoCol ? 0.68 : 0.54;
      if (dmin < ylo) ylo = std::max(0., dmin - 0.05*(yhi - ylo));
      if (dmax > ylo + below*(yhi - ylo)) yhi = ylo + (dmax - ylo)/below;
      TH1D *frame = tdrHist(Form("fr%d%d", ie, m), MetricLabel(metrics[m]), ylo, yhi,
                            "p_{T}^{link}, p_{T}^{linkz} or p_{T}^{pure} [GeV]", xlo, kPtDrawHi);
      std::unique_ptr<TCanvas> c(tdrCanvas(Form("c_r%d%d", ie, m), frame, 8, 11, kSquare));
      c->SetLogx();
      if (metrics[m] == kMMedian || metrics[m] == kMSigRL || metrics[m] == kMFCore) UnitLine(xlo, kPtDrawHi);
      TLegend *leg = twoCol ? tdrLeg(0.42, 0.70, 0.94, 0.86) : tdrLeg(0.42, 0.58, 0.94, 0.86);
      leg->SetTextSize(twoCol ? 0.027 : 0.026); leg->SetHeader(EtaLabel(ie));
      if (twoCol) leg->SetNColumns(2);
      for (int d = 0; d < kNRDef; ++d) if (HasPoints(h[d])) {
        tdrDraw(h[d], "Pz", RMarker[d], RColour[d], kSolid, -1, kNone, 0, 0.9);
        leg->AddEntry(h[d], twoCol ? RShort(d) : RLabel(d), "PL");
      }
      Thresholds();
      fixOverlay();
      figs.push_back(SavePdf(c.get(), Form("resp_%s_%s", MetricTag(metrics[m]), e.Data())));
      if (ie == 0 && metrics[m] == kMSigCore) {
        // a failed core is "not defined", never 0 (a converged core always has an error)
        const int b = BinAt(any, 5);
        capResp = "At $|\\eta|<0.5$ and 5--6\\,GeV $\\sigma_{\\mathrm{core}}/\\mu$ is ";
        int nc = 0;
        for (int d = 0; d < kNRDef; ++d) if (hr[d] && d != kRecoRS)
          capResp += Form("%s%s (%s)", nc++ ? ", " : "", Val(h[d], b, "%.3f", "not defined").Data(), RTex(d));
        capResp += ".";
      }
    }
    TString who;
    for (int d = 0; d < kNRDef; ++d) if (hr[d]) {
      static const char *name[kNRDef] = {"reco/link", "reco/linkz", "reco/pure by $\\Delta R<0.2$", "RecoJetSeed/pure",
                                         "the RecoJetSeed before PUPPI", "reco/RecoJetSeed of the mutual-dominance partner"};
      who += Form("%s%s (%s)", who.IsNull() ? "" : ", ", name[d], RColTex(d));
    }
    { const Ssiz_t k = who.Last(','); if (k != kNPOS) who.Replace(k, 1, " and"); }
    TexFigure(figs, Form("Response metrics against the denominator $\\pt$ at %s: median, $\\sigma_{\\mathrm{core}}/\\mu_{\\mathrm{core}}$, "
                         "$\\sigma_R/\\sigma_L$, $P(R>1.5)$ and $f_{\\mathrm{core}}$ (left to right, top to bottom) for %s. Bins that "
                         "fail the core-fit protocol are absent from the core metrics; for the left-tailed RecoJetSeed responses "
                         "that is most bins. %s", EtaTex(ie).Data(), who.Data(), ie == 0 ? capResp.Data() : ""),
              Form("resp_%s", e.Data()));
  }
  // the JEC nodes from the median of R = pT_raw / pT_link -> text/jec<tag>.txt
  {
    const TString fjec = Form("text/jec%s.txt", tag);
    std::ofstream ojec(fjec.Data());
    ojec << "# genseed correction: 1 / median(pT_raw / pT_link) per pT_link bin, re-indexed by pT_raw = median x geometric bin centre\n"
         << "# only bins with n >= 200; the correction is flat below the first node and above the last (jec.h)\n"
         << "# ie pT_raw[GeV] correction\n";
    printf("\n=== 4b. JEC nodes from the median of R = pT_raw / pT_link (n >= %.0f) -> %s ===\n", kMinN, fjec.Data());
    std::vector<TGraphErrors*> gjec;
    TString capJec;
    for (int ie = 0; ie < kNEta; ++ie) {
      TH2D *h = Get2(Form(RHist(kLink), EtaTag(ie).Data()));
      if (!h) continue;
      TGraphErrors *g = new TGraphErrors(); g->SetName(Form("gjec%d", ie));
      double lo = 0, hi = 0;
      for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const Shape s = DescribeBin(h, b);
        if (!s.ok || s.med <= 0.05) continue;
        const double ptl = GeoCentre(h->GetXaxis(), b), praw = s.med*ptl;
        const double emed = 1.2533*0.5*(s.q84 - s.q16)/sqrt(s.n);
        g->SetPoint(g->GetN(), praw, 1./s.med); g->SetPointError(g->GetN()-1, 0, emed/(s.med*s.med));
        ojec << Form("%3d %12.4f %12.5f\n", ie, praw, 1./s.med);
        if (lo == 0) lo = praw; hi = praw;
      }
      ojec << Form("# %s: %d nodes, pT_raw %.2f - %.2f GeV, flat below/above\n", EtaTag(ie).Data(), g->GetN(), lo, hi);
      printf("  %-18s %2d nodes, pT_raw %.2f - %.2f GeV\n", EtaLabel(ie).Data(), g->GetN(), lo, hi);
      if (g->GetN()) gjec.push_back(g);
      if (ie == 0 && g->GetN()) { double x, y; g->GetPoint(0, x, y); capJec = Form("At $|\\eta|<0.5$ the first node is %.2f at $\\pt^{\\mathrm{raw}}=%.1f$\\,GeV and the table has %d nodes up to %.0f\\,GeV.", y, x, g->GetN(), hi); }
    }
    ojec.close();
    if (!gjec.empty()) {
      static const int ec[kNEta] = {kBlack, kBlue+2, kBlue-4, kCyan+2, kGreen+2, kSpring-1, kOrange+7, kRed+1, kMagenta+2, kGray+2};
      static const int em[kNEta] = {kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullDiamond,
                                    kOpenCircle, kOpenSquare, kOpenTriangleUp, kOpenDiamond, kOpenCross};
      TH1D *frame = tdrHist("fjec", "Correction 1 / median(p_{T}^{raw} / p_{T}^{link})", 0.7, 2.0, "p_{T}^{raw} [GeV]", 1., kPtDrawHi);
      std::unique_ptr<TCanvas> c(tdrCanvas("c_jec", frame, 8, 11, kSquare));
      c->SetLogx();
      UnitLine(1., kPtDrawHi);
      TLegend *leg = tdrLeg(0.42, 0.56, 0.94, 0.86);
      leg->SetTextSize(0.026); leg->SetNColumns(2);
      for (TGraphErrors *g : gjec) {
        const int ie = atoi(g->GetName() + 4);
        tdrDraw(g, "PLz", em[ie], ec[ie], kSolid, -1, kNone, 0, 0.9);
        leg->AddEntry(g, EtaLabel(ie), "PL");
      }
      Note(0.19, 0.20, "nodes at p_{T}^{raw} = median(R) #times #sqrt{p_{T,lo}^{link} p_{T,hi}^{link}}", 0.026);
      Thresholds();
      fixOverlay();
      SavePdf(c.get(), "jec");
      TexFigure({"jec"}, Form("The jet energy correction of pass 1, $1/\\mathrm{median}(\\pt^{\\mathrm{raw}}/\\pt^{\\mathrm{link}})$ per "
                              "$\\pt^{\\mathrm{link}}$ bin with at least 200 entries, re-indexed to $\\pt^{\\mathrm{raw}}$ at the median "
                              "times the geometric bin centre, one curve per $|\\eta|$ bin; written to \\texttt{text/jec%s.txt} and "
                              "read by \\texttt{jec.h}, flat outside the nodes. %s", TexEscape(tag).Data(), capJec.Data()), "jec");
    }
  }

  // ======================================================================
  // 5. response distributions in representative bins, four definitions
  // ======================================================================
  printf("\n=== 5. response distributions in representative bins ===\n");
  TexSection("Response distributions", "The distributions behind the metrics, unit-normalised, in four $\\pt$ bins at two "
             "$|\\eta|$ bins, with every definition of the file overlaid; the legend quotes $\\sigma_{\\mathrm{core}}/\\mu_{\\mathrm{core}}$ "
             "and $P(R>1.5)$ of each. A right tail is a denominator that is too small for the jet it sits under; the flat "
             "left side of the RecoJetSeed responses is the part of the jet that the reconstruction, or PUPPI, did not keep.");
  {
    const int distEta[2] = {0, 3};
    for (int ie : distEta) {
      std::vector<TString> figs; TString cap; bool have[kNRDef] = {false};
      for (double pt : ladderPt) {
        TH1D *u[kNRDef] = {0}; Shape sh[kNRDef]; double ymax = 0; TString blab, btex; int nu = 0;
        for (int d = 0; d < kNRDef; ++d) {
          TH2D *h = dynamic_cast<TH2D*>(GetH(Form(RHist(d), EtaTag(ie).Data()), true));
          if (!h) continue;
          const int b = BinAt(h, pt);
          if (b < 1 || b > h->GetNbinsX()) continue;
          blab = BinLabel(h->GetXaxis(), b); btex = BinTex(h->GetXaxis(), b);
          TH1D *p = h->ProjectionY(Form("d_%d_%g_%d", ie, pt, d), b, b); p->SetDirectory(0);
          sh[d] = Describe(p);
          if (!sh[d].ok) { delete p; continue; }
          // coarser bins for the eye; the metrics were computed on the fine ones
          const int nb = p->GetNbinsX(), rb = std::max(1, nb/100);
          if (rb > 1 && nb % rb == 0) p->Rebin(rb);
          u[d] = Unit(p, Form("u_%d_%g_%d", ie, pt, d)); ++nu; have[d] = true;
          ymax = std::max(ymax, u[d]->GetMaximum());
        }
        if (!nu) continue;
        // the legend grows with the number of curves; the frame top with it,
        // so that the peaks stay under the legend
        const double ytop = ymax*(nu > 4 ? 2000 : 400), ylegLo = 0.87 - 0.035*(nu + 1);
        TH1D *frame = tdrHist(Form("fd%d%g", ie, pt), "Fraction of jets", ymax*3e-4, ytop, "R", 0., 2.5);
        std::unique_ptr<TCanvas> c(tdrCanvas(Form("c_d%d%g", ie, pt), frame, 8, 11, kSquare));
        c->SetLogy();
        TLegend *leg = tdrLeg(0.37, ylegLo, 0.94, 0.87);   // clear of the CMS label
        leg->SetTextSize(0.021); leg->SetHeader(Form("%s, %s GeV", EtaLabel(ie).Data(), blab.Data()));
        printf("--- %s, %s GeV\n%-22s %s\n", EtaLabel(ie).Data(), blab.Data(), "definition", MetricHeader());
        for (int d = 0; d < kNRDef; ++d) if (u[d]) {
          tdrDraw(u[d], "HIST", kNone, RColour[d], RLine[d], -1, kNone, 0, 0, 2);
          const TString sc = sh[d].core ? Form("%.3f", sh[d].sigC/sh[d].muC) : "fail";
          leg->AddEntry(u[d], Form("%s: #sigma_{c}/#mu %s, P(R>1.5) %.3f", RShort(d), sc.Data(), sh[d].pHi), "L");
          printf("%-22s %s\n", RShort(d), MetricRow(sh[d]).Data());
        }
        TLine l; l.SetLineStyle(kDashed); l.SetLineColor(kGray+1); l.DrawLine(1, ymax*3e-4, 1, ytop);
        Thresholds();
        fixOverlay();
        figs.push_back(SavePdf(c.get(), Form("dist_%s_pt%02.0f", EtaTag(ie).Data(), pt)));
        if (pt == ladderPt[0]) {
          cap = Form("At %s\\,GeV $P(R>1.5)$ is ", btex.Data());
          int nc = 0;
          for (int d = 0; d < kNRDef; ++d) if (u[d]) cap += Form("%s%.3f (%s)", nc++ ? ", " : "", sh[d].pHi, RTex(d));
          const TString rsc = u[kRS] && sh[kRS].core ? Form("%.3f", sh[kRS].sigC/sh[kRS].muC) : "not defined (the fit fails the protocol)";
          cap += Form(", and $\\sigma_{\\mathrm{core}}/\\mu$ of the RecoJetSeed is %s.", rsc.Data());
        }
      }
      TString who;
      for (int d = 0; d < kNRDef; ++d) if (have[d]) {
        static const char *name[kNRDef] = {"reco/link (black)", "reco/linkz (green)", "reco/pure by $\\Delta R<0.2$ (blue)",
                                           "RecoJetSeed/pure (red)", "the RecoJetSeed before PUPPI (dashed orange)",
                                           "reco/RecoJetSeed of the mutual-dominance partner (dash-dotted magenta)"};
        who += Form("%s%s", who.IsNull() ? "" : ", ", name[d]);
      }
      { const Ssiz_t k = who.Last(','); if (k != kNPOS) who.Replace(k, 1, " and"); }
      TexFigure(figs, Form("Response distributions at %s in the bins 5--6, 8--10, 12--15 and 18--21\\,GeV of the denominator, "
                           "unit-normalised: %s. %s", EtaTex(ie).Data(), who.Data(), cap.Data()), Form("dist_%s", EtaTag(ie).Data()), 0.48);
    }
  }

  // ======================================================================
  // 6. the pileup proof and the seed fractions
  // ======================================================================
  printf("\n=== 6. pileup proof: reco/link by N_PU slice ===\n");
  TexSection("Pileup dependence", "The link response sliced by $N_{\\mathrm{PU}}$: if the seed follows what PUPPI kept, the "
             "median and the core width of the slices lie on top of each other. The second row gives the far (not the own "
             "interaction), unlinked (ghost) and orphan (PF without a gen partner) fractions of the seed and of the jet, "
             "against $\\pt^{\\mathrm{link}}$.");
  {
    struct Slice { const char *tag, *label; };
    const Slice puSlice[3] = {{"pu0","N_{PU} < 35"}, {"pu1","35 #leq N_{PU} #leq 55"}, {"pu2","N_{PU} > 55"}};
    const int sliceColour[3] = {kBlue+1, kGreen+2, kRed+1};
    const int sliceMarker[3] = {kFullSquare, kFullTriangleUp, kFullTriangleDown};
    std::vector<TString> figs; TString cap;
    for (int ie : {0, 3}) {
      const TString e = EtaTag(ie);
      TH2D *all = Get2(Form(RHist(kLink), e.Data()));
      if (!all) continue;
      const double xlo = AxLo(all);
      TH1D *up = tdrHist(Form("fpu%d", ie), MetricLabel(kMMedian), 0.55, 1.45, "p_{T}^{link} [GeV]", xlo, kPtDrawHi);
      TH1D *dw = tdrHist(Form("fpd%d", ie), MetricLabel(kMSigCore), 0., 0.79, "p_{T}^{link} [GeV]", xlo, kPtDrawHi);
      up->GetXaxis()->SetLabelSize(0.);
      std::unique_ptr<TCanvas> c(tdrDiCanvas(Form("c_pu%d", ie), up, dw, 8, 11));
      c->cd(1); gPad->SetLogx();
      UnitLine(xlo, kPtDrawHi);
      TLegend *leg = tdrLeg(0.45, 0.54, 0.94, 0.80);
      leg->SetTextSize(0.033); leg->SetHeader(Form("%s, %s", RLabel(kLink), EtaLabel(ie).Data()));
      TH1D *med[4] = {0}, *sig[4] = {0};
      med[3] = MetricVsPt(all, Form("pum%d3", ie), kMMedian);
      sig[3] = MetricVsPt(all, Form("pus%d3", ie), kMSigCore);
      tdrDraw(med[3], "Pz", kFullCircle, kBlack, kSolid, -1, kNone, 0, 0.9);
      leg->AddEntry(med[3], "all", "PL");
      double spread = 0;
      for (int s = 0; s < 3; ++s) {
        TH2D *hs = Get2(Form("resp/hresppu_link_%s_%s", puSlice[s].tag, e.Data()));
        if (!hs) continue;
        med[s] = MetricVsPt(hs, Form("pum%d%d", ie, s), kMMedian);
        sig[s] = MetricVsPt(hs, Form("pus%d%d", ie, s), kMSigCore);
        tdrDraw(med[s], "Pz", sliceMarker[s], sliceColour[s], kSolid, -1, kNone, 0, 0.9);
        leg->AddEntry(med[s], puSlice[s].label, "PL");
        for (int b = 1; b <= all->GetNbinsX(); ++b)
          if (all->GetXaxis()->GetBinLowEdge(b) <= 30 && med[s]->GetBinError(b) > 0 && med[3]->GetBinError(b) > 0)
            spread = std::max(spread, fabs(med[s]->GetBinContent(b) - med[3]->GetBinContent(b)));
      }
      Thresholds(0.94, 0.925, 33, 0.024);
      c->cd(2); gPad->SetLogx();
      for (int s = 3; s >= 0; --s) if (sig[s])
        tdrDraw(sig[s], "Pz", s == 3 ? kFullCircle : sliceMarker[s], s == 3 ? kBlack : sliceColour[s], kSolid, -1, kNone, 0, 0.9);
      fixOverlay();
      figs.push_back(SavePdf(c.get(), Form("pu_%s", e.Data())));
      printf("--- %s: median / sigma_core/mu, all and per N_PU slice; largest |median_slice - median_all| below 30 GeV: %.4f\n",
             EtaLabel(ie).Data(), spread);
      printf("%-10s %13s", "pT link", "all");
      for (int s = 0; s < 3; ++s) printf("%17s", puSlice[s].tag);
      printf("\n");
      for (int b = 1; b <= all->GetNbinsX(); ++b) {
        if (all->GetXaxis()->GetBinLowEdge(b) > 60 || med[3]->GetBinError(b) <= 0) continue;
        auto cell = [&](TH1D *m, TH1D *w) {
          return TString(m && m->GetBinError(b) > 0 ? Form("%.3f", m->GetBinContent(b)) : "    -") + "/" +
                 (w && w->GetBinError(b) > 0 ? Form("%.3f", w->GetBinContent(b)) : "  -  "); };
        printf("%-10s %s", BinLabel(all->GetXaxis(), b).Data(), cell(med[3], sig[3]).Data());
        for (int s = 0; s < 3; ++s) printf("      %s", cell(med[s], sig[s]).Data());
        printf("\n");
      }
      if (ie == 0) cap = Form("At $|\\eta|<0.5$ the medians of the three slices stay within %.3f of the inclusive one below 30\\,GeV.", spread);
    }
    // the far, unlinked and orphan fractions, in the second row
    const int nPu = (int)figs.size();
    for (int ie : {0, 3}) {
      const TString e = EtaTag(ie);
      const char *fp[3] = {"resp/hfar_%s", "resp/hunl_%s", "resp/horph_%s"};
      const char *fl[3] = {"far (not own interaction) fraction of the seed", "unlinked (ghost) fraction of the seed", "orphan fraction of the reco p_{T}"};
      TH1D *fr[3] = {0}; TProfile *pf[3] = {0};     // the profiles say which bins exist (HasP)
      for (int k = 0; k < 3; ++k) if ((pf[k] = GetP(Form(fp[k], e.Data())))) fr[k] = ProfX(pf[k], Form("frac%d%d", ie, k));
      TH1D *any = 0; for (int k = 0; k < 3; ++k) if (fr[k]) { any = fr[k]; break; }
      if (!any) continue;
      const double xlo = AxLo(any);
      TH1D *frame = tdrHist(Form("ffr%d", ie), "Fraction", 0., 1.09, "p_{T}^{link} [GeV]", xlo, kPtDrawHi);
      std::unique_ptr<TCanvas> c(tdrCanvas(Form("c_fr%d", ie), frame, 8, 11, kSquare));
      c->SetLogx();
      TLegend *leg = tdrLeg(0.35, 0.62, 0.94, 0.86);
      leg->SetTextSize(0.028); leg->SetHeader(Form("%s, %s", RLabel(kLink), EtaLabel(ie).Data()));
      for (int k = 0; k < 3; ++k) if (fr[k]) {
        tdrDraw(fr[k], "Pz", TopoMarker[k], TopoColour[k], kSolid, -1, kNone, 0, 0.9);
        leg->AddEntry(fr[k], fl[k], "PL");
      }
      Thresholds();
      fixOverlay();
      figs.push_back(SavePdf(c.get(), Form("fractions_%s", e.Data())));
      printf("--- %s: far / unlinked / orphan fraction vs pT_link\n", EtaLabel(ie).Data());
      for (int b = 1; b <= any->GetNbinsX(); ++b) {
        if (any->GetXaxis()->GetBinLowEdge(b) > 60 || !(HasP(pf[0], b) || HasP(pf[1], b) || HasP(pf[2], b))) continue;
        printf("%-10s", BinLabel(any->GetXaxis(), b).Data());
        for (int k = 0; k < 3; ++k) printf(" %7s", ValP(pf[k], b, "%.3f", 0, "-").Data());
        printf("\n");
      }
      if (ie == 0) { const int b = BinAt(any, 5);
        cap += Form(" At 5--6\\,GeV the far fraction of the seed is %s, the unlinked %s and the orphan fraction of the jet %s.",
                    ValP(pf[0], b).Data(), ValP(pf[1], b).Data(), ValP(pf[2], b).Data()); }
    }
    TexFigure(figs, "Top: median (upper panel) and $\\sigma_{\\mathrm{core}}/\\mu_{\\mathrm{core}}$ (lower panel) of "
              "$\\pt^{\\mathrm{raw}}/\\pt^{\\mathrm{link}}$ for all events (black) and in three $N_{\\mathrm{PU}}$ slices, at "
              "$|\\eta|<0.5$ (left) and $1.5<|\\eta|<2.0$ (right). Bottom: the far, unlinked and orphan fractions against "
              "$\\pt^{\\mathrm{link}}$ in the same two $|\\eta|$ bins. " + cap, "pileup", 0.48, {nPu, (int)figs.size() - nPu});
  }

  // ======================================================================
  // 7. topology fractions: gen side (pure jets) and reco side (kept reco jets)
  // ======================================================================
  printf("\n=== 7. topology fractions, |y| < 1.0 ===\n");
  TexSection("Topology", Form("Each pure generated jet of an interaction with a usable vertex within \\texttt{dzRecoSeed} is "
             "classified by how its $\\pt$ reached the reconstructed jets at that vertex (shares above %s of its "
             "$\\pt^{\\mathrm{pure}}$): one-to-one, split into two or more, merged with another pure jet into one, tangled, "
             "or lost; each kept reconstructed jet likewise by its pure contributors, with \\emph{none} where no pure jet "
             "gives it more than %s of its seed $\\pt^{\\mathrm{link}}$. The gen side is against $\\pt^{\\mathrm{pure}}$, "
             "the reco side against the jet's own %s $\\pt$ (%s). Both are partitions, so the fractions sum to one.",
             shareTex.Data(), shareTex.Data(), RecoState(), PtRecoTex()));
  {
    std::vector<TString> figs; TString cap;
    TexTableBegin("Topology fractions at $|y|<1.0$ in four $\\pt$ bins: gen side (pure jets, against $\\pt^{\\mathrm{pure}}$) "
                  "and reco side (kept reconstructed jets, against " + TString(PtRecoTex()) + ").", "topology", "llrrrrrr",
                  "side & $\\pt$ [GeV] & $N$ & 1$\\leftrightarrow$1 & split & merge & tangle & lost / none");
    for (int side = 0; side < 2; ++side) {
      TH1D *num[kNTopo] = {0}, *den = 0;
      for (int t = 0; t < kNTopo; ++t)
        num[t] = SumY<TH1D>(side == 0 ? Form("topo/h%s_%%s", GTopoTag(t)) : Form("topo/hr%s_%%s", RTopoTag(t)), 0, 1, Form("tp%d%d", side, t));
      if (side == 0) den = SumY<TH1D>("topo/hall_%s", 0, 1, "tpall");
      else for (int t = 0; t < kNTopo; ++t) if (num[t]) { if (!den) { den = (TH1D*)num[t]->Clone("tprall"); den->SetDirectory(0); } else den->Add(num[t]); }
      if (!den) continue;
      TH1D *fr[kNTopo] = {0};
      const double xlo = AxLo(den);
      TH1D *frame = tdrHist(Form("ft%d", side), side == 0 ? "Fraction of pure gen jets" : "Fraction of kept reco jets", 0., 1.59,
                            side == 0 ? "p_{T}^{pure} [GeV]" : Form("%s [GeV]", PtReco()), xlo, kPtDrawHi);
      std::unique_ptr<TCanvas> c(tdrCanvas(Form("c_t%d", side), frame, 8, 11, kSquare));
      c->SetLogx();
      UnitLine(xlo, kPtDrawHi);
      TLegend *leg = tdrLeg(0.40, 0.60, 0.94, 0.86);
      leg->SetTextSize(0.028); leg->SetHeader(side == 0 ? "gen side, |y| < 1.0" : "reco side, |y| < 1.0");
      for (int t = 0; t < kNTopo; ++t) if (num[t]) {
        fr[t] = Fraction(num[t], den, Form("tf%d%d", side, t));
        tdrDraw(fr[t], "Pz", TopoMarker[t], TopoColour[t], kSolid, -1, kNone, 0, 0.9);
        leg->AddEntry(fr[t], side == 0 ? GTopoLabel(t) : RTopoLabel(t), "PL");
      }
      Thresholds();
      fixOverlay();
      figs.push_back(SavePdf(c.get(), side == 0 ? "topo_gen" : "topo_reco"));
      printf("--- %s side\n%-10s %9s", side == 0 ? "gen" : "reco", side == 0 ? "pT pure" : (gJecOn ? "pT corr" : "pT raw"), "n");
      for (int t = 0; t < kNTopo; ++t) printf("%8s", side == 0 ? GTopoTag(t) : RTopoTag(t));
      printf("%8s\n", "sum");
      for (int b = 1; b <= den->GetNbinsX(); ++b) {
        if (den->GetXaxis()->GetBinLowEdge(b) > 60 || den->GetBinContent(b) < kMinDen) continue;
        printf("%-10s %9.0f", BinLabel(den->GetXaxis(), b).Data(), den->GetBinContent(b));
        double sum = 0;
        for (int t = 0; t < kNTopo; ++t) { const double v = fr[t] ? fr[t]->GetBinContent(b) : 0; sum += v; printf("%8.3f", v); }
        printf("%8.3f\n", sum);
      }
      for (double pt : ladderPt) {
        const int b = BinAt(den, pt);
        if (b < 1 || b > den->GetNbinsX()) continue;
        gTexT << (pt == ladderPt[0] ? (side == 0 ? "gen" : "reco") : "") << " & " << BinTex(den->GetXaxis(), b) << Form(" & %.0f", den->GetBinContent(b));
        for (int t = 0; t < kNTopo; ++t) gTexT << " & " << Val(fr[t], b, "%.3f");
        gTexT << " \\\\\n";
      }
      if (side == 0) gTexT << "\\hline\n";
      const int b = BinAt(den, 5);
      cap += Form("%s side at 5--6\\,GeV: 1$\\leftrightarrow$1 %s, split %s, merge %s, tangle %s, %s %s. ", side == 0 ? "Gen" : "Reco",
                  Val(fr[0], b).Data(), Val(fr[1], b).Data(), Val(fr[2], b).Data(), Val(fr[3], b).Data(),
                  side == 0 ? "lost" : "none", Val(fr[4], b).Data());
    }
    TexTableEnd();
    TexFigure(figs, Form("Topology fractions at $|y|<1.0$: gen side, the pure generated jets against $\\pt^{\\mathrm{pure}}$ (left), and "
              "reco side, the kept reconstructed jets against %s (right). %s", PtRecoTex(), cap.Data()), "topology");
  }

  // ======================================================================
  // 8. the bridge: C(pT) = N_pure / N_seed and its decomposition
  // ======================================================================
  printf("\n=== 8. bridge C(pT) = hgen / hseed per |y| bin and half, and the decomposition ===\n");
  TexSection("Bridge", "Whatever PUPPI decided is not in $R$ but in the map from the seed to the pure generated jet: "
             "$C(\\pt)=N_{\\mathrm{pure}}/N_{\\mathrm{seed}}$ per bin, measured in each half of the sample, and decomposed as "
             "$C = M\\,(1-f_{\\mathrm{fake}})/(1-f_{\\mathrm{miss}})$ with $M$ the migration of the mutual pairs "
             "(projection ratio of the pair matrix), $f_{\\mathrm{fake}}$ the fraction of seeds without a pure partner and "
             "$f_{\\mathrm{miss}}$ the fraction of pure jets without a seed. The identities "
             "$N_{\\mathrm{seed}} = N_{\\mathrm{fake}} + \\mathrm{proj}_x$ and $N_{\\mathrm{pure}} = N_{\\mathrm{miss}} + "
             "\\mathrm{proj}_y$ are checked bin by bin.");
  {
    std::vector<TString> figs; TString cap;
    // C per half, |y| < 0.5
    TH1D *pure[2] = {0}, *C[2] = {0};
    for (int h = 0; h < 2; ++h) pure[h] = Get1(Form("spec/hgen_y00_%s", HalfTag(h)));
    if (pure[0] && pure[1]) {
      const double xlo = AxLo(pure[0]);
      TH1D *frame = tdrHist("fbr", "C = N_{pure} / N_{seed} per p_{T} bin", 0.4, 2.2, "p_{T}^{pure} over p_{T}^{link} [GeV]", xlo, kPtSpecHi);
      std::unique_ptr<TCanvas> c(tdrCanvas("c_br", frame, 8, 11, kSquare));
      c->SetLogx();
      UnitLine(xlo, kPtSpecHi);
      TLegend *leg = tdrLeg(0.42, 0.68, 0.94, 0.86);
      leg->SetTextSize(0.028); leg->SetHeader(YLabel(0));
      for (int h = 0; h < 2; ++h) {
        TH1D *seed = Get1(Form("spec/hseed_y00_%s", HalfTag(h)));
        if (!seed) continue;
        C[h] = Like(pure[h], Form("C%d", h));
        for (int b = 1; b <= C[h]->GetNbinsX(); ++b) {
          const double p = pure[h]->GetBinContent(b), s = seed->GetBinContent(b);
          if (s < kMinDen || p <= 0) continue;
          C[h]->SetBinContent(b, p/s); C[h]->SetBinError(b, p/s*sqrt(1./p + 1./s));
        }
        tdrDraw(C[h], "Pz", h == 0 ? kFullCircle : kOpenCircle, kBlack, kSolid, -1, kNone, 0, 0.9);
        leg->AddEntry(C[h], Form("half %s", HalfTag(h)), "PL");
      }
      Thresholds();
      fixOverlay();
      figs.push_back(SavePdf(c.get(), "bridge_y00"));
    }
    // decomposition, a + b, per |y| bin (plot and table for y00)
    for (int iy = 0; iy < kNY; ++iy) {
      const TString y = YTag(iy);
      TH1D *pu = SumAB<TH1D>(Form("spec/hgen_%s_%%s", y.Data()), Form("bp%d", iy));
      TH1D *se = SumAB<TH1D>(Form("spec/hseed_%s_%%s", y.Data()), Form("bs%d", iy));
      TH1D *bf = SumAB<TH1D>(Form("bridge/hbfake_%s_%%s", y.Data()), Form("bf%d", iy));
      TH1D *bm = SumAB<TH1D>(Form("bridge/hbmiss_%s_%%s", y.Data()), Form("bm%d", iy));
      TH2D *br = SumAB<TH2D>(Form("bridge/hbridge_%s_%%s", y.Data()), Form("bb%d", iy));
      if (!pu || !se || !bf || !bm || !br) continue;
      TH1D *px = ProjX(br, Form("bpx%d", iy)), *py = ProjY(br, Form("bpy%d", iy));
      double worstS = 0, worstP = 0;
      for (int b = 1; b <= se->GetNbinsX(); ++b) {
        if (se->GetBinContent(b) > 0) worstS = std::max(worstS, fabs(bf->GetBinContent(b) + px->GetBinContent(b) - se->GetBinContent(b))/se->GetBinContent(b));
        if (pu->GetBinContent(b) > 0) worstP = std::max(worstP, fabs(bm->GetBinContent(b) + py->GetBinContent(b) - pu->GetBinContent(b))/pu->GetBinContent(b));
      }
      printf("%s: identities hseed = hbfake + proj_x worst %.2e, hgen = hbmiss + proj_y worst %.2e; paired %.0f of %.0f seeds, %.0f pure%s\n",
             YLabel(iy).Data(), worstS, worstP, px->Integral(), se->Integral(), pu->Integral(),
             (worstS > 1e-6 || worstP > 1e-6) ? "   IDENTITY BROKEN" : "");
      if (iy != 0) continue;
      TH1D *Cs = Like(pu, "dC"), *M = Like(pu, "dM"), *F = Like(pu, "dF"), *Mi = Like(pu, "dMi");
      printf("%-10s %8s %8s %8s %8s %8s\n", "pT bin", "C", "M", "1-ffake", "1/(1-fm)", "check");
      TexTableBegin(Form("The bridge at $|y|<0.5$: $C=N_{\\mathrm{pure}}/N_{\\mathrm{seed}}$ per $\\pt$ bin in the two halves of the sample, and "
                         "its decomposition on the summed sample into the migration $M$ of the mutual pairs, the fake factor "
                         "$1-f_{\\mathrm{fake}}$ and the miss factor $1/(1-f_{\\mathrm{miss}})$; the identities hold to %.1e and %.1e.", worstS, worstP),
                    "bridge", "lrrrrrr", "$\\pt$ [GeV] & $N_{\\mathrm{pure}}$ & $C_a$ & $C_b$ & $M$ & $1-f_{\\mathrm{fake}}$ & $1/(1-f_{\\mathrm{miss}})$");
      for (int b = 1; b <= pu->GetNbinsX(); ++b) {
        const double p = pu->GetBinContent(b), s = se->GetBinContent(b), x = px->GetBinContent(b), yy = py->GetBinContent(b);
        if (s < kMinDen || p < kMinDen || x <= 0) continue;
        const double cc = p/s, m = yy/x, ff = 1 - bf->GetBinContent(b)/s, fm = 1./(1 - bm->GetBinContent(b)/p);
        Cs->SetBinContent(b, cc); Cs->SetBinError(b, cc*sqrt(1./p + 1./s));
        M->SetBinContent(b, m); M->SetBinError(b, m*sqrt(1./yy + 1./x));
        F->SetBinContent(b, ff); F->SetBinError(b, sqrt(std::max(0., ff*(1-ff)/s)));
        Mi->SetBinContent(b, fm); Mi->SetBinError(b, fm*fm*sqrt(std::max(0., (1/fm)*(1-1/fm)/p)));
        if (pu->GetXaxis()->GetBinLowEdge(b) <= 120) {
          printf("%-10s %8.3f %8.3f %8.3f %8.3f %8.4f\n", BinLabel(pu->GetXaxis(), b).Data(), cc, m, ff, fm, cc - m*ff*fm);
          gTexT << BinTex(pu->GetXaxis(), b) << Form(" & %.0f & %s & %s & %.3f & %.3f & %.3f \\\\\n", p,
                   C[0] && C[0]->GetBinError(b) > 0 ? Form("%.3f", C[0]->GetBinContent(b)) : "--",
                   C[1] && C[1]->GetBinError(b) > 0 ? Form("%.3f", C[1]->GetBinContent(b)) : "--", m, ff, fm);
        }
        if (pu->GetXaxis()->GetBinLowEdge(b) == 5) cap = Form("At 5--6\\,GeV $C=%.3f$, $M=%.3f$, $1-f_{\\mathrm{fake}}=%.3f$, $1/(1-f_{\\mathrm{miss}})=%.3f$; "
                                                              "the identities hold to %.1e (seeds) and %.1e (pure).", cc, m, ff, fm, worstS, worstP);
      }
      TexTableEnd();
      const double xlo = AxLo(pu);
      TH1D *frame = tdrHist("fbd", "Factor", 0.4, 2.2, "p_{T}^{link} or p_{T}^{pure} [GeV]", xlo, kPtSpecHi);
      std::unique_ptr<TCanvas> c(tdrCanvas("c_bd", frame, 8, 11, kSquare));
      c->SetLogx();
      UnitLine(xlo, kPtSpecHi);
      TLegend *leg = tdrLeg(0.40, 0.62, 0.94, 0.86);
      leg->SetTextSize(0.028); leg->SetHeader(Form("%s, a + b", YLabel(0).Data()));
      tdrDraw(Cs, "Pz", kFullCircle, kBlack, kSolid, -1, kNone, 0, 0.9);        leg->AddEntry(Cs, "C = N_{pure} / N_{seed}", "PL");
      tdrDraw(M, "Pz", kFullSquare, kBlue+1, kSolid, -1, kNone, 0, 0.9);        leg->AddEntry(M, "M: paired, proj_{y} / proj_{x}", "PL");
      tdrDraw(F, "Pz", kFullTriangleUp, kRed+1, kSolid, -1, kNone, 0, 0.9);     leg->AddEntry(F, "1 - f_{fake} (unpaired seeds)", "PL");
      tdrDraw(Mi, "Pz", kFullTriangleDown, kGreen+2, kSolid, -1, kNone, 0, 0.9);leg->AddEntry(Mi, "1 / (1 - f_{miss}) (unpaired pure)", "PL");
      Note(0.19, 0.20, "C = M (1 - f_{fake}) / (1 - f_{miss})", 0.030);
      Thresholds();
      fixOverlay();
      figs.push_back(SavePdf(c.get(), "bridge_decomp_y00"));
    }
    TexFigure(figs, "The bridge at $|y|<0.5$. Left: $C=N_{\\mathrm{pure}}/N_{\\mathrm{seed}}$ per $\\pt$ bin, half a (filled) and "
              "half b (open). Right: its decomposition on the summed sample into the migration $M$ of the mutual pairs (blue), "
              "the fake factor $1-f_{\\mathrm{fake}}$ (red) and the miss factor $1/(1-f_{\\mathrm{miss}})$ (green). " + cap, "bridge");
  }

  // ======================================================================
  // 9. the eta-phi maps
  // ======================================================================
  printf("\n=== 9. eta-phi maps ===\n");
  TexSection("Maps", Form("Occupancy and response over the detector: jets per $(\\eta,\\phi)$ cell for the three collections "
             "above 5\\,GeV (the reconstructed jets in %s), the mean response of the link seed and of the RecoJetSeed, the "
             "orphan fraction of the reconstructed $\\pt$, the linked fraction of the generated particles and the far fraction "
             "of the seed. A structure in $\\phi$ is a detector effect; a structure in $\\eta$ alone follows the calorimeter "
             "boundaries.", PtRecoTex()));
  {
    // The z ranges hold every populated cell of the p1 file: COLZ does not
    // paint a cell below the minimum, so a range that starts above the data
    // (0.5 for a RecoJetSeed response of ~0.55 in the barrel and ~0.25 in HF)
    // shows holes that read as "no jets", and one that ends below it (the
    // link response reaches 2 in HF) saturates.
    struct Map { const char *name; TString z; double zlo, zhi; bool prof; const char *tex; };
    const Map maps[8] = {
      {"map/hjet_reco", Form("kept reco jets, %s > 5 GeV", PtReco()), 0, 0, false, "reconstructed jets"},
      {"map/hjet_gen",  "pure gen jets, p_{T}^{pure} > 5 GeV", 0, 0, false, "pure generated jets"},
      {"map/hjet_rs",   "RecoJetSeeds, p_{T}^{rs} > 5 GeV", 0, 0, false, "RecoJetSeeds"},
      {"map/presp_link","#LTp_{T}^{raw} / p_{T}^{link}#GT, p_{T}^{link} > 5 GeV", 0.3, 2.2, true, "link response"},
      {"map/presp_rs",  "#LTp_{T}^{rs} / p_{T}^{pure}#GT, p_{T}^{pure} > 5 GeV", 0., 1., true, "RecoJetSeed response"},
      {"map/porph",     "orphan fraction of the reco p_{T}", 0, 0.1, true, "orphan fraction"},
      {"map/pflink",    "linked fraction of the gen particles", 0, 1, true, "linked fraction"},
      {"map/pfar",      "far fraction of the seed", 0, 1, true, "far fraction of the seed"}};
    std::vector<TString> figs; TString cap; std::vector<int> rows(3, 0);
    for (const Map &m : maps) {
      TH1 *h = GetH(m.name);
      if (!h) continue;
      ++rows[&m - maps < 3 ? 0 : &m - maps < 5 ? 1 : 2];
      const char *base = strrchr(m.name, '/') ? strrchr(m.name, '/') + 1 : m.name;
      std::unique_ptr<TCanvas> c(MapCanvas(Form("c_%s", base), h, m.z.Data(), m.zlo, m.zhi));
      Thresholds(0.12, 0.004, 11, 0.017, true);
      figs.push_back(SavePdf(c.get(), Form("map_%s", base)));
      // the number: the mean over |eta| < 2.5 of a profile, the total of a count map
      double sum = 0, n = 0;
      TH2 *h2 = dynamic_cast<TH2*>(h);
      for (int i = 1; h2 && i <= h2->GetNbinsX(); ++i) {
        if (fabs(h2->GetXaxis()->GetBinCenter(i)) > 2.5) continue;
        for (int j = 1; j <= h2->GetNbinsY(); ++j) {
          if (m.prof) { const double w = ((TProfile2D*)h2)->GetBinEntries(h2->GetBin(i, j)); sum += h2->GetBinContent(i, j)*w; n += w; }
          else { sum += h2->GetBinContent(i, j); n = 1; }
        }
      }
      const TString num = m.prof ? (n > 0 ? TString(Form("%.3f", sum/n)) : TString("--")) : TString(Form("%.0f", sum));
      printf("  %-16s %s |eta| < 2.5: %s\n", base, m.prof ? "mean over" : "total in", num.Data());
      cap += Form("%s%s %s at $|\\eta|<2.5$: %s", cap.IsNull() ? "" : "; ", m.prof ? "mean" : "total", m.tex, num.Data());
    }
    // three rows as the caption says them: occupancy, responses, fractions
    TexFigure(figs, Form("Maps in $(\\eta,\\phi)$: kept reconstructed jets (%s), pure generated jets and RecoJetSeeds above "
              "5\\,GeV (first row), the mean link response and the mean RecoJetSeed response (second row; $z$ from %g to %g "
              "and from %g to %g), the orphan fraction of the reconstructed $\\pt$, the linked fraction of the generated "
              "particles and the far fraction of the seed (third row). %s.", PtRecoTex(), maps[3].zlo, maps[3].zhi,
              maps[4].zlo, maps[4].zhi, cap.Data()), "maps", 0.32, rows);
  }

  // ======================================================================
  // 10. linker diagnostics
  // ======================================================================
  printf("\n=== 10. linker: linked fraction of the generated pT, per class and region ===\n");
  TexSection("Linker diagnostics", "Everything above depends on the gen-particle to PF linker. Its diagnostics: the linked "
             "fraction of the generated $\\pt$ per particle class and detector region; the parallax slope, "
             "$\\Delta\\eta$ of the PF photon against $(z_{\\mathrm{gen}}-z_{\\mathrm{PV}})/(r_{\\mathrm{ECAL}}\\cosh\\eta)$, "
             "of magnitude 1 if PF points from the primary vertex and the particle left its own one (v4 measured $+0.98$); "
             "the bend slope of $q\\,\\Delta\\phi$ against $1/\\pt$, $-0.735$\\,rad\\,GeV at the ECAL front face for $B=3.8$\\,T "
             "(v4 measured $-0.97$: the deposit sits deeper than the front face); the $\\Delta z$ of the track links; the "
             "single-particle response $\\pt^{\\mathrm{PF}}/G$ per PF class for single-partner candidates; and the orphan "
             "fraction of the PF $\\pt$ against $N_{\\mathrm{PU}}$.");
  {
    std::vector<TString> figs; TString cap;
    // hlinkfrac is filled with weight pT: its entries are not particles, and
    // a cell in which every particle was linked has error 0 - so a cell is
    // asked for its effective entries (kMinEffLink), never for its error
    const double kMinEffLink = 10;
    for (int r = 0; r < kNReg; ++r) {
      TH1D *fr[kNCls] = {0}; TProfile *pr[kNCls] = {0};
      for (int c = 0; c < kNCls; ++c)
        if ((pr[c] = GetP(Form("link/hlinkfrac_%s_%s", ClsTag(c), RegTag(r))))) fr[c] = ProfX(pr[c], Form("lf_%d_%d", c, r));
      TH1D *any = 0; for (int c = 0; c < kNCls; ++c) if (fr[c]) { any = fr[c]; break; }
      if (!any) continue;
      // opened here, after the check, so that it is always closed below
      if (r == 0) TexTableBegin(Form("Linked fraction of the generated $\\pt$ per particle class at $|\\eta|<2.5$, in a few bins "
                                     "of the particle $\\pt$ (the profile is $\\pt$-weighted; a cell needs %.0f effective "
                                     "particles, else --).", kMinEffLink),
                                "linkfrac", "lrrrrrr", "$\\pt^{\\mathrm{gen}}$ [GeV] & $\\pi^\\pm$ & $K^\\pm$, p & $\\gamma$ & n, $K^0_L$ & $V^0$ & e, $\\mu$");
      const double xlo = std::max(AxLo(any), 0.05), xhi = any->GetXaxis()->GetXmax();
      TH1D *frame = tdrHist(Form("flf%d", r), "Linked fraction of generated p_{T}", 0., 1.39, "p_{T}^{gen particle} [GeV]", xlo, xhi);
      std::unique_ptr<TCanvas> c(tdrCanvas(Form("c_lf%d", r), frame, 8, 11, kSquare));
      c->SetLogx();
      UnitLine(xlo, xhi);
      TLegend *leg = tdrLeg(0.55, 0.17, 0.94, 0.47);
      leg->SetTextSize(0.030); leg->SetHeader(RegLabel(r));
      for (int c2 = 0; c2 < kNCls; ++c2) if (fr[c2]) {
        tdrDraw(fr[c2], "Pz", ClsMarker[c2], ClsColour[c2], kSolid, -1, kNone, 0, 0.9);
        leg->AddEntry(fr[c2], ClsLabel(c2), "PL");
      }
      Thresholds();
      fixOverlay();
      figs.push_back(SavePdf(c.get(), Form("linkfrac_%s", RegTag(r))));
      printf("  %-12s", RegTag(r));
      for (int c2 = 0; c2 < kNCls; ++c2) printf("%9s", ClsTag(c2));
      printf("\n");
      // one row per bin of the profile axis that holds a probe value, labelled
      // by the bin it is (the average is over the bin, not at the probe)
      int last = -1;
      for (double pt : {0.3, 0.5, 1., 2., 5., 10., 20.}) {
        const int b = any->FindFixBin(pt);
        if (b < 1 || b > any->GetNbinsX() || b == last) continue;
        last = b;
        bool anyCell = false;                        // a row of "--" says nothing
        for (int c2 = 0; c2 < kNCls; ++c2) anyCell = anyCell || HasP(pr[c2], b, kMinEffLink);
        if (!anyCell) continue;
        printf("  %-12s", BinLabel(any->GetXaxis(), b).Data());
        if (r == 0) gTexT << BinTex(any->GetXaxis(), b);
        for (int c2 = 0; c2 < kNCls; ++c2) {
          printf("%9s", ValP(pr[c2], b, "%.3f", kMinEffLink, "-").Data());
          if (r == 0) gTexT << " & " << ValP(pr[c2], b, "%.3f", kMinEffLink);
        }
        printf("\n");
        if (r == 0) gTexT << " \\\\\n";
      }
      if (r == 0) { TexTableEnd();
        const int b = any->FindFixBin(1.);
        cap = Form("At $|\\eta|<2.5$ and %s\\,GeV the linked fraction is %s for $\\pi^\\pm$ and %s for $\\gamma$.",
                   BinTex(any->GetXaxis(), b).Data(), ValP(pr[0], b, "%.2f", kMinEffLink).Data(),
                   ValP(pr[2], b, "%.2f", kMinEffLink).Data()); }
    }
    TexFigure(figs, "Linked fraction of the generated $\\pt$ per particle class against the particle $\\pt$, at $|\\eta|<2.5$, "
              "$2.5<|\\eta|<3.0$ and $|\\eta|>3.0$. " + cap, "linkfrac");

    // parallax and bend slopes: the two geometric checks of the calo targets
    figs.clear(); cap = "";
    struct SlopeCheck { const char *hist, *xlab, *ylab, *pdf; double expect; };
    const SlopeCheck slopes[2] = {
      {"link/hparallax", "(z_{gen} - z_{PV}) / (r_{ECAL} cosh#eta_{gen})", "#eta_{PF} - #eta_{gen}", "parallax", 1.00},
      {"link/hbend", "1 / p_{T}^{gen} [GeV^{-1}]", "q (#phi_{PF} - #phi_{gen})", "bend", -0.735}};
    for (const SlopeCheck &sc : slopes) {
      TH2D *h2 = Get2(sc.hist);
      if (!h2 || h2->GetEntries() < 10) continue;
      const TAxis *ax = h2->GetXaxis(), *ay = h2->GetYaxis();
      TH1D *frame = tdrHist(Form("f_%s", sc.pdf), sc.ylab, ay->GetXmin(), ay->GetXmax(), sc.xlab, ax->GetXmin(), ax->GetXmax());
      std::unique_ptr<TCanvas> c(tdrCanvas(Form("c_%s", sc.pdf), frame, 8, 11, kSquare));
      h2->Draw("COL SAME");
      TH1D *pp = MeanY(h2, Form("pp_%s", sc.pdf));
      TF1 f1(Form("f1_%s", sc.pdf), "pol1", ax->GetXmin(), ax->GetXmax());
      TF1 f0(Form("f0_%s", sc.pdf), "[0]*x", ax->GetXmin(), ax->GetXmax());
      const int st1 = pp->Fit(&f1, "QRN0"), st0 = pp->Fit(&f0, "QRN0");
      tdrDraw(pp, "Pz", kFullCircle, kBlack, kSolid, -1, kNone, 0, 0.7);
      f1.SetLineColor(kRed+1); f1.SetLineWidth(2); f1.Draw("SAME");
      printf("\n=== %s slope: pol1 %.4f +- %.4f (offset %.4f +- %.4f, chi2/ndf %.1f/%d), through origin %.4f +- %.4f; expect %s%.3f%s\n",
             sc.pdf, f1.GetParameter(1), f1.GetParError(1), f1.GetParameter(0), f1.GetParError(0), f1.GetChisquare(), f1.GetNDF(),
             f0.GetParameter(0), f0.GetParError(0), sc.expect > 0 ? "|slope| = " : "", sc.expect, (st1 || st0) ? "  (a fit FAILED)" : "");
      Note(0.93, 0.83, Form("slope %.3f #pm %.3f  (expect %s%.3f)", f1.GetParameter(1), f1.GetParError(1), sc.expect > 0 ? "#pm" : "", sc.expect), 0.030, 31);
      Note(0.93, 0.79, Form("offset %.4f #pm %.4f", f1.GetParameter(0), f1.GetParError(0)), 0.030, 31);
      Thresholds();
      fixOverlay();
      figs.push_back(SavePdf(c.get(), sc.pdf));
      cap += Form("The %s slope is $%.3f\\pm%.3f$ (expected $%s%.3f$). ", sc.pdf, f1.GetParameter(1), f1.GetParError(1), sc.expect > 0 ? "\\pm" : "", sc.expect);
    }
    // dz of the track links, per region
    {
      TH1D *u[kNReg] = {0}; double xlo = -1, xhi = 1, ymax = 0;
      for (int r = 0; r < kNReg; ++r) {
        TH1D *h = Get1(Form("link/hlinkdz_%s", RegTag(r)));
        if (!h || h->Integral() <= 0) continue;
        u[r] = Unit(h, Form("udz%d", r));
        xlo = u[r]->GetXaxis()->GetXmin(); xhi = u[r]->GetXaxis()->GetXmax();
        ymax = std::max(ymax, u[r]->GetMaximum());
      }
      if (ymax > 0) {
        TH1D *frame = tdrHist("fdz", "Fraction of track links", ymax*1e-4, ymax*30, "z_{PF} - z_{gen} [cm]", xlo, xhi);
        std::unique_ptr<TCanvas> c(tdrCanvas("c_dz", frame, 8, 11, kSquare));
        c->SetLogy();
        TLegend *leg = tdrLeg(0.62, 0.66, 0.94, 0.84);
        leg->SetTextSize(0.028);
        printf("\n=== track-link dz: rms and core sigma [cm] per region ===\n");
        int nn = 0;
        for (int r = 0; r < kNReg; ++r) if (u[r]) {
          double mu = u[r]->GetMean(), sig = u[r]->GetRMS();
          TF1 fg(Form("fgdz%d", r), "gaus", xlo, xhi);
          for (int it = 0; it < 5 && sig > 0; ++it) {
            fg.SetRange(mu - 2*sig, mu + 2*sig); fg.SetParameters(u[r]->GetMaximum(), mu, sig);
            if (u[r]->Fit(&fg, "QRN0") != 0) break;
            mu = fg.GetParameter(1); sig = fabs(fg.GetParameter(2));
          }
          printf("  %-6s rms %.4f  core %.4f (mu %.4f)\n", RegTag(r), u[r]->GetRMS(), sig, mu);
          tdrDraw(u[r], "HIST", kNone, ClsColour[r], kSolid, -1, kNone, 0, 0, 2);
          leg->AddEntry(u[r], RegLabel(r), "L");
          Note(0.19, 0.78 - 0.035*nn++, Form("#sigma_{core} = %.3f cm", sig), 0.026, 11, ClsColour[r]);
          if (r == 0) cap += Form("The core of $\\Delta z$ of the track links is %.3f\\,cm at $|\\eta|<2.5$. ", sig);
        }
        Thresholds();
        fixOverlay();
        figs.push_back(SavePdf(c.get(), "linkdz"));
      }
    }
    // single-particle response pT_PF / G per PF class
    {
      TH1D *u[kNPF] = {0}; double ymax = 0, xlo = 0, xhi = 3;
      printf("\n=== pT_PF / G for single-partner PF: median, peak and core width per PF class ===\n");
      for (int c = 0; c < kNPF; ++c) {
        TH1D *h = Get1(Form("link/hpfratio_%s", PFTag(c)));
        if (!h || h->Integral() < kMinN) continue;
        u[c] = Unit(h, Form("upf%d", c));
        xlo = u[c]->GetXaxis()->GetXmin(); xhi = u[c]->GetXaxis()->GetXmax();
        ymax = std::max(ymax, u[c]->GetMaximum());
        double mu = h->GetBinCenter(h->GetMaximumBin()), sig = 0.15;
        TF1 fg(Form("fgpf%d", c), "gaus", xlo, xhi);
        for (int it = 0; it < 5; ++it) {
          fg.SetRange(mu - 1.5*sig, mu + 1.5*sig); fg.SetParameters(h->GetMaximum(), mu, sig);
          if (h->Fit(&fg, "QRN0") != 0) break;
          mu = fg.GetParameter(1); sig = fabs(fg.GetParameter(2));
        }
        double q, pr = 0.5; h->GetQuantiles(1, &q, &pr);
        printf("  %-6s n %8.0f  median %.3f  peak %.3f  sigma %.3f\n", PFTag(c), h->Integral(0, h->GetNbinsX()+1), q, mu, sig);
        if (c == 0) cap += Form("The charged-hadron $\\pt^{\\mathrm{PF}}/G$ peaks at %.3f with a core width of %.3f. ", mu, sig);
      }
      if (ymax > 0) {
        TH1D *frame = tdrHist("fpf", "Fraction of PF candidates", 0, ymax*1.6, "p_{T}^{PF} / G_{T} (routed gen p_{T})", xlo, xhi);
        std::unique_ptr<TCanvas> c(tdrCanvas("c_pf", frame, 8, 11, kSquare));
        TLegend *leg = tdrLeg(0.60, 0.58, 0.94, 0.86);
        leg->SetTextSize(0.028);
        for (int c2 = 0; c2 < kNPF; ++c2) if (u[c2]) {
          tdrDraw(u[c2], "HIST", kNone, ClsColour[c2], kSolid, -1, kNone, 0, 0, 2);
          leg->AddEntry(u[c2], PFLabel(c2), "L");
        }
        TLine l; l.SetLineStyle(kDashed); l.SetLineColor(kGray+1); l.DrawLine(1, 0, 1, ymax*1.6);
        Thresholds();
        fixOverlay();
        figs.push_back(SavePdf(c.get(), "pfratio"));
      }
    }
    // orphan fraction against N_PU, per PF class, one canvas per region
    printf("\n=== orphan pT fraction at N_PU = 30 / 50 / 70, per PF class and region ===\n%-6s", "");
    for (int c = 0; c < kNPF; ++c) printf("%18s", PFTag(c));
    printf("\n");
    for (int r = 0; r < kNReg; ++r) {
      TH1D *o[kNPF] = {0}; TH1D *any = 0;
      printf("%-6s", RegTag(r));
      for (int c = 0; c < kNPF; ++c) {
        TProfile *p = GetP(Form("link/horphan_%s_%s", PFTag(c), RegTag(r)));
        if (!p) { printf("%18s", "-"); continue; }
        o[c] = ProfX(p, Form("orph%d%d", c, r)); any = o[c];
        printf("   %.3f/%.3f/%.3f", o[c]->GetBinContent(o[c]->FindFixBin(30)), o[c]->GetBinContent(o[c]->FindFixBin(50)), o[c]->GetBinContent(o[c]->FindFixBin(70)));
        if (r == 0 && c == 0) cap += Form("The charged-hadron orphan fraction at $|\\eta|<2.5$ is %.3f at $N_{\\mathrm{PU}}=50$.", o[c]->GetBinContent(o[c]->FindFixBin(50)));
      }
      printf("\n");
      if (!any) continue;
      TH1D *frame = tdrHist(Form("forph%d", r), "Orphan fraction of the PF p_{T}", 0., 0.59, "N_{PU}", any->GetXaxis()->GetXmin(), any->GetXaxis()->GetXmax());
      std::unique_ptr<TCanvas> c(tdrCanvas(Form("c_orph%d", r), frame, 8, 11, kSquare));
      TLegend *leg = tdrLeg(0.55, 0.58, 0.94, 0.86);
      leg->SetTextSize(0.028); leg->SetHeader(RegLabel(r));
      for (int c2 = 0; c2 < kNPF; ++c2) if (o[c2]) {
        tdrDraw(o[c2], "Pz", ClsMarker[c2], ClsColour[c2], kSolid, -1, kNone, 0, 0.9);
        leg->AddEntry(o[c2], PFLabel(c2), "PL");
      }
      Thresholds();
      fixOverlay();
      figs.push_back(SavePdf(c.get(), Form("orphan_%s", RegTag(r))));
    }
    // two rows as the caption says them: the geometry, then the responses
    int nGeo = 0;
    for (const TString &f : figs) if (f == "parallax" || f == "bend" || f == "linkdz") ++nGeo;
    TexFigure(figs, "Linker geometry and response. First row: the parallax check, $\\eta_{\\mathrm{PF}}-\\eta_{\\mathrm{gen}}$ of "
              "single-source photon links against $(z_{\\mathrm{gen}}-z_{\\mathrm{PV}})/(r_{\\mathrm{ECAL}}\\cosh\\eta)$, and the "
              "bend check, $q\\,(\\phi_{\\mathrm{PF}}-\\phi_{\\mathrm{gen}})$ of untracked charged links against $1/\\pt$, each "
              "with its profile and a linear fit; the $\\Delta z$ of the track links per region. Second row: the "
              "single-particle response $\\pt^{\\mathrm{PF}}/G$ per PF class, and the orphan fraction of the PF $\\pt$ against "
              "$N_{\\mathrm{PU}}$ per PF class in the three regions. " + cap, "linker", 0.24, {nGeo, (int)figs.size() - nGeo});
  }

  gTexP.close(); gTexT.close();
  printf("\n=== core-fit protocol: %d distributions with n >= %.0f, %d converged, %d failed as sigma_core > 1.1 rms, "
         "%d as mu_core outside [0.5, 1.3] median\n", gNDesc, kMinN, gNConv, gNWide, gNMu);
  printf("\ndrawGenSeed: wrote %d PDFs to plots/ (%d cropped), doc/plots.tex, doc/tables.tex and text/jec%s.txt; %d histogram(s) missing\n",
         gNPdf, gNCrop, tag, gMissing);
}
