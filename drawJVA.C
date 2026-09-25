// drawJVA.C - every plot, table and generated LaTeX of the JVA study (jet-vertex
// association of the forward jets by global MET minimisation), from the
// histogram file that JVA.C writes.
//
//   root -l -b -q -e 'gSystem->SetBuildDir("'$TMPDIR'/aclic_drawjva",true);' \
//        'drawJVA.C+("rootfiles/JVA_v1.root","_v1")'
//
// An optional third argument, jerMerge (default 2), is the number of fine
// pT_avg bins merged into one point of the JER (section 5; 1 = the fine axis).
//
// Writes plots/jva_*.pdf (cropped), doc/jva_plots.tex and doc/jva_tables.tex
// (\input by doc/jva.tex; \plotdir and \pt are \providecommand'ed), and prints
// every number that goes into a table.  The tag names the production in the
// generated files; the output names are fixed.  Build it outside the source
// tree (the SetBuildDir above): ACLiC otherwise drops its .so/.d/.pcm next to
// the macro, in a synced directory.
//
// WHY THIS EXISTS.  Per-vertex PUPPI cannot tell the vertices apart beyond the
// tracker: every candidate at |eta| > 2.5 has the same weight at every vertex,
// so every vertex carries the same forward jets.  Those jets are mostly
// combinations of several interactions (the 2.5 < |eta| < 3.0 jet spike), they
// put ~20 GeV of MET on every vertex, and a dijet whose probe is one of them
// measures nothing.  JVA hands each forward cluster to the one vertex whose
// recoil it balances best - a global minimisation of the summed MET chi2 over
// all vertices, with a price per assignment and the in-cone tracks narrowing the
// candidates - or to nobody.  Seven ways of distributing the vertex-blind energy
// are compared on the same crossings: dup (every vertex, PUPPI as it is; at the
// PV this is CMS PUPPI), lv (the PV only), trk (the candidate vertex with the
// most track evidence), jva (this work), two oracles from the generator links,
// ojet (whole clusters to the vertex of their dominant interaction) and opart
// (every candidate to its own interaction: the upper bound), and none (nobody:
// every forward cluster dropped, the reference any association has to beat;
// drawn black and dashed wherever it appears).  At the PV, dup and
// lv are the same event, the standard leading-vertex reconstruction; the
// question this macro answers is how much of the other vertices' events becomes
// usable, and down to which pT.
//
// WHAT IS MEASURED.  (1) The MET of every vertex: the width of its components
// per method and vertex class, against N_PU and against the scalar sum S_v of
// the vertex-resolved part, and the residual against the true visible
// imbalance of the owner interaction; the residual width against S_v is the
// noise model the optimiser assumes, sigma_v^2 = sigma0^2 + sigmaK^2 S_v, and
// that curve is drawn on top of it.  (2) The forward spike: jets per owned
// vertex against |eta| over the pure GenJets of the owner interaction, the pT
// spectra in three regions, the fraction of a jet's linked generated pT that
// is the owner's, and the dominant-interaction share.  (3) The assignment
// itself, cluster by cluster: right vertex, wrong vertex, single nulled,
// combination nulled (the right thing to do with a combination) and combination
// kept, with and without track candidates; "correct" is right + combination
// nulled, and the table carries the reference of dropping every cluster,
// which scores the combination fraction: the none method's own outcomes
// (fwd/none), or for a file without them the combinations of the first method
// present (single or combination is a property of the cluster, not of the
// method).  (4) L2Res from every vertex: R_DB and R_MPF against
// probe |eta| at five pT_avg bins, closed against r_true, the median of the
// per-event truth relative response of the same dijets; with alpha < 0.3 and,
// as a second set of figures (_a1), without the alpha cut, which is the one
// with points on a small sample.  (5) The dijet JER by the DB-bisector
// and MPFX methods, JER^2 = (RMS_par^2 - RMS_perp^2)/2, against the truth core
// width of pT^corr/pT^linkz, and the lowest pT_avg that each method reaches
// within 10%.  (6) Summary tables.
//
// TESTED ON A TOY with a known answer (a synthetic file of the dijet
// histograms only: Gaussian responses of known width, an isotropic soft
// imbalance, a probe scale c(|eta|) = 1 - 0.02 |eta|): R_DB, R_MPF and the
// median r_true all return c to <1% at 8-30 GeV; the bisector JER from the
// histograms equals the directly computed sqrt((<par^2> - <perp^2>)/2) of
// the same dijets, for the +-v fill of JVA.C and for a pT-ordered fill alike.
//
// THE PROTOCOLS.  Two, and only two.  A RESPONSE (pT^corr/pT^linkz, around 1
// with a right tail) gets drawGenSeed.C's core-fit protocol unchanged: the
// median first, then a Gaussian iterated to convergence in the asymmetric
// window [mu - 2 sigma, mu + sigma], failed (never replaced by the rms) below
// 200 entries, without convergence, when sigma_core > 1.1 rms or when mu_core
// leaves [0.5, 1.3] x median.  The one change: the likelihood needs counts, so
// a histogram filled with one common weight w (JVA.C's 1/2) is fitted as
// counts, content/w, and one of mixed weights by chi2 (ROOT 6.28's weighted
// likelihood crashes with a fixed parameter).  A ZERO-CENTRED distribution (a
// MET component, a MET residual, a bisector projection) has no median to scale
// by and no preferred tail, so its Gaussian is iterated in the symmetric
// window mu +- 2 sigma, the same way, and
// fails below 200 entries, without convergence, when |mu_core - median| >
// 0.5 sigma_core, or when the core is wider than the distribution - where "the
// distribution" is what a Gaussian of that width would show inside the axis
// (its truncated rms must not exceed 1.1 x the observed rms).  Without that
// qualification the bisector projections at pT_avg ~ 4 GeV, which fill the
// whole -1..1 axis, would fail exactly where a Gaussian fit is the only width
// estimate that the axis does not bias.  The plain rms is always the in-axis
// rms, and a truncated axis biases it low: that is what the core is for.
//
// THE BISECTOR PROJECTIONS ARE SYMMETRISED.  n = (u1 - u2)/|u1 - u2| changes
// sign when the two jets are swapped, and so do both projections.  If the
// analyzer takes the jets in pT order, (p1 + p2).n = (pT1 - pT2) cos(theta) >= 0
// and the DB-par distribution is folded at zero (and MET.n is mostly negative);
// its rms about its mean is not a width.  The mirror sum h(x) + h(-x) is the
// distribution under a random order, whatever order was used, and its rms
// (about zero) and its core are the widths; for a random order it only halves
// the noise.  JVA.C already fills every projection at +v and -v with weight
// 1/2: such a histogram is recognised as symmetric and used as it is.  Either
// way every dijet is in the fitted histogram twice, so the 200-entry threshold
// and the errors count the independent dijets (the fit errors are scaled by
// sqrt(neff/n_dijets) = sqrt(2)).  An axis 0..x is taken as already folded
// (|v|) and mirrored.  The JER per jet is then sigma = sqrt(P^2 - Q^2)/sqrt(2),
// P and Q the widths of par and perp: P^2 holds both jets' resolution plus the
// imbalance along the axis, Q^2 the imbalance across it.  A bin with Q >= P
// has no JER.
//
// THE TRUTH.  JER_truth is sigma_core/mu_core of pT^corr/pT^linkz of both jets
// of the same same-bin dijets (hrlink), in the same pT_avg bin: the bisector
// projections are normalised to pT_avg in the corrected scale, so the truth is
// the relative width in that scale.  That is the truth of the core-fit JER.
// The plain-rms JER is compared with rms/mean of the same hrlink, like with
// like: on a Gaussian toy the plain-rms JER equals the truth rms/mean to 1%
// above 12 GeV, and the truth core is 5-10% narrower than its rms because a
// pT_avg bin mixes jets of different true pT.  A bin where both widths exist
// but RMS_perp >= RMS_par has no JER; it is a FAILED bin, drawn as a cross on
// the bottom edge of the ratio panel and counted.  The REACH of a method is
// the lower edge of the lowest pT_avg bin whose ratio JER/JER_truth is within
// 10% of 1 and whose next two measured bins (fewer where the measurement ends)
// are not significantly outside the band (|ratio - 1| < 0.1 + 2 sigma, and not
// failed): one lucky low bin is not a reach, and neither is a range with a
// failed bin in it.  The JER needs 200 same-bin dijets per pT_avg bin, and
// the pT_avg axis of section 5 merges two fine bins per point from 3 GeV
// (3-5, 5-8, 8-12, ...; the third argument).  L2Res closure is the mean of
// |R/r_true - 1| over the probe |eta| bins where both are measured (at least
// 50 effective dijets), with the number of such bins; R = (1 + <A>)/(1 - <A>).
//
// ROBUSTNESS.  A missing (or mistyped) histogram is reported once by name -
// the first 25, then per directory at the end - and the curve, or the
// section, it belongs to is skipped; a number that does not exist is "--" in
// every table and "not defined" in every caption, never 0.  A figure lists only
// the PDFs that were saved, a table is written only if it got a row, a section
// only if it produced a figure, and a reference to a figure or table that was
// not produced says so in words (no "??").  Axes are cloned from the histograms
// (Like()), never assumed; names are tried as given and under hist/.  The probe
// |eta| edges, the method, class and outcome names are the contract with JVA.C
// and are written here once.
//
// genseed - Copyright 2026 Mikko Voutilainen - Apache License 2.0
#include "tdrstyle_mod22.C"

#include <TFile.h>
#include <TH1D.h>
#include <TH2.h>
#include <TProfile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TSystem.h>
#include <TError.h>
#include <TMath.h>
#include <TNamed.h>
#include <fstream>
#include <memory>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <cstring>
#include <algorithm>

namespace {
  // -------------------------------------------------------------------
  // the contract with JVA.C: methods, vertex classes, outcomes, bins
  const int kNM = 7;
  enum { kDup = 0, kLv, kTrk, kJva, kOJet, kOPart, kNobody };
  const char* MTag(int m)   { static const char *s[kNM] = {"dup","lv","trk","jva","ojet","opart","none"}; return s[m]; }
  const char* MShort(int m) { static const char *s[kNM] = {"dup","lv","trk","jva","oracle jet","oracle part.","none (drop)"}; return s[m]; }
  const char* MTex(int m)   { static const char *s[kNM] = {"\\texttt{dup}","\\texttt{lv}","\\texttt{trk}","\\texttt{jva}",
                                                           "\\texttt{ojet}","\\texttt{opart}","\\texttt{none}"}; return s[m]; }
  // the oracles are open markers and broken lines, so they read as references;
  // so is none, the drop-all reference, black and dashed
  const int MColour[kNM] = {kGray+2, kBlue+1, kGreen+2, kRed+1, kOrange+2, kMagenta+1, kBlack};
  const int MMarker[kNM] = {kFullCircle, kFullSquare, kFullTriangleDown, kFullTriangleUp, kOpenDiamond, kOpenCircle, kOpenSquare};
  const int MLine[kNM]   = {kSolid, kSolid, kSolid, kSolid, kDashed, kDotted, kDashed};
  const char* MColTex(int m) { static const char *s[kNM] = {"grey","blue","green","red","open orange","open magenta","open black, dashed"}; return s[m]; }
  const double kMidM = 0.5*(kNM - 1);     // the middle of the methods' side-by-side offsets
  TString MWho(const bool *have) {        // "dup (grey), lv (blue) and jva (red)"
    TString who;
    for (int m = 0; m < kNM; ++m) if (have[m]) who += Form("%s%s (%s)", who.IsNull() ? "" : ", ", MTex(m), MColTex(m));
    const Ssiz_t k = who.Last(','); if (k != kNPOS) who.Replace(k, 1, " and");
    return who;
  }

  // vertex classes: MET in all / PV (key 0) / the others; dijets in all / PV
  const int kNC = 3, kNDC = 2;
  const char* CTag(int c)   { static const char *s[kNC] = {"all","pv","pu"}; return s[c]; }
  const char* CLabel(int c) { static const char *s[kNC] = {"all vertices","primary vertex","other vertices"}; return s[c]; }
  const char* CTex(int c)   { static const char *s[kNC] = {"all vertices","PV","other vertices"}; return s[c]; }

  // the four methods that assign whole clusters (fwd/<m>/hcat_*) and the outcomes
  const int kNA = 4;
  const int AMethod[kNA] = {kLv, kTrk, kJva, kOJet};
  const int kNOut = 5;
  enum { kRight = 0, kWrong, kNullS, kNullC, kKeptC };
  const char* OutTag(int o)   { static const char *s[kNOut] = {"right","wrong","nulled_single","nulled_combo","kept_combo"}; return s[o]; }
  const char* OutLabel(int o) { static const char *s[kNOut] = {"single: right vertex","single: wrong vertex","single: nulled",
                                                               "combination: nulled","combination: kept"}; return s[o]; }
  const int OutColour[kNOut] = {kGreen+2, kRed+1, kOrange+2, kBlue+1, kMagenta+1};
  const int OutMarker[kNOut] = {kFullCircle, kFullSquare, kFullTriangleDown, kFullTriangleUp, kFullDiamond};
  const char* SelTag(int t)   { return t == 0 ? "trk" : "notrk"; }
  const char* SelLabel(int t) { return t == 0 ? "with track candidates" : "without track candidates"; }

  // the probe |eta| bins of CMS L2Res, e00..e17; the tag region
  const int kNP = 18;
  const double kProbeEta[kNP+1] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.930, 2.172,
                                   2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};
  TString PTag(int i)   { return Form("e%02d", i); }
  TString PTex(int i)   { return Form("$%g<|\\eta|<%g$", kProbeEta[i], kProbeEta[i+1]); }
  // The JER bins: the 18 same-bin |eta| bins of the contract, then four coarse
  // regions made of whole fine bins (a sum of same-bin dijets is still a
  // same-bin sample): barrel, endcap, the transition and HF.  The fine bins
  // are thin - JVA.C's 300-crossing test has 0.08 same-bin dijets per crossing
  // over all 18 of them - and the regions are what the summary quotes.
  const int kNJR = 4;
  const int JRLo[kNJR] = {0, 5, 11, 15}, JRHi[kNJR] = {4, 10, 14, 17};
  const int kNJB = kNP + kNJR;
  int JBLo(int k) { return k < kNP ? k : JRLo[k - kNP]; }
  int JBHi(int k) { return k < kNP ? k : JRHi[k - kNP]; }
  TString JBTag(int k)   { return k < kNP ? PTag(k) : TString(Form("r%d", k - kNP)); }
  TString JBLabel(int k) { return Form("%g < |#eta| < %g", kProbeEta[JBLo(k)], kProbeEta[JBHi(k)+1]); }
  TString JBTex(int k)   { return Form("$%g<|\\eta|<%g$", kProbeEta[JBLo(k)], kProbeEta[JBHi(k)+1]); }
  TString JBTxt(int k)   { return Form("%g-%g", kProbeEta[JBLo(k)], kProbeEta[JBHi(k)+1]); }

  // jet regions (jets/<m>/hjetpt_r<r>, hdom_r<r>) and the thresholds of hjeteta
  const int kNReg = 3;
  const char* RegLabel(int r) { static const char *s[kNReg] = {"|#eta| < 2.5","2.5 < |#eta| < 3.0","3.0 < |#eta| < 5.0"}; return s[r]; }
  const int kNThr = 3;
  const int kThr[kNThr] = {5, 10, 20};

  // dijets: alpha cuts a0 (alpha < 0.3), a1 (none); the L2Res pT_avg bins
  const int kNAl = 2;
  const char* AlLabel(int a) { return a == 0 ? "#alpha < 0.3" : "no #alpha cut"; }
  const int kNL2 = 5;
  const double kL2Pt[kNL2] = {5, 8, 12, 18, 28};
  enum { kDB = 0, kMPF, kNRV };
  const char* RVTag(int v)   { return v == kDB ? "db" : "mpf"; }
  const char* RVLabel(int v) { return v == kDB ? "R_{DB}" : "R_{MPF}"; }
  const char* RVTex(int v)   { return v == kDB ? "$R_{\\mathrm{DB}}$" : "$R_{\\mathrm{MPF}}$"; }
  // JER estimators: DB-bisector and MPFX; width variants: core fit, plain rms
  enum { kJDB = 0, kJMPFX, kNJ };
  const char* JTag(int j)   { return j == kJDB ? "db" : "mpfx"; }
  const char* JLabel(int j) { return j == kJDB ? "DB bisector" : "MPFX"; }
  const char* JHist(int j, bool par) { return j == kJDB ? (par ? "hdbpar" : "hdbperp") : (par ? "hmpfpar" : "hmpfperp"); }

  const double kMinN = 200.;      // entries a distribution needs to be described
  const double kMinDen = 20.;     // denominator a fraction needs to be drawn
  const double kMinMean = 50.;    // effective entries a mean or a median needs
  const double kPtDrawHi = 100.;  // upper edge of the pT axes of the metric plots
  const double kJerTol = 0.10;    // |JER/JER_truth - 1| of the reach
  const double kSigmaInel_mb = 80.0;  // 13.6 TeV minimum bias, the pileup library's value

  // -------------------------------------------------------------------
  // getters: a name is tried as given and under hist/; a missing histogram
  // (or one of the wrong type) is reported once by name and counted - the
  // first kMaxReport by name, the rest per directory at the end (a file from
  // an older analyzer, or an empty one, would otherwise print thousands)
  TFile *gF = 0;
  int    gMissing = 0;
  const int kMaxReport = 25;
  std::set<std::string> gReported;
  std::map<std::string, int> gMissingDir;
  void ReportMissing(const char *name, const char *what) {
    ++gMissing;
    if (gMissing <= kMaxReport) printf("drawJVA: %s %s\n", name, what);
    else if (gMissing == kMaxReport + 1) printf("drawJVA: ... (further missing histograms are summarised per directory at the end)\n");
    const char *sl = strrchr(name, '/');
    ++gMissingDir[sl ? std::string(name, sl - name) : std::string("(top)")];
  }
  TH1* GetH(const char *name, bool quiet = false) {
    TH1 *h = gF ? dynamic_cast<TH1*>(gF->Get(name)) : 0;
    if (!h && gF) h = dynamic_cast<TH1*>(gF->Get(Form("hist/%s", name)));
    if (!h) {
      if (!quiet && gReported.insert(name).second) ReportMissing(name, "is missing");
      return 0;
    }
    h->SetDirectory(0);
    return h;
  }
  template<class T> T* GetT(const char *name, bool quiet = false) {
    TH1 *h = GetH(name, quiet);
    if (!h) return 0;
    T *t = dynamic_cast<T*>(h);
    if (!t && gReported.insert(std::string(name) + "#type").second) ReportMissing(name, Form("is a %s, not a %s", h->ClassName(), T::Class_Name()));
    return t;
  }
  TH2*      Get2(const char *n, bool quiet = false) { return GetT<TH2>(n, quiet); }
  // the sum of <base><hist>_eXX over the fine |eta| bins [lo, hi] (one bin: the histogram itself)
  TH2* SumEta2(const TString &base, const char *hist, int lo, int hi, const char *name) {
    TH2 *sum = 0;
    for (int ie = lo; ie <= hi; ++ie) {
      TH2 *h = Get2(base + hist + "_" + PTag(ie));
      if (!h) continue;
      if (!sum) { sum = (TH2*)h->Clone(name); sum->SetDirectory(0); } else sum->Add(h);
    }
    return sum;
  }
  TProfile* GetP(const char *n, bool quiet = false) { return GetT<TProfile>(n, quiet); }
  // the first of several spellings that exists; the first is reported if none
  TH1* GetAny(const std::vector<TString> &names) {
    for (const TString &n : names) if (TH1 *h = GetH(n, true)) return h;
    GetH(names[0]);
    return 0;
  }
  TString GetText(const char *name) {
    TNamed *o = gF ? dynamic_cast<TNamed*>(gF->Get(name)) : 0;
    if (!o && gF) o = dynamic_cast<TNamed*>(gF->Get(Form("hist/%s", name)));
    return o ? TString(o->GetTitle()) : TString("");
  }

  // empty TH1D with the x axis of src: the ONLY way an axis is made here
  TH1D* Like(const TH1 *src, const char *name) {
    const TAxis *a = src->GetXaxis();
    TH1D *h = a->GetXbins()->GetSize() > 0
      ? new TH1D(name, "", a->GetNbins(), a->GetXbins()->GetArray())
      : new TH1D(name, "", a->GetNbins(), a->GetXmin(), a->GetXmax());
    h->SetDirectory(0); h->Sumw2();
    return h;
  }
  // a TH1D copy of any 1D histogram (TH1F, TH1I, ...): the protocols take TH1D
  TH1D* AsD(TH1 *h, const char *name) {
    if (!h) return 0;
    if (TH1D *d = dynamic_cast<TH1D*>(h)) if (!dynamic_cast<TProfile*>(h)) return d;
    TH1D *d = Like(h, name);
    for (int b = 0; b <= h->GetNbinsX()+1; ++b) { d->SetBinContent(b, h->GetBinContent(b)); d->SetBinError(b, h->GetBinError(b)); }
    d->ResetStats();
    return d;
  }
  // where a falling spectrum puts the mean of a bin; the middle of a bin that starts at 0
  double GeoCentre(const TAxis *a, int b) {
    const double lo = a->GetBinLowEdge(b), hi = a->GetBinUpEdge(b);
    return lo > 0 ? sqrt(lo*hi) : 0.5*(lo + hi); }
  TString BinLabel(const TAxis *a, int b) { return Form("%g-%g", a->GetBinLowEdge(b), a->GetBinUpEdge(b)); }
  TString BinTex(const TAxis *a, int b)   { return Form("%g--%g", a->GetBinLowEdge(b), a->GetBinUpEdge(b)); }
  int BinAt(const TH1 *h, double x) { return h->GetXaxis()->FindFixBin(x + 1e-3); }
  double AxLo(const TH1 *h) { return h->GetXaxis()->GetXmin(); }

  // fractions: error floored at 1/d, so that error > 0 is the test of a defined bin
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
  bool Has(const TH1 *h, int b) { return h && b >= 1 && b <= h->GetNbinsX() && h->GetBinError(b) > 0; }
  TString Val(const TH1 *h, int b, const char *fmt = "%.2f", const char *miss = "--") {
    return Has(h, b) ? TString(Form(fmt, h->GetBinContent(b))) : TString(miss); }
  TString Num(bool ok, double v, const char *fmt = "%.3f", const char *miss = "--") { return ok ? TString(Form(fmt, v)) : TString(miss); }
  bool HasP(const TProfile *p, int b, double minEff = 0) {
    return p && b >= 1 && b <= p->GetNbinsX() && p->GetBinEntries(b) > 0 && p->GetBinEffectiveEntries(b) >= minEff; }
  TH1D* ProjX(TH2 *h, const char *name) { TH1D *p = h->ProjectionX(name, 0, -1); p->SetDirectory(0); return p; }
  // unit-normalised copy for shape comparisons, normalised to every entry
  // (under- and overflow included); fold: the under- and overflow are drawn in
  // the edge bins (a share of exactly 1 is the overflow of a 0..1 axis)
  TH1D* Unit(const TH1 *h, const char *name, int rebin = 1, bool fold = false) {
    TH1D *u = Like(h, name);
    const int nb = h->GetNbinsX();
    for (int b = 1; b <= nb; ++b) { u->SetBinContent(b, h->GetBinContent(b)); u->SetBinError(b, h->GetBinError(b)); }
    if (fold) {
      u->SetBinContent(1, u->GetBinContent(1) + h->GetBinContent(0));
      u->SetBinContent(nb, u->GetBinContent(nb) + h->GetBinContent(nb+1));
    }
    if (rebin > 1 && nb % rebin == 0) u->Rebin(rebin);
    const double n = h->Integral(0, nb+1);
    if (n > 0) u->Scale(1./n);
    return u;
  }
  // num/den with Poisson errors, drawn where both are populated
  TH1D* Ratio(const TH1 *num, const TH1 *den, const char *name) {
    TH1D *r = Like(num, name);
    for (int b = 1; b <= r->GetNbinsX(); ++b) {
      const double n = num->GetBinContent(b), d = den->GetBinContent(b);
      if (n <= 0 || d < kMinDen) continue;
      r->SetBinContent(b, n/d); r->SetBinError(b, n/d*sqrt(1./n + 1./d));
    }
    return r;
  }
  // counts per unit x per normalisation unit (owned vertex, crossing)
  TH1D* PerUnit(const TH1 *h, const char *name, double norm) {
    TH1D *d = Like(h, name);
    for (int b = 1; b <= d->GetNbinsX(); ++b) {
      const double n = h->GetBinContent(b), w = h->GetBinWidth(b);
      if (n <= 0) continue;
      d->SetBinContent(b, n/w/norm); d->SetBinError(b, h->GetBinError(b)/w/norm);
    }
    return d;
  }
  // the integral of the bins inside [x1, x2] (edges a hair inside)
  double Sum(const TH1 *h, double x1, double x2) {
    return h ? h->Integral(h->GetXaxis()->FindFixBin(x1 + 1e-6), h->GetXaxis()->FindFixBin(x2 - 1e-6)) : 0; }
  // a TGraphErrors of the defined bins of h, x shifted by a fraction of the
  // bin (linear axis) or by a factor (log axis), so that seven methods do not
  // sit on top of each other
  TGraphErrors* Graph(const TH1 *h, double shift, bool logx, double xhi = 1e30) {
    TGraphErrors *g = new TGraphErrors();
    if (!h) return g;
    const TAxis *a = h->GetXaxis();
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
      if (!Has(h, b) || a->GetBinLowEdge(b) >= xhi) continue;
      const double x = logx ? GeoCentre(a, b)*(1 + shift) : a->GetBinCenter(b) + shift*a->GetBinWidth(b);
      const int n = g->GetN();
      g->SetPoint(n, x, h->GetBinContent(b)); g->SetPointError(n, 0, h->GetBinError(b));
    }
    return g;
  }
  bool HasPoints(const TH1 *h) {
    if (!h) return false;
    for (int b = 1; b <= h->GetNbinsX(); ++b) if (h->GetBinError(b) > 0) return true;
    return false;
  }
  // the extent of the defined points below xhi, for a frame
  // (points with an error above relMax x |value| do not open the frame: one
  // low-statistics point would otherwise flatten everything else)
  void Extent(const TH1 *h, double &dmin, double &dmax, double xhi = 1e30, double relMax = 1e30) {
    for (int b = 1; h && b <= h->GetNbinsX(); ++b) {
      if (!Has(h, b) || h->GetXaxis()->GetBinLowEdge(b) >= xhi) continue;
      if (h->GetBinError(b) > relMax*fabs(h->GetBinContent(b))) continue;
      dmin = std::min(dmin, h->GetBinContent(b) - h->GetBinError(b));
      dmax = std::max(dmax, h->GetBinContent(b) + h->GetBinError(b));
    }
  }
  // The frame of a reco/gen panel: linear from 0 while the points span less
  // than a factor 20, logarithmic otherwise.  A linear 0..8 put dup (reco/gen
  // 50-130 at 2.5 < |eta| < 3.0) off the page and jva and the oracles
  // (0.1-0.3) on the zero line, and that contrast is the result of the spike
  // figures.  Points with an error above half their value do not set it.
  bool RatioFrame(TH1D *const *rt, int n, double xhi, double &lo, double &hi) {
    double pmin = 1e30, pmax = -1e30;
    for (int m = 0; m < n; ++m) for (int b = 1; rt[m] && b <= rt[m]->GetNbinsX(); ++b) {
      if (!Has(rt[m], b) || rt[m]->GetXaxis()->GetBinLowEdge(b) >= xhi) continue;
      const double v = rt[m]->GetBinContent(b), e = rt[m]->GetBinError(b);
      if (v <= 0 || e > 0.5*v) continue;
      pmin = std::min(pmin, v); pmax = std::max(pmax, v + e);
    }
    if (pmin < 1e29 && pmax/pmin > 20) { lo = std::max(1e-3, 0.5*pmin); hi = 3*pmax; return true; }
    lo = 0; hi = std::max(2.4, std::min(8., 1.1*(pmax > 0 ? pmax : 1.)));
    return false;
  }
  // The frame: the default [lo, hi], opened up where the points need it, so
  // that the highest point stays below the fraction 'below' of the frame
  // height (the legend sits above it); never below 'floor'
  void OpenRange(double &lo, double &hi, double dmin, double dmax, double below, double floor = -1e30) {
    if (dmin > dmax) return;
    if (dmin < lo) lo = std::max(floor, dmin - 0.05*(hi - lo));
    if (dmax > lo + below*(hi - lo)) hi = lo + (dmax - lo)/below;
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
  TString gThr1, gThr2;
  void Thresholds(double x = 0.94, double y = 0.925, int align = 33, double size = 0.019, bool up = false) {
    TLatex t; t.SetNDC(); t.SetTextSize(size); t.SetTextAlign(align);
    t.DrawLatex(x, up ? y + 1.3*size : y, gThr1);
    t.DrawLatex(x, up ? y : y - 1.3*size, gThr2);
  }
  // the same in the upper pad of a tdrDiCanvas, whose frame ends at 0.92 of
  // the pad: drawGenSeed.C's 0.925 puts the first line across the frame edge
  void ThresholdsDi() { Thresholds(0.94, 0.905, 33, 0.022); }
  void Note(double x, double y, const char *s, double size = 0.026, int align = 11, int colour = kBlack) {
    TLatex t; t.SetNDC(); t.SetTextSize(size); t.SetTextAlign(align); t.SetTextColor(colour); t.DrawLatex(x, y, s); }
  void HLine(double x1, double x2, double y = 1., int colour = kGray+1, int style = kDashed) {
    TLine l; l.SetLineStyle(style); l.SetLineColor(colour); l.DrawLine(x1, y, x2, y); }
  TString LumiText(double pbinv) {
    if (pbinv >= 1)    return Form("%.3g pb^{-1}", pbinv);
    if (pbinv >= 1e-3) return Form("%.3g nb^{-1}", pbinv*1e3);
    return Form("%.3g #mub^{-1}", pbinv*1e6);
  }
  bool gJecOn = false;
  const char* PtCorr() { return gJecOn ? "p_{T}^{corr}" : "p_{T}^{raw}"; }

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
    if (gGs.IsNull()) printf("drawJVA: ghostscript (gs) not found, the PDFs are left uncropped\n");
    else printf("drawJVA: cropping with %s\n", gGs.Data());
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

  // -------------------------------------------------------------------
  // THE CORE-FIT PROTOCOL FOR A RESPONSE (drawGenSeed.C, from drawSeed4.C).
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
  // The fits are Poisson likelihoods, which need counts.  A histogram filled
  // with one common weight w (JVA.C fills both assignments of a barrel pair,
  // and both signs of a bisector projection, at 1/2) is counts x w, and is
  // fitted as counts, content/w, which moves no parameter; the amplitude is
  // scaled back.  A histogram of mixed weights gets a chi2 fit instead, said
  // once: ROOT 6.28's weighted likelihood ("WL") crashes when a parameter is
  // fixed, and the centred fits fix one.
  double CommonWeight(const TH1 *p) {      // 1: counts; w: one common weight; 0: mixed
    double w = -1;
    for (int i = 0; i <= p->GetNbinsX()+1; ++i) {
      const double c = p->GetBinContent(i), e = p->GetBinError(i);
      if (c <= 0) continue;
      const double wb = e*e/c;
      if (w < 0) w = wb; else if (fabs(wb - w) > 1e-6*w) return 0;
    }
    return w <= 0 ? 1. : w;
  }
  bool gSaidChi2 = false;
  // the histogram to fit (p itself, or a scaled copy owned by keep), its weight and the fit option
  TH1D* FitCopy(TH1D *p, std::unique_ptr<TH1D> &keep, double &w, const char *&opt) {
    w = CommonWeight(p);
    opt = w > 0 ? "QRN0L" : "QRN0";
    if (w <= 0) {
      if (!gSaidChi2) { printf("drawJVA: %s has mixed weights: chi2 core fits for such histograms\n", p->GetName()); gSaidChi2 = true; }
      w = 1; return p;
    }
    if (fabs(w - 1) < 1e-6) return p;
    keep.reset((TH1D*)p->Clone(Form("%s_cnt", p->GetName()))); keep->SetDirectory(0); keep->Scale(1./w);
    return keep.get();
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
    std::unique_ptr<TH1D> pcnt; double wcnt = 1; const char *opt = "QRN0L";
    TH1D *pf = FitCopy(p, pcnt, wcnt, opt);
    for (s.npass = 0; s.npass < 20; ++s.npass) {
      fg.SetRange(mu - 2*sig, mu + sig);
      fg.SetParameters(pf->GetMaximum(), mu, sig);
      if (pf->Fit(&fg, opt) != 0) break;
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
    s.fcore = wcnt*fg.GetParameter(0)*sig*sqrt(2*M_PI)/p->GetBinWidth(1)/s.n;
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
  // sigma_core/mu_core of a (pT, R) histogram against pT, on its own axis; and,
  // into rmsOut if given, rms/mean (the in-axis moments) of the same bins
  TH1D* SigCoreVsPt(TH2 *h, const char *name, TH1D **rmsOut = 0) {
    if (rmsOut) *rmsOut = 0;
    if (!h) return 0;
    TH1D *o = Like(h, name), *r = rmsOut ? Like(h, Form("%s_rms", name)) : 0;
    for (int b = 1; b <= o->GetNbinsX(); ++b) {
      const Shape s = DescribeBin(h, b);
      if (r && s.ok && s.mean > 0) {
        // the error of rms/mean: that of an rms, sqrt(1/2n) relative
        const double v = s.rms/s.mean;
        r->SetBinContent(b, v); r->SetBinError(b, std::max(v/sqrt(2*s.n), 1e-9));
      }
      if (!s.ok || !s.core || s.muC <= 0) continue;
      const double v = s.sigC/s.muC, e = sqrt(pow(s.sigCerr/s.muC,2) + pow(v*s.muCerr/s.muC,2));
      o->SetBinContent(b, v); o->SetBinError(b, std::max(e, 1e-9));
    }
    if (rmsOut) *rmsOut = r;
    return o;
  }

  // -------------------------------------------------------------------
  // THE PROTOCOL FOR A ZERO-CENTRED DISTRIBUTION (see the header): the same
  // iteration in the symmetric window mu +- 2 sigma, clipped to the axis.
  // nIndep > 0: the number of independent entries behind p when p holds each
  // of them more than once (a mirror sum, or a histogram filled at +x and -x);
  // the 200-entry threshold and the errors are then those of nIndep: the fit
  // errors, which count p's effective entries, are scaled by
  // sqrt(neff(p)/nIndep).
  struct SymShape {
    double n, neff, ovf;              // original entries (incl. under/overflow), effective entries, outside the axis
    double med, mean, rms, hw;        // median, in-axis mean and rms, (q84 - q16)/2
    double amp, muC, sigC, muCerr, sigCerr, tail;   // core; tail: beyond +-2.5 sigma_core
    int    npass; bool conv, ok, core, centred;
  };
  double NormCdf(double x) { return 0.5*(1 + TMath::Erf(x/M_SQRT2)); }
  double NormPdf(double x) { return exp(-0.5*x*x)/sqrt(2*M_PI); }
  // the rms of a Gaussian (mu, sigma) seen only inside [lo, hi]
  double TruncRms(double mu, double sig, double lo, double hi) {
    const double a = (lo - mu)/sig, b = (hi - mu)/sig, z = NormCdf(b) - NormCdf(a);
    if (z <= 1e-12) return sig;
    const double m1 = (NormPdf(a) - NormPdf(b))/z, v = 1 + (a*NormPdf(a) - b*NormPdf(b))/z - m1*m1;
    return sig*sqrt(std::max(v, 0.));
  }
  int gNDescS = 0, gNConvS = 0, gNWideS = 0, gNMuS = 0;
  SymShape DescribeSym(TH1D *p, bool centred = false, double nIndep = 0)
  {
    SymShape s; memset(&s, 0, sizeof(s));
    s.centred = centred;
    if (!p) return s;
    const int nb = p->GetNbinsX();
    const double all = p->Integral(0, nb+1), neffP = p->GetEffectiveEntries();
    s.n = nIndep > 0 ? nIndep : all;
    if (s.n < kMinN || all <= 0 || neffP <= 0) return s;
    s.neff = nIndep > 0 ? nIndep : neffP;
    const double eScale = nIndep > 0 ? sqrt(neffP/nIndep) : 1.;
    s.ovf = (p->GetBinContent(0) + p->GetBinContent(nb+1))/all;
    double pr[3] = {0.16, 0.5, 0.84}, q[3];
    p->GetQuantiles(3, q, pr);
    s.med = q[1]; s.hw = 0.5*(q[2] - q[0]);
    s.mean = p->GetMean(); s.rms = p->GetRMS();
    s.ok = true;
    const double xlo = p->GetXaxis()->GetXmin(), xhi = p->GetXaxis()->GetXmax();
    double mu = centred ? 0 : s.med, sig = s.hw > 0 ? s.hw : s.rms;
    if (sig <= 0) return s;
    TF1 fg(Form("fs_%s", p->GetName()), "gaus", xlo, xhi);
    std::unique_ptr<TH1D> pcnt; double wcnt = 1; const char *opt = "QRN0L";
    TH1D *pf = FitCopy(p, pcnt, wcnt, opt);
    for (s.npass = 0; s.npass < 20; ++s.npass) {
      fg.SetRange(std::max(xlo, mu - 2*sig), std::min(xhi, mu + 2*sig));
      fg.SetParameters(pf->GetMaximum(), mu, sig);
      if (centred) fg.FixParameter(1, 0.);
      if (pf->Fit(&fg, opt) != 0) break;
      const double mu1 = centred ? 0. : fg.GetParameter(1), sig1 = fabs(fg.GetParameter(2));
      if (sig1 <= 0) break;
      const bool c = fabs(sig1 - sig)/sig < 0.01;
      mu = c ? mu1 : 0.5*(mu + mu1); sig = c ? sig1 : 0.5*(sig + sig1);
      if (c) { s.conv = true; break; }
    }
    s.amp = wcnt*fg.GetParameter(0); s.muC = mu; s.sigC = sig;
    s.muCerr = centred ? 0. : fg.GetParError(1)*eScale; s.sigCerr = fg.GetParError(2)*eScale;
    // wider than the distribution: what a Gaussian of this width would show
    // inside the axis is wider than what is there
    const bool wide = TruncRms(mu, sig, xlo, xhi) > 1.1*s.rms, off = !centred && fabs(mu - s.med) > 0.5*sig;
    s.core = s.conv && !wide && !off;
    ++gNDescS; if (s.conv) ++gNConvS;
    if (s.conv && wide) ++gNWideS;
    if (s.conv && off) ++gNMuS;
    if (s.core) s.tail = 1. - (FracAbove(p, mu - 2.5*sig) - FracAbove(p, mu + 2.5*sig))/all;
    return s;
  }
  // the mirror sum h(x) + h(-x) (see the header).  An axis that starts at 0
  // with uniform bins is taken as already folded (|x| was filled) and is
  // mirrored onto -x..x, every entry twice as in the mirror sum; 0 for any
  // other axis that is not symmetric about zero
  bool gWarnedFold = false;
  TH1D* Symmetrise(const TH1D *p, const char *name, const char *what = "") {
    const TAxis *a = p->GetXaxis(); const int nb = a->GetNbins();
    const double tol = 1e-6*(a->GetXmax() - a->GetXmin());
    if (fabs(a->GetXmin()) < tol && a->GetXbins()->GetSize() == 0) {
      if (!gWarnedFold) { printf("drawJVA: WARNING %s has a y axis 0..%g: taken as folded (|x|) and mirrored\n", what, a->GetXmax()); gWarnedFold = true; }
      TH1D *s = new TH1D(name, "", 2*nb, -a->GetXmax(), a->GetXmax()); s->SetDirectory(0); s->Sumw2();
      for (int b = 1; b <= nb; ++b) {
        s->SetBinContent(nb + b, p->GetBinContent(b));     s->SetBinError(nb + b, p->GetBinError(b));
        s->SetBinContent(nb + 1 - b, p->GetBinContent(b)); s->SetBinError(nb + 1 - b, p->GetBinError(b));
      }
      s->SetBinContent(0, p->GetBinContent(nb+1)); s->SetBinContent(2*nb+1, p->GetBinContent(nb+1));
      s->SetBinError(0, p->GetBinError(nb+1));     s->SetBinError(2*nb+1, p->GetBinError(nb+1));
      s->ResetStats();
      return s;
    }
    for (int b = 1; b <= nb; ++b) if (fabs(a->GetBinLowEdge(b) + a->GetBinUpEdge(nb+1-b)) > tol) return 0;
    TH1D *s = Like(p, name);
    for (int b = 0; b <= nb+1; ++b) {
      const int m = nb + 1 - b;
      s->SetBinContent(b, p->GetBinContent(b) + p->GetBinContent(m));
      s->SetBinError(b, sqrt(pow(p->GetBinError(b), 2) + pow(p->GetBinError(m), 2)));
    }
    s->ResetStats();
    return s;
  }
  // Already symmetric: every entry was filled at +x and -x (JVA.C fills the
  // bisector projections both ways at weight 1/2).  Such a histogram IS the
  // mirror sum, but its effective entries count every dijet twice.  The
  // tolerance allows for the odd fill that rounding put one bin off.
  bool IsSymmetric(const TH1D *p) {
    const TAxis *a = p->GetXaxis(); const int nb = a->GetNbins();
    const double tol = 1e-6*(a->GetXmax() - a->GetXmin());
    for (int b = 1; b <= nb; ++b) if (fabs(a->GetBinLowEdge(b) + a->GetBinUpEdge(nb+1-b)) > tol) return false;
    double sum = 0, dif = 0;
    for (int b = 0; b <= nb+1; ++b) { sum += fabs(p->GetBinContent(b)); dif += fabs(p->GetBinContent(b) - p->GetBinContent(nb+1-b)); }
    return sum > 0 && dif < 2e-3*sum;
  }
  // the symmetrised y distribution of one x bin, described: as it is when it
  // was filled symmetrically (nIndep = half its effective entries), else its
  // mirror sum (nIndep = the projection's effective entries); the projection
  // itself, uncentred, when the axis can be neither (said once)
  bool gWarnedMirror = false, gSaidSym = false;
  SymShape SymBin(TH2 *h, int b, TH1D **keep = 0) {
    SymShape s; memset(&s, 0, sizeof(s));
    if (!h || b < 1 || b > h->GetNbinsX()) return s;
    std::unique_ptr<TH1D> p(h->ProjectionY(Form("_sy_%s_%d", h->GetName(), b), b, b)); p->SetDirectory(0);
    if (p->Integral(0, p->GetNbinsX()+1) <= 0) return s;
    const double neffP = p->GetEffectiveEntries();
    TH1D *m = 0;
    if (IsSymmetric(p.get())) {
      if (!gSaidSym) { printf("drawJVA: %s is filled symmetrically (+x and -x): used as it is, with half its effective entries\n", h->GetName()); gSaidSym = true; }
      m = (TH1D*)p->Clone(Form("_sm_%s_%d", h->GetName(), b)); m->SetDirectory(0);
      s = DescribeSym(m, true, 0.5*neffP);
    } else if ((m = Symmetrise(p.get(), Form("_sm_%s_%d", h->GetName(), b), h->GetName()))) {
      s = DescribeSym(m, true, neffP);
    } else {
      if (!gWarnedMirror) { printf("drawJVA: WARNING %s has a y axis that cannot be mirrored; its rms is about the mean\n", h->GetName()); gWarnedMirror = true; }
      s = DescribeSym(p.get(), false);
    }
    if (keep) *keep = m ? m : (TH1D*)p->Clone(Form("_sk_%s_%d", h->GetName(), b));
    else delete m;
    if (keep && *keep) (*keep)->SetDirectory(0);
    return s;
  }
  // the width of a zero-centred distribution vs the x axis of a TH2: rms and
  // core sigma (not symmetrised: a MET component is its own distribution)
  void WidthVsX(TH2 *h, const char *name, TH1D *&hrms, TH1D *&hcore) {
    hrms = hcore = 0;
    if (!h) return;
    hrms = Like(h, Form("%s_rms", name)); hcore = Like(h, Form("%s_core", name));
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
      std::unique_ptr<TH1D> p(h->ProjectionY(Form("_wx_%s_%d", name, b), b, b)); p->SetDirectory(0);
      const SymShape s = DescribeSym(p.get());
      if (!s.ok) continue;
      hrms->SetBinContent(b, s.rms); hrms->SetBinError(b, std::max(s.rms/sqrt(2*std::max(s.neff, 1.)), 1e-9));
      if (s.core) { hcore->SetBinContent(b, s.sigC); hcore->SetBinError(b, std::max(s.sigCerr, 1e-9)); }
    }
  }

  // -------------------------------------------------------------------
  // dijet numbers
  struct Pt { bool ok; double v, e, n; };
  // <A> of one pT_avg bin and R = (1 + <A>)/(1 - <A>)
  Pt RFromA(TH2 *h, int b) {
    Pt r = {false, 0, 0, 0};
    if (!h || b < 1 || b > h->GetNbinsX()) return r;
    std::unique_ptr<TH1D> p(h->ProjectionY(Form("_a_%s_%d", h->GetName(), b), b, b)); p->SetDirectory(0);
    r.n = p->GetEffectiveEntries();
    if (r.n < kMinMean) return r;
    const double a = p->GetMean(), ea = p->GetRMS()/sqrt(r.n);
    if (a <= -1 || a >= 1) return r;
    r.v = (1 + a)/(1 - a); r.e = std::max(2*ea/((1 - a)*(1 - a)), 1e-9); r.ok = true;
    return r;
  }
  // the median of one pT_avg bin, with the Gaussian-equivalent error
  Pt MedianBin(TH2 *h, int b) {
    Pt r = {false, 0, 0, 0};
    if (!h || b < 1 || b > h->GetNbinsX()) return r;
    std::unique_ptr<TH1D> p(h->ProjectionY(Form("_m_%s_%d", h->GetName(), b), b, b)); p->SetDirectory(0);
    r.n = p->GetEffectiveEntries();
    if (r.n < kMinMean) return r;
    double pr[3] = {0.16, 0.5, 0.84}, q[3];
    p->GetQuantiles(3, q, pr);
    r.v = q[1]; r.e = std::max(1.2533*0.5*(q[2] - q[0])/sqrt(r.n), 1e-9); r.ok = r.v > 0;
    return r;
  }
  // sigma = sqrt(P^2 - Q^2)/sqrt(2) into bin b of o.  Returns false when both
  // widths exist but Q >= P: the imbalance across the axis is as wide as the
  // one along it, the bin HAS no JER, and it is flagged in fail (content 1) -
  // a measured failure, which the reach must see, not an unmeasured bin that it
  // may skip.
  bool FillJer(TH1D *o, TH1D *fail, int b, double P, double dP, double Q, double dQ) {
    if (!(P > Q) || Q < 0) { if (fail) fail->SetBinContent(b, 1); return false; }
    const double v = sqrt((P*P - Q*Q)/2), e = sqrt(P*P*dP*dP + Q*Q*dQ*dQ)/(2*v);
    o->SetBinContent(b, v); o->SetBinError(b, std::max(e, 1e-9));
    return true;
  }
  // The reach (see the header): the lower edge of the lowest bin with
  // |ratio - 1| < kJerTol whose next two measured bins (fewer where the
  // measurement ends) are not significantly outside the band,
  // |ratio - 1| < kJerTol + 2 sigma.  A bin with Q >= P where the truth exists
  // (fail) is a measured bin that fails.  kReachNone if bins were measured but
  // none qualifies, kReachNA if nothing was measured.  The first version asked
  // the next two bins for |ratio - 1| < kJerTol strictly, and skipped the
  // Q >= P bins: on a toy with 1000-5000 dijets per bin and a 5% offset that
  // scattered the reach of one and the same truth between 5 and 28 GeV, and a
  // Q >= P bin in the middle of the range went unnoticed.
  const double kReachNone = -1, kReachNA = -2;
  double Reach(const TH1D *ratio, const TH1D *fail = 0, const TH1D *truth = 0, double xhi = kPtDrawHi) {
    if (!ratio) return kReachNA;
    std::vector<int> bins; std::vector<bool> failed;
    for (int b = 1; b <= ratio->GetNbinsX(); ++b) {
      if (ratio->GetXaxis()->GetBinLowEdge(b) >= xhi) continue;
      if (Has(ratio, b)) { bins.push_back(b); failed.push_back(false); }
      else if (fail && fail->GetBinContent(b) > 0 && Has(truth, b)) { bins.push_back(b); failed.push_back(true); }
    }
    if (bins.empty()) return kReachNA;
    auto dev = [&](size_t k) { return fabs(ratio->GetBinContent(bins[k]) - 1); };
    for (size_t k = 0; k < bins.size(); ++k) {
      if (failed[k] || !(dev(k) < kJerTol)) continue;
      bool ok = true;
      for (size_t j = k + 1; j < std::min(bins.size(), k + 3); ++j)
        ok = ok && !failed[j] && dev(j) < kJerTol + 2*ratio->GetBinError(bins[j]);
      if (ok) return ratio->GetXaxis()->GetBinLowEdge(bins[k]);
    }
    return kReachNone;
  }
  // "none": measured, never within the tolerance; "--": nothing to measure
  TString ReachStr(double r) { return r > 0 ? TString(Form("%g", r)) : r == kReachNone ? TString("none") : TString("--"); }
  // The pT_avg axis of the JER, merged: from the first edge at or above
  // 'from', every 'merge' fine bins make one (the bins below 'from' are kept).
  // The JER needs 200 same-bin dijets per bin, and the kPtBinsExt steps of
  // ~20% leave nearly every bin below that on the samples of this study; two
  // fine bins per point (3-5, 5-8, 8-12, 12-18, ...) are still finer than the
  // resolution changes.  merge <= 1 returns h itself.
  TH2* MergeX(TH2 *h, int merge, const char *name, double from = 3.) {
    if (!h || merge <= 1) return h;
    const TAxis *a = h->GetXaxis();
    std::vector<double> e;
    int b = 1;
    for (; b <= a->GetNbins() && a->GetBinLowEdge(b) < from - 1e-6; ++b) e.push_back(a->GetBinLowEdge(b));
    for (; b <= a->GetNbins(); b += merge) e.push_back(a->GetBinLowEdge(b));
    e.push_back(a->GetXmax());
    const TAxis *y = h->GetYaxis();
    TH2D *m = y->GetXbins()->GetSize() > 0
      ? new TH2D(name, h->GetTitle(), (int)e.size() - 1, e.data(), y->GetNbins(), y->GetXbins()->GetArray())
      : new TH2D(name, h->GetTitle(), (int)e.size() - 1, e.data(), y->GetNbins(), y->GetXmin(), y->GetXmax());
    m->SetDirectory(0); m->Sumw2();
    for (int i = 0; i <= a->GetNbins() + 1; ++i) {
      const int im = i == 0 ? 0 : i > a->GetNbins() ? m->GetNbinsX() + 1 : m->GetXaxis()->FindFixBin(a->GetBinCenter(i));
      for (int j = 0; j <= y->GetNbins() + 1; ++j) {
        m->SetBinContent(im, j, m->GetBinContent(im, j) + h->GetBinContent(i, j));
        m->SetBinError(im, j, sqrt(pow(m->GetBinError(im, j), 2) + pow(h->GetBinError(i, j), 2)));
      }
    }
    m->SetEntries(h->GetEntries());
    return m;
  }

  // -------------------------------------------------------------------
  // the LaTeX writers.  A figure lists only PDFs that were actually saved.
  std::ofstream gTexP, gTexT;
  TString TexEscape(const char *s) { TString t(s); t.ReplaceAll("_", "\\_"); return t; }
  // A section is written when it ends (at the next TexSection or at the end),
  // text first and its figures after it, and only if it produced a figure (or
  // was forced: the summary, whose content is in the tables).  References are
  // written as tokens, @fig:<label>@ and @tab:<label>@, and resolved then: to
  // "figure~\ref{fig:jva_<label>}" if that figure was written, and to words
  // saying it was not otherwise, so that a missing histogram leaves no "??".
  struct Sec { TString title, text, figs; bool open, force; };
  Sec gSec = {"", "", "", false, false};
  std::set<std::string> gFigDone, gTabDone;
  TString Resolve(TString t) {
    for (const char *kind : {"fig", "tab"}) {
      const TString key = Form("@%s:", kind);
      Ssiz_t i;
      while ((i = t.Index(key)) != kNPOS) {
        const Ssiz_t j = t.Index("@", i + key.Length());
        if (j == kNPOS) break;
        const TString lab = t(i + key.Length(), j - i - key.Length());
        const bool done = TString(kind) == "fig" ? gFigDone.count(lab.Data()) > 0 : gTabDone.count(lab.Data()) > 0;
        t.Replace(i, j - i + 1, done ? Form("%s~\\ref{%s:jva_%s}", TString(kind) == "fig" ? "figure" : "table", kind, lab.Data())
                                     : (TString(kind) == "fig" ? "a figure that was not produced" : "a table that was not produced"));
      }
    }
    return t;
  }
  void FlushSection() {
    if (!gSec.open) return;
    if (gSec.figs.IsNull() && !gSec.force) printf("drawJVA: section '%s' produced no figure; its text is not written\n", gSec.title.Data());
    else gTexP << "\n\\subsection{" << gSec.title << "}\n" << Resolve(gSec.text) << "\n" << Resolve(gSec.figs);
    gSec = {"", "", "", false, false};
  }
  void TexSection(const char *title, const char *text, bool force = false) {
    FlushSection();
    gSec = {title, text, "", true, force};
  }
  // rows: optional sizes of the row groups over the INPUT list, so that a
  // caption can say "first row: ..." (drawGenSeed.C's TexFigure)
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
    TString f = "\\begin{figure}[htbp]\\centering\n";
    size_t inRow = 0;
    for (size_t k = 0; k < n; ++k) {
      ++inRow;
      const bool brk = k + 1 < n && (grp[k+1] != grp[k] || inRow == perRow);
      if (brk) inRow = 0;
      f += Form("\\includegraphics[width=%.2f\\textwidth]{\\plotdir %s.pdf}%s\n", width, ok[k].Data(),
                k + 1 == n ? "" : brk ? "\\\\" : "\\hfil");
    }
    f += "\\caption{" + caption + "}\n\\label{fig:jva_" + label + "}\n\\end{figure}\n\n";
    gFigDone.insert(label);
    if (gSec.open) gSec.figs += f; else gTexP << Resolve(f);
  }
  // A table: LaTeX to doc/jva_tables.tex, the same cells to stdout under a
  // plain-text header of the same column width.  Buffered, and written only
  // if it got a row.
  int gTabW = 9, gTabRows = 0;
  TString gTab, gTabLabel;
  void TexTableBegin(const char *caption, const char *label, const char *cols, const char *header,
                     const std::vector<TString> &txtCols, int w = 9, const char *size = "\\small") {
    gTabW = w; gTabRows = 0; gTabLabel = label;
    gTab = Form("\\begin{table}[htbp]\\centering%s\n\\caption{", size);
    gTab += TString(caption) + "}\n\\label{tab:jva_" + label + "}\n\\begin{tabular}{" + cols + "}\\hline\n" + header + " \\\\ \\hline\n";
    printf("--- table tab:jva_%s\n  %-24s", label, "");
    for (const TString &c : txtCols) printf(" %*s", w, c.Data());
    printf("\n");
  }
  void TexRow(const TString &texLabel, const TString &txtLabel, const std::vector<TString> &cells, int w = -1) {
    if (w < 0) w = gTabW;
    gTab += texLabel;
    for (const TString &c : cells) gTab += " & " + c;
    gTab += " \\\\\n";
    ++gTabRows;
    printf("  %-24s", txtLabel.Data());
    for (const TString &c : cells) printf(" %*s", w, c.Data());
    printf("\n");
  }
  void TexHline() { gTab += "\\hline\n"; }
  void TexTableEnd() {
    if (!gTabRows) { printf("  (no rows: table tab:jva_%s not written)\n", gTabLabel.Data()); return; }
    // a trailing \hline from the last group is not doubled
    if (gTab.EndsWith("\\hline\n")) gTab.Remove(gTab.Length() - 7);
    gTexT << gTab << "\\hline\n\\end{tabular}\n\\end{table}\n\n";
    gTabDone.insert(gTabLabel.Data());
  }
}

void drawJVA(const char *fname = "rootfiles/JVA_v1.root", const char *tag = "_v1", int jerMerge = 2)
{
  gF = TFile::Open(fname);
  if (!gF || gF->IsZombie()) { printf("drawJVA: cannot open %s\n", fname); return; }
  setTDRStyle();
  gErrorIgnoreLevel = kWarning;   // no "pdf file created" between the tables
  gSystem->mkdir("plots", kTRUE); gSystem->mkdir("doc", kTRUE);
  FindGs();

  // ---- counters, header, knobs ------------------------------------------
  double events = 0, npu = 0, nvtx = 0, nown = 0;
  if (TH1 *hc = GetH("hcount")) {
    auto at = [hc](const char *l) { const int b = hc->GetXaxis()->FindFixBin(l); return (b > 0 && b <= hc->GetNbinsX()) ? hc->GetBinContent(b) : 0.; };
    if (at("jobs") > 0) gJobs = at("jobs");
    events = at("events_a") + at("events_b"); npu = at("npu_a") + at("npu_b");
    nvtx = at("vertices"); nown = at("owned");
    printf("=== hcount: %.0f + %.0f events, %.0f entries, %g jobs, sum N_PU %.0f, <PU> = %.1f, L = %.3g ub-1; "
           "%.0f usable vertices (%.2f per crossing), %.0f owned (%.3f)\n",
           at("events_a"), at("events_b"), at("entries"), gJobs, npu, events > 0 ? npu/events : 0,
           npu/(kSigmaInel_mb*1e9)*1e6, nvtx, events > 0 ? nvtx/events : 0, nown, nvtx > 0 ? nown/nvtx : 0);
  }
  extraText = "Simulation"; extraText2 = "";
  lumi_136TeV = Form("%s, #LTPU#GT = %.0f", LumiText(npu/(kSigmaInel_mb*1e9)).Data(), events > 0 ? npu/events : 0.);
  gOpts = GetH("hopts");
  if (!gOpts) printf("drawJVA: WARNING no hopts in this file: every knob in the labels is its default\n");
  if (gOpts) {
    printf("=== knobs (hopts / jobs):");
    for (int b = 1; b <= gOpts->GetNbinsX(); ++b) printf(" %s=%g", gOpts->GetXaxis()->GetBinLabel(b), gOpts->GetBinContent(b)/gJobs);
    printf("\n");
  }
  const double R = Knob("R", 0.4), ptMinCluster = Knob("ptMinCluster", 1.), etaFwd = Knob("etaFwd", 2.5), trkMinPt = Knob("trkMinPt", 0.5);
  const double sigma0 = Knob("sigma0", 1.), sigmaK = Knob("sigmaK", 1.), tauTrk = Knob("tauTrk", 4.), tauAll = Knob("tauAll", 12.);
  const double ptMinJet = Knob("ptMinJet", 3.), t1Min = Knob("t1Min", 3.), dphiMin = Knob("dphiMin", 2.7);
  const double dzOracle = Knob("dzOracle", 0.2), dzVertexGen = Knob("dzVertexGen", 0.05);
  const TString puppiVersion = GetText("puppiVersion"), jecFile = GetText("jecFile"), backend = GetText("backend");
  gJecOn = !jecFile.IsNull() && jecFile != "none";
  printf("=== jets in %s%s%s; PUPPI %s; backend %s\n", gJecOn ? "p_T^corr, JEC " : "p_T^raw (no JEC loaded: every 'corr' is raw)",
         gJecOn ? jecFile.Data() : "", "", puppiVersion.IsNull() ? "?" : puppiVersion.Data(), backend.IsNull() ? "?" : backend.Data());
  if (!gJecOn) printf("drawJVA: WARNING no JEC in this file: the dijet part (sections 4 and 5) is uncorrected\n");
  gThr1 = Form("anti-k_{T} R = %g from %g GeV, |#eta| > %g vertex-blind, tracks > %g GeV", R, ptMinCluster, etaFwd, trkMinPt);
  gThr2 = Form("#sigma_{v}^{2} = %g^{2} + %g^{2} S_{v}, #tau = %g / %g; dijets %s > %g GeV, #Delta#phi > %g",
               sigma0, sigmaK, tauTrk, tauAll, PtCorr(), ptMinJet, dphiMin);
  const double perX = events > 0 ? 1./events : 0;

  gTexP.open("doc/jva_plots.tex"); gTexT.open("doc/jva_tables.tex");
  gTexP << "% doc/jva_plots.tex - generated by drawJVA.C from " << fname << " (tag " << tag << "); do not edit\n"
        << "\\providecommand{\\plotdir}{../plots/}\n\\providecommand{\\pt}{p_{\\mathrm{T}}}\n";
  gTexT << "% doc/jva_tables.tex - generated by drawJVA.C from " << fname << " (tag " << tag << "); do not edit\n"
        << "\\providecommand{\\pt}{p_{\\mathrm{T}}}\n";
  TString thrTex = Form("The sample has %.0f crossings at $\\langle N_{\\mathrm{PU}}\\rangle=%.0f$ with %.1f usable vertices per "
                        "crossing, %.0f\\%% of them owned by an interaction (within %g\\,cm). Anti-$k_{\\mathrm{T}}$ $R=%g$ from "
                        "%g\\,GeV; candidates beyond $|\\eta|=%g$ are vertex-blind; track candidates need %g\\,GeV in the cone; "
                        "the optimiser uses $\\sigma_v^2=%g^2+%g^2\\,S_v$\\,GeV$^2$ and the prices $\\tau_{\\mathrm{trk}}=%g$, "
                        "$\\tau_{\\mathrm{all}}=%g$; the oracles accept a vertex within %g\\,cm.",
                        events, events > 0 ? npu/events : 0., events > 0 ? nvtx/events : 0., nvtx > 0 ? 100*nown/nvtx : 0.,
                        dzVertexGen, R, ptMinCluster, etaFwd, trkMinPt, sigma0, sigmaK, tauTrk, tauAll, dzOracle);
  thrTex += gJecOn ? Form(" Jets are corrected with \\texttt{%s}.", TexEscape(jecFile).Data())
                   : TString(" No jet energy correction was loaded, so every $\\pt^{\\mathrm{corr}}$ below is the uncorrected $\\pt^{\\mathrm{raw}}$.");
  if (!puppiVersion.IsNull()) thrTex += Form(" PUPPI weights from \\texttt{puppi.h} version \\texttt{%s}.", TexEscape(puppiVersion).Data());

  // the numbers the summary collects on the way (NaN-free: a flag says whether they exist)
  bool   haveMet[kNM] = {false}, haveMetRes[kNM] = {false}, haveMetCore[kNM] = {false}, haveResCore[kNM] = {false};
  double metRms[kNM] = {0}, metCore[kNM] = {0}, resCore[kNM] = {0}, resRms[kNM] = {0};
  bool   haveSpike[kNM] = {false}, haveOwn[kNM] = {false};
  double spike10[kNM] = {0}, own1[kNM] = {0};
  bool   haveCorrect[kNM] = {false};
  double correct[kNM] = {0}, dropAllRef = -1;
  bool   haveDij[kNM][kNDC] = {{false}};
  double dijPerX[kNM][kNDC] = {{0}};
  bool   haveL2[kNM][kNDC] = {{false}}, haveL2a1[kNM][kNDC] = {{false}};
  double l2mpf[kNM][kNDC] = {{0}}, l2mpfa1[kNM][kNDC] = {{0}};
  int    l2n[kNM][kNDC] = {{0}}, l2na1[kNM][kNDC] = {{0}};
  double reachMpfx[kNM][kNDC][kNJB], reachDb[kNM][kNDC][kNJB];
  for (int m = 0; m < kNM; ++m) for (int c = 0; c < kNDC; ++c) for (int k = 0; k < kNJB; ++k) reachMpfx[m][c][k] = reachDb[m][c][k] = kReachNA;

  // ======================================================================
  // 1. MET per method and vertex class
  // ======================================================================
  printf("\n=== 1. MET per method and vertex class ===\n");
  TexSection("Missing transverse momentum", Form("The missing transverse momentum of every usable vertex, "
             "$\\vec{p}_{\\mathrm{T}}^{\\,\\mathrm{miss}}=-\\sum w_i\\,\\vec{p}_{\\mathrm{T},i}$ over the candidates of its per-vertex "
             "event under each method, for all vertices, the primary vertex and the others (@fig:met_dist@); "
             "its width against $N_{\\mathrm{PU}}$ and against $S_v$, the scalar $\\pt$ sum of the vertex-resolved part "
             "(@fig:met_width@); and the residual against $-\\sum\\vec{p}_{\\mathrm{T}}$ of the generated "
             "particles of the owner interaction with $|\\eta|<5$, the imbalance a perfect detector would see. The rms is the "
             "in-axis rms; the core is a Gaussian iterated to convergence in $\\mu\\pm2\\sigma$, failed below %.0f entries, "
             "without convergence, when it is wider than the distribution or off its median by more than $\\sigma/2$. At the primary vertex \\texttt{dup} and \\texttt{lv} are the same "
             "event, the standard reconstruction. %s", kMinN, thrTex.Data()));
  {
    std::vector<TString> figs; std::vector<int> rows;
    SymShape sxy[kNM][kNC], sres[kNM][kNC];
    double mabs[kNM][kNC]; bool habs[kNM][kNC];
    memset(sxy, 0, sizeof(sxy)); memset(sres, 0, sizeof(sres)); memset(mabs, 0, sizeof(mabs)); memset(habs, 0, sizeof(habs));
    const char *kind[3]  = {"hmetxy", "hmet", "hmetres"}, *kpdf[3] = {"metxy", "metabs", "metres"};
    const char *kxlab[3] = {"p_{T,x}^{miss}, p_{T,y}^{miss} [GeV]", "|p_{T}^{miss}| [GeV]", "p_{T}^{miss} - p_{T}^{miss,true} (x, y) [GeV]"};
    const char *kylab[3] = {"Fraction of components", "Fraction of vertices", "Fraction of components"};
    for (int k = 0; k < 3; ++k) {
      for (int c = 0; c < kNC; ++c) {
        TH1D *h[kNM] = {0}, *u[kNM] = {0}; double ymax = 0; int nu = 0;
        for (int m = 0; m < kNM; ++m) {
          h[m] = AsD(GetT<TH1>(Form("met/%s/%s_%s", MTag(m), kind[k], CTag(c))), Form("met%d%d%d", k, c, m));
          if (!h[m] || h[m]->Integral(0, h[m]->GetNbinsX()+1) <= 0) { h[m] = 0; continue; }
          if (k == 0) sxy[m][c] = DescribeSym(h[m]);
          if (k == 1) { habs[m][c] = true; mabs[m][c] = h[m]->GetMean(); }
          if (k == 2) sres[m][c] = DescribeSym(h[m]);
          const int nb = h[m]->GetNbinsX();
          u[m] = Unit(h[m], Form("umet%d%d%d", k, c, m), nb >= 400 ? 4 : nb >= 200 ? 2 : 1);
          ymax = std::max(ymax, u[m]->GetMaximum()); ++nu;
        }
        if (!nu) continue;
        TH1D *any = 0; for (int m = 0; m < kNM; ++m) if (u[m]) { any = u[m]; break; }
        const double xlo = k == 1 ? 0. : std::max(-100., AxLo(any)), xhi = std::min(k == 1 ? 150. : 100., any->GetXaxis()->GetXmax());
        // the legend (a header and seven entries) needs the top ~40% of the pad
        TH1D *frame = tdrHist(Form("fmet%d%d", k, c), kylab[k], ymax*2e-5, ymax*5000, kxlab[k], xlo, xhi);
        frame->GetXaxis()->SetNdivisions(505);
        std::unique_ptr<TCanvas> cv(tdrCanvas(Form("c_met%d%d", k, c), frame, 8, 11, kSquare));
        cv->SetLogy();
        TLegend *leg = tdrLeg(0.40, 0.87 - 0.036*(nu + 1), 0.94, 0.87);
        leg->SetTextSize(0.024); leg->SetHeader(CLabel(c));
        for (int m = 0; m < kNM; ++m) if (u[m]) {
          tdrDraw(u[m], "HIST", kNone, MColour[m], MLine[m], -1, kNone, 0, 0, m == kJva ? 3 : 2);
          TString lab;
          if (k == 1) lab = Form("%s: #LT|p_{T}^{miss}|#GT %.1f GeV", MShort(m), mabs[m][c]);
          else {
            const SymShape &s = k == 0 ? sxy[m][c] : sres[m][c];
            lab = s.ok ? Form("%s: rms %.1f, #sigma_{c} %s GeV", MShort(m), s.rms, s.core ? Form("%.1f", s.sigC) : "fail")
                       : Form("%s: rms %.1f GeV (%.0f entries)", MShort(m), h[m]->GetRMS(), h[m]->Integral(0, h[m]->GetNbinsX()+1));
          }
          leg->AddEntry(u[m], lab, "L");
        }
        Thresholds();
        fixOverlay();
        figs.push_back(SavePdf(cv.get(), Form("jva_%s_%s", kpdf[k], CTag(c))));
      }
      rows.push_back(kNC);   // one row per kind: the group is the three classes of the input list below
    }
    // input list order is k-major, three classes each: rows {3,3,3} of the input
    std::vector<TString> order;
    for (int k = 0; k < 3; ++k) for (int c = 0; c < kNC; ++c) order.push_back(Form("jva_%s_%s", kpdf[k], CTag(c)));
    // the table
    TexTableBegin("Missing transverse momentum per method and vertex class: vertices (components / 2), in-axis rms and core "
                  "$\\sigma$ of the components, mean $|\\vec{p}_{\\mathrm{T}}^{\\,\\mathrm{miss}}|$, and rms and core $\\sigma$ of the "
                  "residual against the owner interaction's true imbalance (owned vertices); in GeV. \\emph{fail}: the core "
                  "fails the protocol.", "met", "llrrrrrr",
                  "method & vertices & $N$ & rms & $\\sigma_{\\mathrm{core}}$ & $\\langle|p_{\\mathrm{T}}^{\\mathrm{miss}}|\\rangle$ & "
                  "rms$_{\\mathrm{res}}$ & $\\sigma_{\\mathrm{core,res}}$",
                  {"N", "rms", "core", "<|MET|>", "rms_res", "core_res"});
    for (int m = 0; m < kNM; ++m) {
      bool any = false;
      for (int c = 0; c < kNC; ++c) any = any || sxy[m][c].ok || sres[m][c].ok || habs[m][c];
      if (!any) continue;
      for (int c = 0; c < kNC; ++c) {
        const SymShape &a = sxy[m][c], &r = sres[m][c];
        TexRow(Form("%s & %s", c == 0 ? MTex(m) : "", CTex(c)), Form("%s/%s", MTag(m), CTag(c)),
               {Num(a.ok, a.n/2, "%.0f"), Num(a.ok, a.rms, "%.2f"), a.ok ? Num(a.core, a.sigC, "%.2f", "fail") : TString("--"),
                Num(habs[m][c], mabs[m][c], "%.2f"), Num(r.ok, r.rms, "%.2f"), r.ok ? Num(r.core, r.sigC, "%.2f", "fail") : TString("--")});
      }
      if (m + 1 < kNM) TexHline();
      haveMet[m] = sxy[m][0].ok; metRms[m] = sxy[m][0].rms; haveMetCore[m] = sxy[m][0].core; metCore[m] = sxy[m][0].sigC;
      haveMetRes[m] = sres[m][0].ok; resRms[m] = sres[m][0].rms; haveResCore[m] = sres[m][0].core; resCore[m] = sres[m][0].sigC;
    }
    TexTableEnd();
    TString cap;
    for (int m : {kDup, kLv, kJva, kOPart, kNobody}) if (haveMet[m])
      cap += Form("%s%s %.1f\\,GeV", cap.IsNull() ? "Over all vertices the rms of a component is " : ", ", MTex(m), metRms[m]);
    if (!cap.IsNull()) cap += ".";
    TString capr;
    for (int m : {kDup, kLv, kJva, kOPart, kNobody}) if (haveMetRes[m])
      capr += Form("%s%s %s", capr.IsNull() ? " The core width of the residual against the truth is " : ", ", MTex(m),
                   haveResCore[m] ? Form("%.1f\\,GeV", resCore[m]) : "not defined");
    if (!capr.IsNull()) capr += ".";
    bool have[kNM]; for (int m = 0; m < kNM; ++m) have[m] = haveMet[m] || haveMetRes[m];
    TexFigure(order, Form("Missing transverse momentum per vertex, unit-normalised: first row the $x$ and $y$ components, second "
              "row the magnitude, third row the residual of the components against the owner interaction's true imbalance "
              "(owned vertices), for all vertices, the primary vertex and the other vertices (left to right), for %s. The legend "
              "quotes the rms and the core $\\sigma$ (or the mean magnitude). %s%s", MWho(have).Data(), cap.Data(), capr.Data()),
              "met_dist", 0.32, rows);

    // widths against N_PU and S_v
    std::vector<TString> wfigs; TString capw;
    struct WPlot { const char *hist, *pdf, *xlab, *what; bool logx; int rebin; bool model; };
    const WPlot wp[3] = {{"hmetxy_npu", "jva_metw_npu", "N_{PU}", "p_{T}^{miss} component", false, 5, false},
                         {"hmetxy_s", "jva_metw_s", "S_{v} [GeV]", "p_{T}^{miss} component", true, 1, false},
                         {"hmetres_s", "jva_metresw_s", "S_{v} [GeV]", "residual against truth", true, 1, true}};
    for (int k = 0; k < 3; ++k) {
      TH1D *wr[kNM] = {0}, *wc[kNM] = {0}; int nw = 0;
      for (int m = 0; m < kNM; ++m) {
        TH2 *h = Get2(Form("met/%s/%s", MTag(m), wp[k].hist));
        if (!h) continue;
        TH2 *hh = h;
        if (wp[k].rebin > 1 && h->GetNbinsX() % wp[k].rebin == 0 && h->GetXaxis()->GetXbins()->GetSize() == 0)
          hh = (TH2*)h->RebinX(wp[k].rebin, Form("%s_rb%d", h->GetName(), m));
        WidthVsX(hh, Form("w%d%d", k, m), wr[m], wc[m]);
        if (HasPoints(wr[m])) ++nw; else wr[m] = wc[m] = 0;
      }
      if (!nw) continue;
      // x range: where any method has a width
      double x1 = 1e30, x2 = -1e30, r1 = 1e30, r2 = -1e30, c1 = 1e30, c2 = -1e30;
      for (int m = 0; m < kNM; ++m) if (wr[m]) {
        for (int b = 1; b <= wr[m]->GetNbinsX(); ++b) if (Has(wr[m], b)) {
          x1 = std::min(x1, wr[m]->GetXaxis()->GetBinLowEdge(b)); x2 = std::max(x2, wr[m]->GetXaxis()->GetBinUpEdge(b)); }
        Extent(wr[m], r1, r2); Extent(wc[m], c1, c2);
      }
      if (wp[k].logx) x1 = std::max(x1, 1.);
      double ylo = 0, yhi = 30, zlo = 0, zhi = 30;
      OpenRange(ylo, yhi, r1, r2, 0.62, 0.); OpenRange(zlo, zhi, c1, c2, 0.95, 0.);
      TH1D *up = tdrHist(Form("fwu%d", k), "rms [GeV]", ylo, yhi, wp[k].xlab, x1, x2);
      TH1D *dw = tdrHist(Form("fwd%d", k), "#sigma_{core} [GeV]", zlo, zhi, wp[k].xlab, x1, x2);
      up->GetXaxis()->SetLabelSize(0.);
      std::unique_ptr<TCanvas> cv(tdrDiCanvas(Form("c_w%d", k), up, dw, 8, 11));
      cv->cd(1); if (wp[k].logx) gPad->SetLogx();
      TLegend *leg = tdrLeg(0.42, 0.58, 0.94, 0.84);
      leg->SetTextSize(0.033); leg->SetNColumns(2); leg->SetHeader(Form("%s, all vertices", wp[k].what));
      TF1 *model = 0;
      if (wp[k].model) {
        model = new TF1(Form("model%d", k), "sqrt([0]*[0] + [1]*[1]*x)", std::max(x1, 0.), x2);
        model->SetParameters(sigma0, sigmaK); model->SetLineColor(kBlack); model->SetLineStyle(kDashed); model->SetLineWidth(2);
      }
      int iw = 0;
      for (int m = 0; m < kNM; ++m) if (wr[m]) {
        TGraphErrors *g = Graph(wr[m], wp[k].logx ? 0.03*(iw - kMidM) : 0.08*(iw - kMidM), wp[k].logx);
        if (g->GetN()) { tdrDraw(g, "Pz", MMarker[m], MColour[m], kSolid, -1, kNone, 0, 0.8); leg->AddEntry(g, MShort(m), "PL"); }
        ++iw;
      }
      if (model) { model->Draw("SAME"); leg->AddEntry(model, "#sqrt{#sigma_{0}^{2} + #sigma_{K}^{2} S_{v}}", "L"); }
      ThresholdsDi();
      cv->cd(2); if (wp[k].logx) gPad->SetLogx();
      iw = 0;
      for (int m = 0; m < kNM; ++m) if (wr[m]) {
        TGraphErrors *g = Graph(wc[m], wp[k].logx ? 0.03*(iw - kMidM) : 0.08*(iw - kMidM), wp[k].logx);
        if (g->GetN()) tdrDraw(g, "Pz", MMarker[m], MColour[m], kSolid, -1, kNone, 0, 0.8);
        ++iw;
      }
      if (model) model->Draw("SAME");
      fixOverlay();
      wfigs.push_back(SavePdf(cv.get(), wp[k].pdf));
      // the numbers: rms and core per x bin, for the printout
      printf("--- %s: rms / core sigma [GeV] per bin of %s\n  %-10s", wp[k].hist, wp[k].xlab, "bin");
      for (int m = 0; m < kNM; ++m) if (wr[m]) printf(" %13s", MTag(m));
      printf("\n");
      TH1D *anyw = 0; for (int m = 0; m < kNM; ++m) if (wr[m]) { anyw = wr[m]; break; }
      for (int b = 1; b <= anyw->GetNbinsX(); ++b) {
        bool anyb = false; for (int m = 0; m < kNM; ++m) anyb = anyb || Has(wr[m], b);
        if (!anyb) continue;
        printf("  %-10s", BinLabel(anyw->GetXaxis(), b).Data());
        for (int m = 0; m < kNM; ++m) if (wr[m]) printf(" %6s/%-6s", Val(wr[m], b, "%.1f", "-").Data(), Val(wc[m], b, "%.1f", "-").Data());
        printf("\n");
      }
      if (k == 2) {
        // the model against the oracle at S_v = 30 GeV: the residual of opart
        // is the resolution the noise model should describe
        const int b = anyw->GetXaxis()->FindFixBin(30.);
        capw = Form(" At $S_v$ in %s\\,GeV the model gives %.1f\\,GeV, the residual core of \\texttt{opart} is %s and of "
                    "\\texttt{jva} %s.", BinTex(anyw->GetXaxis(), b).Data(), sqrt(sigma0*sigma0 + sigmaK*sigmaK*anyw->GetBinCenter(b)),
                    wc[kOPart] ? Val(wc[kOPart], b, "%.1f\\,GeV", "not defined").Data() : "not defined",
                    wc[kJva] ? Val(wc[kJva], b, "%.1f\\,GeV", "not defined").Data() : "not defined");
      }
    }
    TexFigure(wfigs, "Width of the missing transverse momentum of all vertices: the components against $N_{\\mathrm{PU}}$ "
              "and against $S_v$, and the residual against the owner interaction's true imbalance against $S_v$ "
              "(owned vertices), with the noise model of the optimiser, "
              "$\\sqrt{\\sigma_0^2+\\sigma_K^2 S_v}$ (dashed). Upper panels the rms, lower panels the core $\\sigma$; methods "
              "as in @fig:met_dist@. The model describes the vertex-resolved part against the owner's true central recoil; the residual "
              "here is against the true imbalance over $|\\eta|<5$ and so also holds the part of the owner's forward recoil that no "
              "forward cluster carries (PUPPI keeps little of it), which is why it lies above the model at large $S_v$." + capw, "met_width");
  }

  // ======================================================================
  // 2. the forward spike and the purity of the jets
  // ======================================================================
  printf("\n=== 2. forward spike: jets per owned vertex, spectra, own fraction, dominant share ===\n");
  const double normV = nown > 0 ? nown : (nvtx > 0 ? nvtx : 1.);
  const char *perWhat = nown > 0 ? "owned vertex" : (nvtx > 0 ? "vertex" : "file");
  if (nown <= 0) printf("drawJVA: WARNING hcount 'owned' is missing or 0; jets are normalised per %s\n", perWhat);
  TexSection("The forward spike", Form("Jets of the per-vertex events of the owned vertices under each method, against the pure "
             "generated jets of the owner interaction (filled once per owned vertex): $\\mathrm{d}N/\\mathrm{d}|\\eta|$ per owned "
             "vertex above three thresholds in $\\pt^{\\mathrm{raw}}$ (@fig:spike@), the %s "
             "spectra in three regions (@fig:spectra@), and the purity of the jets above 5\\,GeV: the fraction "
             "of a jet's linked generated $\\pt$ that comes from the owner interaction, and the share of its dominant "
             "interaction (@fig:purity@). A jet made of several interactions has a dominant share well below 1; "
             "the spike at $2.5<|\\eta|<3.0$ is made of them.", gJecOn ? "$\\pt^{\\mathrm{corr}}$" : "$\\pt^{\\mathrm{raw}}$"));
  {
    // dN/d|eta| per owned vertex, three thresholds
    std::vector<TString> figs; TString cap;
    double nr1[kNM+1][kNThr]; bool hr1[kNM+1][kNThr];     // jets per owned vertex at 2.5 < |eta| < 3.0; index kNM = gen
    memset(nr1, 0, sizeof(nr1)); memset(hr1, 0, sizeof(hr1));
    for (int t = 0; t < kNThr; ++t) {
      TH1 *g = GetAny({Form("jets/gen/hjeteta_pt%d", kThr[t]), Form("jets/gen/hjeteta_pt%02d", kThr[t])});
      TH1 *h[kNM] = {0}; int nh = 0;
      for (int m = 0; m < kNM; ++m) {
        h[m] = GetAny({Form("jets/%s/hjeteta_pt%d", MTag(m), kThr[t]), Form("jets/%s/hjeteta_pt%02d", MTag(m), kThr[t])});
        if (h[m]) { ++nh; hr1[m][t] = true; nr1[m][t] = Sum(h[m], 2.5, 3.0)/normV; }
      }
      if (g) { hr1[kNM][t] = true; nr1[kNM][t] = Sum(g, 2.5, 3.0)/normV; }
      if (!nh && !g) continue;
      TH1 *ref = g ? g : h[0]; for (int m = 0; !ref && m < kNM; ++m) ref = h[m];
      TH1D *dg = g ? PerUnit(g, Form("dsg%d", t), normV) : 0, *dn[kNM] = {0}, *rt[kNM] = {0};
      double ymax = dg ? dg->GetMaximum() : 0, ymin = 1e30;
      for (int m = 0; m < kNM; ++m) if (h[m]) {
        dn[m] = PerUnit(h[m], Form("dsr%d%d", t, m), normV); ymax = std::max(ymax, dn[m]->GetMaximum());
        for (int b = 1; b <= dn[m]->GetNbinsX(); ++b) if (dn[m]->GetBinContent(b) > 0) ymin = std::min(ymin, dn[m]->GetBinContent(b));
        if (g) rt[m] = Ratio(h[m], g, Form("dsq%d%d", t, m));
      }
      if (dg) for (int b = 1; b <= dg->GetNbinsX(); ++b) if (dg->GetBinContent(b) > 0) ymin = std::min(ymin, dg->GetBinContent(b));
      if (ymax <= 0) continue;
      if (ymin > ymax) ymin = ymax*1e-3;
      const double xhi = std::min(5., ref->GetXaxis()->GetXmax());
      const double flo = ymin*0.3, fhi = flo*pow(ymax/flo, 1/0.55);
      double rlo, rhi; const bool rlog = RatioFrame(rt, kNM, xhi, rlo, rhi);
      TH1D *up = tdrHist(Form("fsu%d", t), Form("Jets / |#eta| per %s", perWhat), flo, fhi, "|#eta|", 0., xhi);
      TH1D *dw = tdrHist(Form("fsd%d", t), "Reco / gen", rlo, rhi, "|#eta|", 0., xhi);
      up->GetXaxis()->SetLabelSize(0.);
      std::unique_ptr<TCanvas> cv(tdrDiCanvas(Form("c_sp%d", t), up, dw, 8, 11));
      cv->cd(1); gPad->SetLogy();
      TLegend *leg = tdrLeg(0.42, 0.52, 0.94, 0.82);
      leg->SetTextSize(0.033); leg->SetNColumns(2); leg->SetHeader(Form("p_{T}^{raw} > %d GeV (gen: p_{T} > %d GeV)", kThr[t], kThr[t]));
      if (dg) { tdrDraw(dg, "HIST", kNone, kBlack, kSolid, -1, kNone, 0, 0, 2); leg->AddEntry(dg, "pure gen, owner", "L"); }
      for (int m = 0; m < kNM; ++m) if (dn[m]) { tdrDraw(dn[m], "Pz", MMarker[m], MColour[m], kSolid, -1, kNone, 0, 0.7); leg->AddEntry(dn[m], MShort(m), "PL"); }
      ThresholdsDi();
      cv->cd(2); if (rlog) gPad->SetLogy();
      HLine(0, xhi);
      { TLine l; l.SetLineStyle(kDotted); l.SetLineColor(kGray+1); l.DrawLine(2.5, rlo, 2.5, rhi); l.DrawLine(3.0, rlo, 3.0, rhi); }
      for (int m = 0; m < kNM; ++m) if (rt[m]) tdrDraw(rt[m], "Pz", MMarker[m], MColour[m], kSolid, -1, kNone, 0, 0.7);
      fixOverlay();
      figs.push_back(SavePdf(cv.get(), Form("jva_spike_pt%02d", kThr[t])));
    }
    // the table: jets per owned vertex at 2.5 < |eta| < 3.0, reco / gen
    for (int m = 0; m < kNM; ++m) if (hr1[m][1] && hr1[kNM][1] && nr1[kNM][1] > 0) { haveSpike[m] = true; spike10[m] = nr1[m][1]/nr1[kNM][1]; }

    // pT spectra in the three regions
    std::vector<TString> sfigs; TString caps;
    for (int r = 0; r < kNReg; ++r) {
      TH1 *g = GetT<TH1>(Form("jets/gen/hjetpt_r%d", r));
      TH1 *h[kNM] = {0}; int nh = 0;
      for (int m = 0; m < kNM; ++m) if ((h[m] = GetT<TH1>(Form("jets/%s/hjetpt_r%d", MTag(m), r)))) ++nh;
      if (!nh) continue;
      TH1 *ref = g; for (int m = 0; !ref && m < kNM; ++m) ref = h[m];
      TH1D *dg = g ? PerUnit(g, Form("ptg%d", r), normV) : 0, *dn[kNM] = {0}, *rt[kNM] = {0};
      double ymax = dg ? dg->GetMaximum() : 0;
      for (int m = 0; m < kNM; ++m) if (h[m]) {
        dn[m] = PerUnit(h[m], Form("ptr%d%d", r, m), normV); ymax = std::max(ymax, dn[m]->GetMaximum());
        if (g) rt[m] = Ratio(h[m], g, Form("ptq%d%d", r, m));
      }
      if (ymax <= 0) continue;
      const double xlo = AxLo(ref);
      double rlo, rhi; const bool rlog = RatioFrame(rt, kNM, kPtDrawHi, rlo, rhi);
      TH1D *up = tdrHist(Form("fpu%d", r), "dN / dp_{T} [GeV^{-1}]", ymax*1e-7, ymax*1e-7*pow(1e7, 1/0.6), "p_{T} [GeV]", xlo, kPtDrawHi);
      TH1D *dw = tdrHist(Form("fpd%d", r), "Reco / gen", rlo, rhi, Form("%s (reco), p_{T} (gen) [GeV]", PtCorr()), xlo, kPtDrawHi);
      up->GetXaxis()->SetLabelSize(0.);
      std::unique_ptr<TCanvas> cv(tdrDiCanvas(Form("c_pt%d", r), up, dw, 8, 11));
      cv->cd(1); gPad->SetLogx(); gPad->SetLogy();
      TLegend *leg = tdrLeg(0.50, 0.52, 0.94, 0.82);
      leg->SetTextSize(0.033); leg->SetNColumns(2); leg->SetHeader(Form("%s, per %s", RegLabel(r), perWhat));
      if (dg) { tdrDraw(dg, "HIST", kNone, kBlack, kSolid, -1, kNone, 0, 0, 2); leg->AddEntry(dg, "pure gen, owner", "L"); }
      for (int m = 0; m < kNM; ++m) if (dn[m]) { tdrDraw(dn[m], "Pz", MMarker[m], MColour[m], kSolid, -1, kNone, 0, 0.7); leg->AddEntry(dn[m], MShort(m), "PL"); }
      ThresholdsDi();
      cv->cd(2); gPad->SetLogx(); if (rlog) gPad->SetLogy();
      HLine(xlo, kPtDrawHi);
      int iw = 0;
      for (int m = 0; m < kNM; ++m) if (rt[m]) {
        TGraphErrors *gr = Graph(rt[m], 0.025*(iw++ - kMidM), true, kPtDrawHi);
        if (gr->GetN()) tdrDraw(gr, "Pz", MMarker[m], MColour[m], kSolid, -1, kNone, 0, 0.7);
      }
      fixOverlay();
      sfigs.push_back(SavePdf(cv.get(), Form("jva_spectra_r%d", r)));
      if (r == 1 && g) {
        printf("--- 2.5 < |eta| < 3.0: reco/gen of the %s spectrum at 10-12 / 15-18 / 21-24 GeV\n", gJecOn ? "pT^corr" : "pT^raw");
        for (int m = 0; m < kNM; ++m) if (rt[m])
          printf("  %-6s %6s %6s %6s\n", MTag(m), Val(rt[m], BinAt(rt[m], 10), "%.2f", "-").Data(),
                 Val(rt[m], BinAt(rt[m], 15), "%.2f", "-").Data(), Val(rt[m], BinAt(rt[m], 21), "%.2f", "-").Data());
        caps = " At $2.5<|\\eta|<3.0$ and 10--12\\,GeV reco/gen is ";
        int nc = 0;
        for (int m = 0; m < kNM; ++m) if (rt[m]) caps += Form("%s%s (%s)", nc++ ? ", " : "", Val(rt[m], BinAt(rt[m], 10), "%.2f", "not defined").Data(), MTex(m));
        caps += ".";
      }
    }

    // purity: own fraction and dominant share
    std::vector<TString> pfigs; TString capp;
    double dlo[kNM] = {0}, dhi[kNM] = {0}; bool hdm[kNM] = {false};
    {
      TProfile *p[kNM] = {0}; TH1D *pr[kNM] = {0}; int np = 0;
      for (int m = 0; m < kNM; ++m) if ((p[m] = GetP(Form("jets/%s/pown_eta", MTag(m))))) {
        pr[m] = Like(p[m], Form("own%d", m));
        for (int b = 1; b <= p[m]->GetNbinsX(); ++b) if (HasP(p[m], b, 10)) {
          pr[m]->SetBinContent(b, p[m]->GetBinContent(b)); pr[m]->SetBinError(b, std::max(p[m]->GetBinError(b), 1e-4)); }
        // the mean over 2.5 < |eta| < 3.0, weighted by the jets of each bin
        double sw = 0, sv = 0;
        for (int b = p[m]->GetXaxis()->FindFixBin(2.5 + 1e-6); b <= p[m]->GetXaxis()->FindFixBin(3.0 - 1e-6); ++b) {
          sw += p[m]->GetBinEntries(b); sv += p[m]->GetBinEntries(b)*p[m]->GetBinContent(b); }
        if (sw > 0) { haveOwn[m] = true; own1[m] = sv/sw; }
        ++np;
      }
      if (np) {
        TProfile *ref = 0; for (int m = 0; m < kNM; ++m) if (p[m]) { ref = p[m]; break; }
        const double xhi = std::min(5., ref->GetXaxis()->GetXmax());
        TH1D *frame = tdrHist("fown", "Own fraction of the linked p_{T}", 0., 1.49, "|#eta|", 0., xhi);
        std::unique_ptr<TCanvas> cv(tdrCanvas("c_own", frame, 8, 11, kSquare));
        HLine(0, xhi);
        { TLine l; l.SetLineStyle(kDotted); l.SetLineColor(kGray+1); l.DrawLine(2.5, 0, 2.5, 1.0); l.DrawLine(3.0, 0, 3.0, 1.0); }
        TLegend *leg = tdrLeg(0.40, 0.70, 0.94, 0.86);
        leg->SetTextSize(0.028); leg->SetNColumns(2); leg->SetHeader("jets p_{T}^{raw} > 5 GeV, owned vertices");
        for (int m = 0; m < kNM; ++m) if (pr[m] && HasPoints(pr[m])) {
          tdrDraw(pr[m], "Pz", MMarker[m], MColour[m], kSolid, -1, kNone, 0, 0.7); leg->AddEntry(pr[m], MShort(m), "PL"); }
        Thresholds();
        fixOverlay();
        pfigs.push_back(SavePdf(cv.get(), "jva_own"));
      }
    }
    for (int r = 0; r < kNReg; ++r) {
      TH1D *u[kNM] = {0}; int nu = 0; double ymax = 0;
      for (int m = 0; m < kNM; ++m) {
        TH1 *h = GetT<TH1>(Form("jets/%s/hdom_r%d", MTag(m), r));
        if (!h || h->Integral(0, h->GetNbinsX()+1) <= 0) continue;
        u[m] = Unit(h, Form("dom%d%d", r, m), 1, true); ymax = std::max(ymax, u[m]->GetMaximum()); ++nu;
        if (r == 1) {
          // below 0.4 and above 0.6, the edges hit as bin edges where they are
          hdm[m] = true;
          const TAxis *a = u[m]->GetXaxis();
          dlo[m] = u[m]->Integral(1, a->FindFixBin(0.4 - 1e-6));
          dhi[m] = u[m]->Integral(a->FindFixBin(0.6 + 1e-6), u[m]->GetNbinsX());
        }
      }
      if (!nu) continue;
      TH1D *frame = tdrHist(Form("fdom%d", r), "Fraction of jets", 0., ymax*1.9, "Dominant share of the linked p_{T}", 0., 1.);
      std::unique_ptr<TCanvas> cv(tdrCanvas(Form("c_dom%d", r), frame, 8, 11, kSquare));
      TLegend *leg = tdrLeg(0.20, 0.54, 0.70, 0.76);
      leg->SetTextSize(0.028); leg->SetNColumns(2); leg->SetHeader(Form("%s, p_{T}^{raw} > 5 GeV", RegLabel(r)));
      for (int m = 0; m < kNM; ++m) if (u[m]) { tdrDraw(u[m], "HIST", kNone, MColour[m], MLine[m], -1, kNone, 0, 0, m == kJva ? 3 : 2); leg->AddEntry(u[m], MShort(m), "L"); }
      Thresholds();
      fixOverlay();
      pfigs.push_back(SavePdf(cv.get(), Form("jva_dom_r%d", r)));
    }
    // the spike table
    TexTableBegin(Form("The forward spike: jets per %s at $2.5<|\\eta|<3.0$ above 5, 10 and 20\\,GeV (reconstructed: "
                       "$\\pt^{\\mathrm{raw}}$) and their ratio to the pure generated jets of the owner interaction; the mean "
                       "own-interaction fraction of the linked $\\pt$ there, and the fraction of those jets whose dominant "
                       "interaction has less than 0.4, and more than 0.6, of the linked $\\pt$ ($\\pt^{\\mathrm{raw}}>5$\\,GeV).", perWhat),
                  "spike", "lrrrrrrrrr",
                  TString("method & \\multicolumn{3}{c}{jets per ") + perWhat + ", $\\pt>$} & \\multicolumn{3}{c}{reco/gen, $\\pt>$} & own & "
                  "\\multicolumn{2}{c}{dominant share} \\\\\n & 5 & 10 & 20 & 5 & 10 & 20 & & $<0.4$ & $>0.6$",
                  {"N>5", "N>10", "N>20", "r>5", "r>10", "r>20", "own", "dom<0.4", "dom>0.6"}, 7);
    for (int m = 0; m <= kNM; ++m) {
      bool any = false; for (int t = 0; t < kNThr; ++t) any = any || hr1[m][t];
      if (m < kNM) any = any || haveOwn[m] || hdm[m];
      if (!any) continue;
      std::vector<TString> cells;
      for (int t = 0; t < kNThr; ++t) cells.push_back(Num(hr1[m][t], nr1[m][t], "%.4f"));
      for (int t = 0; t < kNThr; ++t) cells.push_back(m < kNM ? Num(hr1[m][t] && hr1[kNM][t] && nr1[kNM][t] > 0, nr1[m][t]/std::max(nr1[kNM][t], 1e-30), "%.2f") : TString("1"));
      cells.push_back(m < kNM ? Num(haveOwn[m], own1[m], "%.3f") : TString("1"));
      cells.push_back(m < kNM ? Num(hdm[m], dlo[m], "%.3f") : TString("--"));
      cells.push_back(m < kNM ? Num(hdm[m], dhi[m], "%.3f") : TString("--"));
      TexRow(m < kNM ? TString(MTex(m)) : TString("pure gen"), m < kNM ? TString(MTag(m)) : TString("gen"), cells);
    }
    TexTableEnd();
    TString c1;
    for (int m = 0; m < kNM; ++m) if (haveSpike[m]) c1 += Form("%s%.2f (%s)", c1.IsNull() ? "" : ", ", spike10[m], MTex(m));
    if (!c1.IsNull()) cap = " At $2.5<|\\eta|<3.0$ above 10\\,GeV the reconstructed jets per owned vertex are " + c1 + " times the generated ones.";
    bool haveS[kNM]; for (int m = 0; m < kNM; ++m) haveS[m] = hr1[m][0] || hr1[m][1] || hr1[m][2];
    TexFigure(figs, Form("Jets per owned vertex against $|\\eta|$ above 5, 10 and 20\\,GeV (left to right; reconstructed jets in "
              "$\\pt^{\\mathrm{raw}}$): the pure generated jets of the owner interaction (black line) and the reconstructed jets of the "
              "owned vertex's event under %s, with the ratio to the generated jets below (on a logarithmic scale where the methods "
              "span more than a factor 20); the dotted lines mark $2.5<|\\eta|<3.0$.%s", MWho(haveS).Data(), cap.Data()), "spike");
    TexFigure(sfigs, TString("$\\pt$ spectra per owned vertex of the pure generated jets of the owner interaction and of the reconstructed "
              "jets (") + (gJecOn ? "$\\pt^{\\mathrm{corr}}$, corrected to the genseed link scale" : "$\\pt^{\\mathrm{raw}}$: no correction loaded") + ") under each method, at $|\\eta|<2.5$, $2.5<|\\eta|<3.0$ and $3.0<|\\eta|<5.0$ (left to "
              "right), with the ratio to the generated spectrum below (on a logarithmic scale where the methods span more than a "
              "factor 20)." + caps, "spectra");
    for (int m = 0; m < kNM; ++m) if (hdm[m])
      capp += Form("%s%s %.2f / %.2f", capp.IsNull() ? " Fraction of the jets at $2.5<|\\eta|<3.0$ with a dominant share below 0.4 / above 0.6: " : ", ",
                   MTex(m), dlo[m], dhi[m]);
    if (!capp.IsNull()) capp += ".";
    TString capo;
    for (int m = 0; m < kNM; ++m) if (haveOwn[m]) capo += Form("%s%s %.2f", capo.IsNull() ? " The mean own fraction at $2.5<|\\eta|<3.0$ is " : ", ", MTex(m), own1[m]);
    if (!capo.IsNull()) capo += ".";
    TexFigure(pfigs, "Purity of the jets above 5\\,GeV at the owned vertices: the fraction of the linked generated $\\pt$ that comes "
              "from the owner interaction against $|\\eta|$ (the profile), and the distribution of the dominant interaction's share of "
              "the linked $\\pt$ at $|\\eta|<2.5$, $2.5<|\\eta|<3.0$ and $3.0<|\\eta|<5.0$ (a share of exactly 1 is in the last bin)." +
              capo + capp, "purity", 0.48);
  }

  // ======================================================================
  // 3. the forward-cluster assignment
  // ======================================================================
  printf("\n=== 3. forward clusters: spectra, share, candidates, assignment outcomes ===\n");
  TexSection("Forward-cluster assignment", Form("Every cluster of the vertex-blind candidates ($|\\eta|>%g$, anti-$k_{\\mathrm{T}}$ "
             "$R=%g$) is classified by the dominant interaction of its linked generated $\\pt$: \\emph{single} when that share "
             "exceeds 0.5 and the interaction has a usable vertex within %g\\,cm, \\emph{combination} otherwise (unlinked "
             "clusters included). Its fate under each method that assigns whole clusters is then \\emph{right} (a single "
             "cluster sent to the vertex of its interaction), \\emph{wrong} (sent elsewhere), \\emph{single nulled}, "
             "\\emph{combination nulled} (the right thing to do with a combination) or \\emph{combination kept}; "
             "\\emph{correct} is right plus combination nulled. The clusters with track candidates (charged candidates of a "
             "vertex above %g\\,GeV in the cone at $|\\eta|<%g$) and those without are shown apart: without tracks the "
             "optimiser can place a cluster at any vertex, at the price $\\tau_{\\mathrm{all}}=%g$.",
             etaFwd, R, dzOracle, trkMinPt, etaFwd, tauAll));
  {
    std::vector<TString> figs; TString cap;
    TH1 *hcl[2] = {GetT<TH1>("fwd/hclpt_trk"), GetT<TH1>("fwd/hclpt_notrk")};
    // spectra with and without track candidates, and the fraction with
    if (hcl[0] || hcl[1]) {
      TH1 *ref = hcl[0] ? hcl[0] : hcl[1];
      TH1D *d[2] = {0}; double ymax = 0;
      for (int t = 0; t < 2; ++t) if (hcl[t]) { d[t] = PerUnit(hcl[t], Form("clpt%d", t), events > 0 ? events : 1.); ymax = std::max(ymax, d[t]->GetMaximum()); }
      TH1D *fr = 0;
      if (hcl[0] && hcl[1]) { TH1D *all = (TH1D*)AsD(hcl[0], "clall0")->Clone("clall"); all->SetDirectory(0); all->Add(hcl[1]); fr = Fraction(hcl[0], all, "clfr"); }
      if (ymax > 0) {
        const double xlo = AxLo(ref);
        TH1D *up = tdrHist("fclu", "Clusters / GeV per crossing", ymax*1e-6, ymax*50, "cluster p_{T} [GeV]", xlo, kPtDrawHi);
        TH1D *dw = tdrHist("fcld", "With tracks", 0., 1.19, "cluster p_{T} [GeV]", xlo, kPtDrawHi);
        up->GetXaxis()->SetLabelSize(0.);
        std::unique_ptr<TCanvas> cv(tdrDiCanvas("c_clpt", up, dw, 8, 11));
        cv->cd(1); gPad->SetLogx(); gPad->SetLogy();
        TLegend *leg = tdrLeg(0.45, 0.62, 0.94, 0.82);
        leg->SetTextSize(0.033); leg->SetHeader(Form("forward clusters, |#eta| > %g", etaFwd));
        for (int t = 0; t < 2; ++t) if (d[t]) { tdrDraw(d[t], "Pz", t == 0 ? kFullCircle : kOpenSquare, t == 0 ? kBlack : kRed+1, kSolid, -1, kNone, 0, 0.8); leg->AddEntry(d[t], SelLabel(t), "PL"); }
        ThresholdsDi();
        cv->cd(2); gPad->SetLogx();
        HLine(xlo, kPtDrawHi);
        if (fr) tdrDraw(fr, "Pz", kFullCircle, kBlack, kSolid, -1, kNone, 0, 0.8);
        fixOverlay();
        figs.push_back(SavePdf(cv.get(), "jva_fwd_clpt"));
        printf("--- clusters per crossing with / without track candidates, and the fraction with, above 2 / 5 / 10 GeV\n");
        for (double x : {2., 5., 10.}) {
          const double a = hcl[0] ? hcl[0]->Integral(BinAt(hcl[0], x), hcl[0]->GetNbinsX()+1) : 0;
          const double b = hcl[1] ? hcl[1]->Integral(BinAt(hcl[1], x), hcl[1]->GetNbinsX()+1) : 0;
          printf("  > %4.0f GeV: %8.3f %8.3f  fraction %.3f\n", x, a*perX, b*perX, a + b > 0 ? a/(a + b) : 0.);
          if (x == 5.) cap = Form("Above 5\\,GeV there are %.2f clusters per crossing, %.0f\\%% of them with track candidates.", (a + b)*perX, a + b > 0 ? 100*a/(a + b) : 0.);
        }
      }
    }
    // the dominant share against cluster pT
    if (TH2 *hs = Get2("fwd/hshare")) {
      TH1D *all = ProjX(hs, "shall");
      const int ny = hs->GetNbinsY();
      TH1D *n5 = Like(hs, "sh5"), *n6 = Like(hs, "sh6"), *n4 = Like(hs, "sh4");
      for (int b = 0; b <= hs->GetNbinsX()+1; ++b) {
        double a5 = 0, a6 = 0, a4 = 0;
        for (int j = 0; j <= ny+1; ++j) {
          const double y = j == 0 ? -1 : j == ny+1 ? 2 : hs->GetYaxis()->GetBinCenter(j), v = hs->GetBinContent(b, j);
          if (y > 0.5) a5 += v;
          if (y > 0.6) a6 += v;
          if (y < 0.4) a4 += v;
        }
        n5->SetBinContent(b, a5); n6->SetBinContent(b, a6); n4->SetBinContent(b, a4);
      }
      TH1D *f5 = Fraction(n5, all, "shf5"), *f6 = Fraction(n6, all, "shf6"), *f4 = Fraction(n4, all, "shf4");
      const double xlo = AxLo(all);
      TH1D *frame = tdrHist("fsh", "Fraction of forward clusters", 0., 1.39, "cluster p_{T} [GeV]", xlo, kPtDrawHi);
      std::unique_ptr<TCanvas> cv(tdrCanvas("c_sh", frame, 8, 11, kSquare));
      cv->SetLogx();
      HLine(xlo, kPtDrawHi);
      TLegend *leg = tdrLeg(0.40, 0.66, 0.94, 0.86);
      leg->SetTextSize(0.028); leg->SetHeader("dominant-interaction share of the linked p_{T}");
      tdrDraw(f5, "Pz", kFullCircle, kBlack, kSolid, -1, kNone, 0, 0.8);        leg->AddEntry(f5, "share > 0.5 (single)", "PL");
      tdrDraw(f6, "Pz", kFullSquare, kBlue+1, kSolid, -1, kNone, 0, 0.8);       leg->AddEntry(f6, "share > 0.6", "PL");
      tdrDraw(f4, "Pz", kFullTriangleUp, kRed+1, kSolid, -1, kNone, 0, 0.8);    leg->AddEntry(f4, "share < 0.4", "PL");
      Thresholds();
      fixOverlay();
      figs.push_back(SavePdf(cv.get(), "jva_fwd_share"));
      printf("--- fraction of forward clusters with share > 0.5 / > 0.6 / < 0.4 per cluster pT\n");
      for (int b = 1; b <= all->GetNbinsX(); ++b) if (Has(f5, b) && all->GetXaxis()->GetBinLowEdge(b) < 60)
        printf("  %-10s %9.0f %6.3f %6.3f %6.3f\n", BinLabel(all->GetXaxis(), b).Data(), all->GetBinContent(b), f5->GetBinContent(b), f6->GetBinContent(b), f4->GetBinContent(b));
      const int b10 = BinAt(all, 10);
      cap += Form(" At 10--12\\,GeV %s of the clusters are single (share $>0.5$) and %s have a share below 0.4.",
                  Val(f5, b10, "%.2f", "not defined").Data(), Val(f4, b10, "%.2f", "not defined").Data());
    }
    // the number of track candidates
    if (TH1 *hc = GetT<TH1>("fwd/hcand")) if (hc->Integral(0, hc->GetNbinsX()+1) > 0) {
      TH1D *u = Unit(hc, "ucand", 1, true);
      int last = 1; for (int b = 1; b <= hc->GetNbinsX(); ++b) if (hc->GetBinContent(b) > 0) last = b;
      const double cx2 = std::min(hc->GetXaxis()->GetXmax(), std::max(hc->GetXaxis()->GetXmin() + 10, hc->GetXaxis()->GetBinUpEdge(last) + 2));
      TH1D *frame = tdrHist("fcand", "Fraction of forward clusters", 0., u->GetMaximum()*1.6, "Number of track candidates",
                            hc->GetXaxis()->GetXmin(), cx2);
      std::unique_ptr<TCanvas> cv(tdrCanvas("c_cand", frame, 8, 11, kSquare));
      tdrDraw(u, "HIST", kNone, kBlack, kSolid, -1, kNone, 0, 0, 2);
      const double p0 = hc->Integral(0, hc->GetXaxis()->FindFixBin(0.5 - 1e-6)) / hc->Integral(0, hc->GetNbinsX()+1);
      Note(0.93, 0.80, Form("mean %.2f, none %.2f", hc->GetMean(), p0), 0.030, 31);
      Thresholds();
      fixOverlay();
      figs.push_back(SavePdf(cv.get(), "jva_fwd_ncand"));
      printf("--- track candidates per forward cluster: mean %.3f, fraction with none %.3f\n", hc->GetMean(), p0);
      cap += Form(" A cluster has %.2f track candidates on average; %.0f\\%% have none.", hc->GetMean(), 100*p0);
    }
    TexFigure(figs, "The forward clusters: their $\\pt$ spectrum per crossing with and without track candidates and the fraction with "
              "tracks, the fraction whose dominant interaction has more than 0.5, more than 0.6 and less than 0.4 of the linked "
              "$\\pt$ (an unlinked cluster has share 0), and the number of track candidates. " + cap, "fwd");

    // outcomes per method, with and without tracks
    std::vector<TString> ofigs, cfigs;
    double tot[kNA][2], frac[kNA][2][kNOut]; bool hof[kNA][2];
    memset(tot, 0, sizeof(tot)); memset(frac, 0, sizeof(frac)); memset(hof, 0, sizeof(hof));
    TH1D *corr[kNA][2] = {{0}};
    // The reference that "correct" has to beat: dropping every cluster nulls
    // every combination, so it scores the combination fraction.  It is the
    // none method's own outcomes (fwd/none, the dashed line) where the file
    // has them.  A file without them gets it derived, as before: single or
    // combination is a property of the cluster, not of the method, so the
    // reference follows from the first method with an hcat.
    TH1D *dropFr[2] = {0, 0}; double dropTot[2] = {0, 0}, dropFrac[2][kNOut];
    memset(dropFrac, 0, sizeof(dropFrac));
    bool dropNone[2] = {false, false};
    double worstId = 0; bool warnedLab = false, idChecked = false;
    auto outBins = [&](TH2 *h, int *ob) {
      for (int o = 0; o < kNOut; ++o) {
        ob[o] = h->GetYaxis()->FindFixBin(OutTag(o));
        if (ob[o] < 1 || ob[o] > h->GetNbinsY()) {
          ob[o] = o + 1;
          if (!warnedLab) { printf("drawJVA: WARNING %s has no outcome label '%s'; outcomes taken in bin order\n", h->GetName(), OutTag(o)); warnedLab = true; }
        }
      }
    };
    // the identity sum over outcomes = hclpt, bin by bin incl. under/overflow
    auto idCheck = [&](TH1D *den, int t) {
      if (!hcl[t]) return;
      idChecked = true;
      for (int b = 0; b <= den->GetNbinsX()+1; ++b) {
        const double x = hcl[t]->GetBinContent(b), y = den->GetBinContent(b);
        if (x > 0 || y > 0) worstId = std::max(worstId, fabs(x - y)/std::max(x, 1.));
      }
    };
    for (int t = 0; t < 2; ++t) {
      TH2 *h = Get2(Form("fwd/%s/hcat_%s", MTag(kNobody), SelTag(t)), true);
      if (!h) continue;
      int ob[kNOut]; outBins(h, ob);
      TH1D *den = ProjX(h, Form("catdenN%d", t));
      idCheck(den, t);
      const int b5 = BinAt(den, 5);
      dropTot[t] = den->Integral(b5, den->GetNbinsX()+1);
      if (!(dropTot[t] > 0)) { dropTot[t] = 0; continue; }
      TH1D *num[kNOut];
      for (int o = 0; o < kNOut; ++o) {
        num[o] = h->ProjectionX(Form("catnN%d%d", t, o), ob[o], ob[o]); num[o]->SetDirectory(0);
        dropFrac[t][o] = num[o]->Integral(b5, num[o]->GetNbinsX()+1)/dropTot[t];
      }
      TH1D *cn = (TH1D*)num[kRight]->Clone(Form("catcN%d", t)); cn->SetDirectory(0); cn->Add(num[kNullC]);
      dropFr[t] = Fraction(cn, den, Form("drop%d", t));
      dropNone[t] = true;
    }
    if (!dropNone[0] && !dropNone[1]) printf("--- no fwd/%s/hcat_*: the drop-all reference is derived from the combinations of the first method\n", MTag(kNobody));
    for (int t = 0; t < 2; ++t) for (int a = 0; a < kNA; ++a) {
      const int m = AMethod[a];
      TH2 *h = Get2(Form("fwd/%s/hcat_%s", MTag(m), SelTag(t)));
      if (!h) continue;
      int ob[kNOut]; outBins(h, ob);
      TH1D *den = ProjX(h, Form("catden%d%d", a, t));
      idCheck(den, t);
      TH1D *num[kNOut];
      for (int o = 0; o < kNOut; ++o) { num[o] = h->ProjectionX(Form("catn%d%d%d", a, t, o), ob[o], ob[o]); num[o]->SetDirectory(0); }
      TH1D *fr[kNOut];
      for (int o = 0; o < kNOut; ++o) fr[o] = Fraction(num[o], den, Form("catf%d%d%d", a, t, o));
      TH1D *cn = (TH1D*)num[kRight]->Clone(Form("catc%d%d", a, t)); cn->SetDirectory(0); cn->Add(num[kNullC]);
      corr[a][t] = Fraction(cn, den, Form("corr%d%d", a, t));
      // the table: clusters above 5 GeV
      const int b5 = BinAt(den, 5);
      tot[a][t] = den->Integral(b5, den->GetNbinsX()+1); hof[a][t] = tot[a][t] > 0;
      for (int o = 0; o < kNOut; ++o) frac[a][t][o] = hof[a][t] ? num[o]->Integral(b5, num[o]->GetNbinsX()+1)/tot[a][t] : 0;
      if (!dropFr[t] && hof[a][t]) {                  // the fallback: derived
        TH1D *cb = (TH1D*)num[kNullC]->Clone(Form("catcb%d", t)); cb->SetDirectory(0); cb->Add(num[kKeptC]);
        dropFr[t] = Fraction(cb, den, Form("drop%d", t));
        dropTot[t] = tot[a][t];
        dropFrac[t][kNullC] = frac[a][t][kNullC] + frac[a][t][kKeptC];
        dropFrac[t][kNullS] = frac[a][t][kRight] + frac[a][t][kWrong] + frac[a][t][kNullS];
      }
      const double xlo = AxLo(den);
      // the legend below the CMS label, the points below 1 = 0.56 of the frame
      TH1D *frame = tdrHist(Form("fcat%d%d", a, t), "Fraction of forward clusters", 0., 1.79, "cluster p_{T} [GeV]", xlo, kPtDrawHi);
      std::unique_ptr<TCanvas> cv(tdrCanvas(Form("c_cat%d%d", a, t), frame, 8, 11, kSquare));
      cv->SetLogx();
      HLine(xlo, kPtDrawHi);
      TLegend *leg = tdrLeg(0.20, 0.60, 0.94, 0.79);
      leg->SetTextSize(0.028); leg->SetNColumns(2); leg->SetHeader(Form("%s, %s", MTag(m), SelLabel(t)));
      for (int o = 0; o < kNOut; ++o) if (HasPoints(fr[o])) {
        tdrDraw(fr[o], "Pz", OutMarker[o], OutColour[o], kSolid, -1, kNone, 0, 0.8); leg->AddEntry(fr[o], OutLabel(o), "PL"); }
      Thresholds();
      fixOverlay();
      ofigs.push_back(SavePdf(cv.get(), Form("jva_cat_%s_%s", MTag(m), SelTag(t))));
    }
    if (idChecked) printf("--- identity sum over outcomes of hcat = hclpt: worst relative violation %.2e%s\n", worstId, worstId > 1e-6 ? "   IDENTITY BROKEN" : "");
    else printf("--- identity sum over outcomes of hcat = hclpt: not checked (no hcat with its hclpt)\n");
    // correct = right + combination nulled, per method
    for (int t = 0; t < 2; ++t) {
      int nc = 0; for (int a = 0; a < kNA; ++a) if (HasPoints(corr[a][t])) ++nc;
      if (!nc) continue;
      TH1D *any = 0; for (int a = 0; a < kNA; ++a) if (corr[a][t]) { any = corr[a][t]; break; }
      const double xlo = AxLo(any);
      TH1D *frame = tdrHist(Form("fcorr%d", t), "Correct fraction", 0., 1.59, "cluster p_{T} [GeV]", xlo, kPtDrawHi);
      std::unique_ptr<TCanvas> cv(tdrCanvas(Form("c_corr%d", t), frame, 8, 11, kSquare));
      cv->SetLogx();
      HLine(xlo, kPtDrawHi);
      TLegend *leg = tdrLeg(0.40, 0.66, 0.94, 0.86);
      leg->SetTextSize(0.028); leg->SetNColumns(2); leg->SetHeader(Form("forward clusters %s", SelLabel(t)));
      TGraphErrors *gd = Graph(dropFr[t], 0, true, kPtDrawHi);
      if (gd->GetN() > 1) { tdrDraw(gd, "LX", kNone, MColour[kNobody], MLine[kNobody], -1, kNone, 0, 0, 2);
                            leg->AddEntry(gd, dropNone[t] ? MShort(kNobody) : "drop all (derived)", "L"); }
      int iw = 0;
      for (int a = 0; a < kNA; ++a) if (HasPoints(corr[a][t])) {
        const int m = AMethod[a];
        TGraphErrors *g = Graph(corr[a][t], 0.03*(iw++ - 1.5), true, kPtDrawHi);
        if (g->GetN()) { tdrDraw(g, "Pz", MMarker[m], MColour[m], kSolid, -1, kNone, 0, 0.8); leg->AddEntry(g, MShort(m), "PL"); }
      }
      Thresholds();
      fixOverlay();
      cfigs.push_back(SavePdf(cv.get(), Form("jva_correct_%s", SelTag(t))));
    }
    // the table
    TexTableBegin("Fate of the forward clusters above 5\\,GeV under the methods that assign whole clusters, with and without track "
                  "candidates: the number of clusters and the fraction of each outcome; \\emph{correct} is right plus "
                  "combination nulled. The last rows are the reference of dropping every cluster (\\texttt{none}, or derived "
                  "from the combinations where the file has no \\texttt{none}), whose correct fraction is the fraction of "
                  "combinations.", "fwd", "llrrrrrrr",
                  "method & tracks & $N$ & right & wrong & single nulled & combo nulled & combo kept & correct",
                  {"N", "right", "wrong", "nulledS", "nulledC", "keptC", "correct"});
    for (int a = 0; a < kNA; ++a) {
      const int m = AMethod[a];
      if (!hof[a][0] && !hof[a][1]) continue;
      for (int t = 0; t < 2; ++t) {
        std::vector<TString> cells = {Num(hof[a][t], tot[a][t], "%.0f")};
        for (int o = 0; o < kNOut; ++o) cells.push_back(Num(hof[a][t], frac[a][t][o], "%.3f"));
        cells.push_back(Num(hof[a][t], frac[a][t][kRight] + frac[a][t][kNullC], "%.3f"));
        TexRow(Form("%s & %s", t == 0 ? MTex(m) : "", t == 0 ? "yes" : "no"), Form("%s/%s", MTag(m), SelTag(t)), cells);
      }
      if (a + 1 < kNA) TexHline();
      const double n = tot[a][0] + tot[a][1];
      if (n > 0) { haveCorrect[m] = true;
        correct[m] = (tot[a][0]*(frac[a][0][kRight] + frac[a][0][kNullC]) + tot[a][1]*(frac[a][1][kRight] + frac[a][1][kNullC]))/n; }
    }
    // the drop-all reference rows (see dropFr above): none as measured, or derived
    double &dropAll = dropAllRef;
    if (dropTot[0] > 0 || dropTot[1] > 0) {
      TexHline();
      const bool meas = dropNone[0] || dropNone[1];
      for (int t = 0; t < 2; ++t) {
        const bool h = dropTot[t] > 0;
        std::vector<TString> cells = {Num(h, dropTot[t], "%.0f")};
        for (int o = 0; o < kNOut; ++o) cells.push_back(Num(h, dropFrac[t][o], "%.3f"));
        cells.push_back(Num(h, dropFrac[t][kRight] + dropFrac[t][kNullC], "%.3f"));
        TexRow(Form("%s & %s", t == 0 ? (meas ? MTex(kNobody) : "drop all (ref.)") : "", t == 0 ? "yes" : "no"),
               Form("%s/%s", meas ? MTag(kNobody) : "dropall", SelTag(t)), cells);
      }
      const double n = dropTot[0] + dropTot[1];
      if (n > 0) dropAll = (dropTot[0]*(dropFrac[0][kRight] + dropFrac[0][kNullC]) + dropTot[1]*(dropFrac[1][kRight] + dropFrac[1][kNullC]))/n;
      if (n > 0 && meas) { haveCorrect[kNobody] = true; correct[kNobody] = dropAll; }
    }
    TexTableEnd();
    TString capo;
    for (int m = 0; m < kNM; ++m) if (haveCorrect[m]) capo += Form("%s%s %.2f", capo.IsNull() ? " Above 5\\,GeV the correct fraction is " : ", ", MTex(m), correct[m]);
    if (!capo.IsNull()) capo += " (clusters with and without tracks together)";
    if (!capo.IsNull() && dropAll >= 0 && !haveCorrect[kNobody]) capo += Form("; dropping every cluster would score %.2f, the fraction of combinations", dropAll);
    else if (!capo.IsNull() && haveCorrect[kNobody]) capo += Form(", where %s, dropping every cluster, scores the fraction of combinations", MTex(kNobody));
    if (!capo.IsNull()) capo += ".";
    if (idChecked) capo += Form(" The outcomes add up to the cluster spectrum in every bin (largest relative difference %.1g).", worstId);
    // input order t-major: {lv, trk, jva, ojet} with tracks, then without
    std::vector<TString> order; std::vector<int> rows(2, 0);
    for (int t = 0; t < 2; ++t) for (int a = 0; a < kNA; ++a) {
      order.push_back(Form("jva_cat_%s_%s", MTag(AMethod[a]), SelTag(t))); ++rows[t]; }
    TexFigure(order, "Outcome fractions of the forward clusters against the cluster $\\pt$ for \\texttt{lv}, \\texttt{trk}, "
              "\\texttt{jva} and \\texttt{ojet} (left to right), for the clusters with track candidates (first row) and without "
              "(second row): single cluster at the right vertex (green), at a wrong vertex (red), nulled (orange); combination "
              "nulled (blue) or kept (magenta)." + capo, "outcome", 0.24, rows);
    TString capc;
    for (int t = 0; t < 2; ++t) {
      TString c;
      for (int a = 0; a < kNA; ++a) if (hof[a][t]) c += Form("%s%s %.3f", c.IsNull() ? "" : ", ", MTex(AMethod[a]), frac[a][t][kRight] + frac[a][t][kNullC]);
      if (dropTot[t] > 0) c += Form("%s%s %.3f", c.IsNull() ? "" : "; ", dropNone[t] ? MTex(kNobody) : "dropping all",
                                    dropFrac[t][kRight] + dropFrac[t][kNullC]);
      if (!c.IsNull()) capc += Form(" Above 5\\,GeV %s tracks: %s.", t == 0 ? "with" : "without", c.Data());
    }
    TexFigure(cfigs, "The correct fraction (right vertex plus combination nulled) of the forward clusters against the cluster "
              "$\\pt$, with track candidates (left) and without (right), per method, and the reference of dropping every cluster "
              "(\\texttt{none}, black dashed), which scores the fraction of combinations: a method that places clusters helps only "
              "where it lies above that line." + capc, "correct");
  }

  // ======================================================================
  // 4. L2Res closure from every vertex
  // ======================================================================
  printf("\n=== 4. dijets: selection, and L2Res closure R_DB, R_MPF against r_true ===\n");
  TexSection("Relative response from every vertex", Form("Dijets are selected from the per-vertex event of every usable vertex: "
             "the two leading jets in %s above %g\\,GeV at $|\\eta|<5.191$, $\\Delta\\phi>%g$, a tag at $|\\eta|<1.305$ and the "
             "other jet as the probe (both assignments at weight 1/2 when both are in the barrel), $\\alpha=\\pt^{(3)}/\\pt^{\\mathrm{avg}}$. "
             "$R_{\\mathrm{DB}}=(1+\\langle A_{\\mathrm{DB}}\\rangle)/(1-\\langle A_{\\mathrm{DB}}\\rangle)$ with "
             "$A_{\\mathrm{DB}}=(\\pt^{\\mathrm{probe}}-\\pt^{\\mathrm{tag}})/(\\pt^{\\mathrm{probe}}+\\pt^{\\mathrm{tag}})$, and the same "
             "for $A_{\\mathrm{MPF}}=\\vec{p}_{\\mathrm{T}}^{\\,\\mathrm{miss,T1}}\\cdot\\hat{p}_{\\mathrm{T}}^{\\,\\mathrm{tag}}/(2\\pt^{\\mathrm{avg}})$ "
             "(type-1 correction above %g\\,GeV). The closure reference $r_{\\mathrm{true}}$ is the median of the per-event ratio of "
             "the probe's to the tag's $\\pt^{\\mathrm{corr}}/\\pt^{\\mathrm{linkz}}$, where $\\pt^{\\mathrm{linkz}}$ is the linked "
             "generated $\\pt$ of the jet. A mean or a median needs %.0f effective dijets.",
             gJecOn ? "$\\pt^{\\mathrm{corr}}$" : "$\\pt^{\\mathrm{raw}}$ (no correction loaded)", ptMinJet, dphiMin, t1Min, kMinMean));
  {
    // the selection counts per crossing
    TH1 *hsel[kNM][kNDC] = {{0}};
    TH1 *selRef = 0;
    for (int m = 0; m < kNM; ++m) for (int c = 0; c < kNDC; ++c)
      if ((hsel[m][c] = GetT<TH1>(Form("dijet/%s/%s/hsel", MTag(m), CTag(c)))) && !selRef) selRef = hsel[m][c];
    if (selRef) {
      const int ns = selRef->GetNbinsX();
      TString head = "method & vertices"; std::vector<TString> txt;
      for (int b = 1; b <= ns; ++b) {
        TString l = selRef->GetXaxis()->GetBinLabel(b), lt = l;
        if (l == ">=2 jets") lt = "$\\geq2$ jets"; else if (l == "alpha<0.3") lt = "$\\alpha<0.3$"; else if (l == "dphi") lt = "$\\Delta\\phi$";
        else if (l == "vertices") lt = "vertex events";
        else lt = TexEscape(l);
        head += " & " + lt; txt.push_back(l);
      }
      TexTableBegin("Dijet selection per crossing and method, for all vertices and for the primary vertex alone: the counts "
                    "after each step divided by the number of crossings, cumulative up to $\\alpha<0.3$ (which requires a tag); "
                    "\\emph{samebin} counts the pairs that pass $\\Delta\\phi$ and $\\alpha<0.3$ with both jets in the same "
                    "$|\\eta|$ bin, tag or not (the sample of the JER), so it can exceed the $\\alpha<0.3$ column.", "sel", Form("ll%s", TString('r', ns).Data()), head, txt);
      for (int m = 0; m < kNM; ++m) {
        if (!hsel[m][0] && !hsel[m][1]) continue;
        for (int c = 0; c < kNDC; ++c) {
          std::vector<TString> cells;
          for (int b = 1; b <= ns; ++b) cells.push_back(hsel[m][c] && b <= hsel[m][c]->GetNbinsX() ? TString(Form("%.4g", hsel[m][c]->GetBinContent(b)*perX)) : TString("--"));
          TexRow(Form("%s & %s", c == 0 ? MTex(m) : "", CTex(c)), Form("%s/%s", MTag(m), CTag(c)), cells);
          if (hsel[m][c]) { haveDij[m][c] = true; dijPerX[m][c] = hsel[m][c]->GetBinContent(ns)*perX; }
        }
        if (m + 1 < kNM) TexHline();
      }
      TexTableEnd();
    }

    // R_DB, R_MPF and r_true per method, class, alpha, pT bin, probe bin
    static Pt rr[kNM][kNDC][kNAl][kNL2][kNRV][kNP], rt[kNM][kNDC][kNAl][kNL2][kNP];
    memset(rr, 0, sizeof(rr)); memset(rt, 0, sizeof(rt));
    TString ptLab[kNL2], ptTex[kNL2]; bool ptHave[kNL2] = {false};
    bool haveM[kNM] = {false};
    for (int m = 0; m < kNM; ++m) for (int c = 0; c < kNDC; ++c) for (int a = 0; a < kNAl; ++a) for (int ie = 0; ie < kNP; ++ie) {
      const TString base = Form("dijet/%s/%s/", MTag(m), CTag(c)), suf = Form("_%s_a%d", PTag(ie).Data(), a);
      TH2 *hdb = Get2(base + "hadb" + suf), *hmp = Get2(base + "hampf" + suf), *htr = Get2(base + "hrtrue" + suf);
      if (!hdb && !hmp && !htr) continue;
      haveM[m] = true;
      for (int ip = 0; ip < kNL2; ++ip) {
        TH2 *ref = hdb ? hdb : hmp ? hmp : htr;
        const int b = BinAt(ref, kL2Pt[ip]);
        if (b < 1 || b > ref->GetNbinsX()) continue;
        if (!ptHave[ip]) { ptHave[ip] = true; ptLab[ip] = BinLabel(ref->GetXaxis(), b); ptTex[ip] = BinTex(ref->GetXaxis(), b); }
        rr[m][c][a][ip][kDB][ie] = RFromA(hdb, BinAt(hdb ? hdb : ref, kL2Pt[ip]));
        rr[m][c][a][ip][kMPF][ie] = RFromA(hmp, BinAt(hmp ? hmp : ref, kL2Pt[ip]));
        rt[m][c][a][ip][ie] = MedianBin(htr, BinAt(htr ? htr : ref, kL2Pt[ip]));
      }
    }
    // closure: mean |R/r_true - 1| over the probe bins where both exist
    double clo[kNM][kNDC][kNAl][kNL2][kNRV]; int cln[kNM][kNDC][kNAl][kNL2][kNRV];
    memset(clo, 0, sizeof(clo)); memset(cln, 0, sizeof(cln));
    for (int m = 0; m < kNM; ++m) for (int c = 0; c < kNDC; ++c) for (int a = 0; a < kNAl; ++a) for (int ip = 0; ip < kNL2; ++ip) for (int v = 0; v < kNRV; ++v) {
      double s = 0; int n = 0;
      for (int ie = 0; ie < kNP; ++ie) {
        const Pt &x = rr[m][c][a][ip][v][ie], &y = rt[m][c][a][ip][ie];
        if (x.ok && y.ok) { s += fabs(x.v/y.v - 1); ++n; }
      }
      cln[m][c][a][ip][v] = n; clo[m][c][a][ip][v] = n ? s/n : 0;
    }
    for (int m = 0; m < kNM; ++m) for (int c = 0; c < kNDC; ++c) {
      const int ip = 1;   // 8-10 GeV
      if (cln[m][c][0][ip][kMPF] > 0) { haveL2[m][c] = true; l2mpf[m][c] = clo[m][c][0][ip][kMPF]; l2n[m][c] = cln[m][c][0][ip][kMPF]; }
      if (cln[m][c][1][ip][kMPF] > 0) { haveL2a1[m][c] = true; l2mpfa1[m][c] = clo[m][c][1][ip][kMPF]; l2na1[m][c] = cln[m][c][1][ip][kMPF]; }
    }
    // the printout of every number of the figures (all vertices; both alpha
    // selections; the probe bins where some method has a number)
    for (int a = 0; a < kNAl; ++a) {
      TString none;
      for (int ip = 0; ip < kNL2; ++ip) {
        if (!ptHave[ip]) continue;
        bool anyp = false;
        for (int ie = 0; ie < kNP; ++ie) for (int m = 0; m < kNM; ++m)
          anyp = anyp || rr[m][0][a][ip][kDB][ie].ok || rr[m][0][a][ip][kMPF][ie].ok || rt[m][0][a][ip][ie].ok;
        if (!anyp) { none += Form(" %s", ptLab[ip].Data()); continue; }
        printf("--- all vertices, %s, pT_avg %s GeV: R_DB / R_MPF / r_true per probe |eta| bin\n  %-16s", a == 0 ? "alpha < 0.3" : "no alpha cut", ptLab[ip].Data(), "probe |eta|");
        for (int m = 0; m < kNM; ++m) if (haveM[m]) printf(" %20s", MTag(m));
        printf("\n");
        for (int ie = 0; ie < kNP; ++ie) {
          bool anye = false;
          for (int m = 0; m < kNM; ++m) anye = anye || rr[m][0][a][ip][kDB][ie].ok || rr[m][0][a][ip][kMPF][ie].ok || rt[m][0][a][ip][ie].ok;
          if (!anye) continue;
          printf("  %-16s", Form("%g-%g", kProbeEta[ie], kProbeEta[ie+1]));
          for (int m = 0; m < kNM; ++m) if (haveM[m])
            printf(" %6s/%6s/%6s", Num(rr[m][0][a][ip][kDB][ie].ok, rr[m][0][a][ip][kDB][ie].v, "%.3f", "-").Data(),
                   Num(rr[m][0][a][ip][kMPF][ie].ok, rr[m][0][a][ip][kMPF][ie].v, "%.3f", "-").Data(),
                   Num(rt[m][0][a][ip][ie].ok, rt[m][0][a][ip][ie].v, "%.3f", "-").Data());
          printf("\n");
        }
      }
      if (!none.IsNull()) printf("--- all vertices, %s: no probe bin with %.0f effective dijets at pT_avg%s GeV\n", a == 0 ? "alpha < 0.3" : "no alpha cut", kMinMean, none.Data());
    }
    // the figures: one per pT bin and estimator, all vertices; alpha < 0.3
    // (the L2Res selection) and, as a second set, no alpha cut: at low pT_avg
    // alpha < 0.3 means "no third jet above ptMinJet" and keeps 5-25% of the
    // tagged dijets, so the uncut set is the one with points on a small sample
    std::vector<TString> figs[kNAl][kNRV]; TString cap[kNAl][kNRV];
    std::vector<int> drawn[kNAl][kNRV];
    for (int a = 0; a < kNAl; ++a) for (int v = 0; v < kNRV; ++v) for (int ip = 0; ip < kNL2; ++ip) {
      if (!ptHave[ip]) continue;
      TGraphErrors *gR[kNM] = {0}, *gT[kNM] = {0}, *gQ[kNM] = {0};
      double r1 = 1e30, r2 = -1e30, q1 = 1e30, q2 = -1e30; int ng = 0;
      for (int m = 0; m < kNM; ++m) {
        if (!haveM[m]) continue;
        gR[m] = new TGraphErrors(); gT[m] = new TGraphErrors(); gQ[m] = new TGraphErrors();
        const double sh = 0.09*(m - 2.5);
        for (int ie = 0; ie < kNP; ++ie) {
          const double xc = 0.5*(kProbeEta[ie] + kProbeEta[ie+1]), w = kProbeEta[ie+1] - kProbeEta[ie], x = xc + sh*w;
          const Pt &p = rr[m][0][a][ip][v][ie], &t = rt[m][0][a][ip][ie];
          if (p.ok) { const int n = gR[m]->GetN(); gR[m]->SetPoint(n, x, p.v); gR[m]->SetPointError(n, 0, p.e);
                      if (p.e < 0.1) { r1 = std::min(r1, p.v - p.e); r2 = std::max(r2, p.v + p.e); } }
          if (t.ok) { const int n = gT[m]->GetN(); gT[m]->SetPoint(n, xc, t.v); r1 = std::min(r1, t.v); r2 = std::max(r2, t.v); }
          if (p.ok && t.ok) {
            const int n = gQ[m]->GetN();
            const double y = p.v/t.v, e = y*sqrt(pow(p.e/p.v, 2) + pow(t.e/t.v, 2));
            gQ[m]->SetPoint(n, x, y); gQ[m]->SetPointError(n, 0, e);
            if (e < 0.1) { q1 = std::min(q1, y - e); q2 = std::max(q2, y + e); }
          }
        }
        if (gR[m]->GetN() || gT[m]->GetN()) ++ng;
      }
      // a figure needs a curve: some method with R in two probe bins (a
      // lone point against |eta| is left to the printout)
      int maxR = 0; for (int m = 0; m < kNM; ++m) if (gR[m]) maxR = std::max(maxR, gR[m]->GetN());
      if (!ng || maxR < 2) continue;
      double ylo = 0.7, yhi = 1.3, zlo = 0.8, zhi = 1.2;
      OpenRange(ylo, yhi, r1, r2, 0.55, 0.); OpenRange(zlo, zhi, q1, q2, 0.95, 0.);
      if (q1 < zlo) zlo = std::max(0., q1 - 0.05);
      TH1D *up = tdrHist(Form("fl2u%d%d%d", a, v, ip), Form("%s, r_{true}", RVLabel(v)), ylo, yhi, "Probe |#eta|", 0., kProbeEta[kNP]);
      TH1D *dw = tdrHist(Form("fl2d%d%d%d", a, v, ip), Form("%s / r_{true}", RVLabel(v)), zlo, zhi, "Probe |#eta|", 0., kProbeEta[kNP]);
      up->GetXaxis()->SetLabelSize(0.);
      std::unique_ptr<TCanvas> cv(tdrDiCanvas(Form("c_l2%d%d%d", a, v, ip), up, dw, 8, 11));
      cv->cd(1);
      HLine(0, kProbeEta[kNP]);
      TLegend *leg = tdrLeg(0.42, 0.56, 0.94, 0.82);
      leg->SetTextSize(0.033); leg->SetNColumns(2);
      leg->SetHeader(Form("p_{T}^{avg} %s GeV, %s, all vertices", ptLab[ip].Data(), AlLabel(a)));
      for (int m = 0; m < kNM; ++m) {
        if (gT[m] && gT[m]->GetN() > 1) tdrDraw(gT[m], "LX", kNone, MColour[m], MLine[m] == kSolid ? kDashed : MLine[m], -1, kNone, 0, 0, 1);
        if (gR[m] && gR[m]->GetN()) { tdrDraw(gR[m], "Pz", MMarker[m], MColour[m], kSolid, -1, kNone, 0, 0.7); leg->AddEntry(gR[m], MShort(m), "PL"); }
      }
      leg->AddEntry((TObject*)0, Form("markers %s, lines r_{true}", RVLabel(v)), "");
      ThresholdsDi();
      cv->cd(2);
      HLine(0, kProbeEta[kNP]);
      HLine(0, kProbeEta[kNP], 0.98, kGray, kDotted); HLine(0, kProbeEta[kNP], 1.02, kGray, kDotted);
      for (int m = 0; m < kNM; ++m) if (gQ[m] && gQ[m]->GetN()) tdrDraw(gQ[m], "Pz", MMarker[m], MColour[m], kSolid, -1, kNone, 0, 0.7);
      fixOverlay();
      figs[a][v].push_back(SavePdf(cv.get(), Form("jva_l2_%s_pt%02.0f%s", RVTag(v), kL2Pt[ip], a == 0 ? "" : "_a1")));
      drawn[a][v].push_back(ip);
    }
    // the closure table: all vertices and PV, alpha < 0.3; rows pT bin x estimator
    TString mh; std::vector<TString> mtxt;
    for (int m = 0; m < kNM; ++m) if (haveM[m]) { mh += TString(" & ") + MTex(m); mtxt.push_back(MTag(m)); }
    int nmh = 0; for (int m = 0; m < kNM; ++m) if (haveM[m]) ++nmh;
    if (nmh) {
      // rows: vertex class x alpha selection x pT bin x estimator, only those
      // where some method has a number (a row of dashes says nothing)
      TexTableBegin(Form("L2Res closure: the mean of $|R/r_{\\mathrm{true}}-1|$ over the probe $|\\eta|$ bins in which both are "
                         "measured (their number in parentheses), per $\\alpha$ selection, $\\pt^{\\mathrm{avg}}$ bin, estimator and "
                         "method, for the dijets of all vertices and of the primary vertex alone; only rows with a number are listed. "
                         "A mean or median needs %.0f effective dijets.", kMinMean), "l2res", Form("llll%s", TString('r', nmh).Data()),
                    "vertices & $\\alpha$ & $\\pt^{\\mathrm{avg}}$ [GeV] & $R$" + mh, mtxt, 10,
                    // seven method columns are 15pt too wide for the text at footnotesize
                    nmh > 6 ? "\\scriptsize\\setlength{\\tabcolsep}{3pt}" : "\\footnotesize\\setlength{\\tabcolsep}{4pt}");
      for (int c = 0; c < kNDC; ++c) {
        bool firstC = true;
        for (int a = 0; a < kNAl; ++a) {
          bool firstA = true;
          for (int ip = 0; ip < kNL2; ++ip) {
            if (!ptHave[ip]) continue;
            bool firstP = true;
            for (int v = 0; v < kNRV; ++v) {
              bool anyr = false; std::vector<TString> cells;
              for (int m = 0; m < kNM; ++m) if (haveM[m]) {
                anyr = anyr || cln[m][c][a][ip][v] > 0;
                cells.push_back(cln[m][c][a][ip][v] > 0 ? TString(Form("%.3f (%d)", clo[m][c][a][ip][v], cln[m][c][a][ip][v])) : TString("--"));
              }
              if (!anyr) continue;
              if (firstC && gTabRows > 0) TexHline();
              TexRow(Form("%s & %s & %s & %s", firstC ? CTex(c) : "", firstA ? (a == 0 ? "$<0.3$" : "none") : "", firstP ? ptTex[ip].Data() : "", RVTex(v)),
                     Form("%s/%s/%s/%s", CTag(c), a == 0 ? "a<0.3" : "a:none", ptLab[ip].Data(), RVTag(v)), cells);
              firstC = firstA = firstP = false;
            }
          }
        }
      }
      TexTableEnd();
    }
    // captions: the pT_avg bins actually drawn, and the closure of every
    // method at each of them (the numbers of tab:l2res)
    for (int a = 0; a < kNAl; ++a) for (int v = 0; v < kNRV; ++v) {
      TString bins, nums;
      for (size_t k = 0; k < drawn[a][v].size(); ++k) {
        const int ip = drawn[a][v][k];
        bins += Form("%s%s", k == 0 ? "" : k + 1 == drawn[a][v].size() ? " and " : ", ", ptTex[ip].Data());
        TString who;
        for (int m = 0; m < kNM; ++m) if (cln[m][0][a][ip][v] > 0)
          who += Form("%s%s %.3f (%d)", who.IsNull() ? "" : ", ", MTex(m), clo[m][0][a][ip][v], cln[m][0][a][ip][v]);
        if (!who.IsNull()) nums += Form("%sat %s\\,GeV %s", nums.IsNull() ? "" : "; ", ptTex[ip].Data(), who.Data());
      }
      const TString paren = drawn[a][v].size() > 1
        ? Form("left to right, top to bottom; of the five L2Res bins, those in which some method has two probe bins with %.0f effective dijets", kMinMean)
        : Form("the only one of the five L2Res bins in which some method has two probe bins with %.0f effective dijets", kMinMean);
      if (!nums.IsNull()) cap[a][v] = Form(" The mean $|%s/r_{\\mathrm{true}}-1|$ over the probe bins where both are measured (their "
                                          "number in parentheses) is, ", v == kDB ? "R_{\\mathrm{DB}}" : "R_{\\mathrm{MPF}}") + nums + ".";
      TexFigure(figs[a][v], Form("%s against the probe $|\\eta|$ for the dijets of all vertices %s, at $\\pt^{\\mathrm{avg}}$ "
                "%s\\,GeV (%s), per method (markers, slightly "
                "shifted in $|\\eta|$ for legibility), with the truth relative response $r_{\\mathrm{true}}$ of the same "
                "dijets (lines, same colour); below, the closure $R/r_{\\mathrm{true}}$, with $\\pm2\\%%$ dotted.%s",
                RVTex(v), a == 0 ? "with $\\alpha<0.3$" : "without the $\\alpha$ cut (the $\\alpha<0.3$ selection is in "
                "@fig:l2_db@ and @fig:l2_mpf@)", bins.Data(), paren.Data(), cap[a][v].Data()), Form("l2_%s%s", RVTag(v), a == 0 ? "" : "_a1"), 0.32);
    }
  }

  // ======================================================================
  // 5. the dijet JER: DB bisector and MPFX against the truth
  // ======================================================================
  printf("\n=== 5. dijet JER: DB bisector and MPFX, core fit and plain rms, against the truth of pT^corr/pT^linkz ===\n");
  if (jerMerge > 1) printf("--- the JER pT_avg axis merges %d fine bins per point from 3 GeV (drawJVA's third argument; 1 = the fine axis)\n", jerMerge);
  TexSection("Jet energy resolution from dijets", Form("Dijets whose two jets fall into the same $|\\eta|$ bin, with $\\alpha<0.3$, are "
             "projected on the bisector axis $\\hat{n}=(\\hat{u}_1-\\hat{u}_2)/|\\hat{u}_1-\\hat{u}_2|$, along which both jets make the "
             "same angle, and on its normal: $(\\vec{p}_{\\mathrm{T},1}+\\vec{p}_{\\mathrm{T},2})\\cdot\\hat{n}/\\pt^{\\mathrm{avg}}$ for the "
             "DB bisector and $\\vec{p}_{\\mathrm{T}}^{\\,\\mathrm{miss,T1}}\\cdot\\hat{n}/\\pt^{\\mathrm{avg}}$ for MPFX. The parallel "
             "width holds both jets' resolution and the imbalance along the axis, the perpendicular one the imbalance across it, so "
             "the resolution per jet is $\\sigma=\\sqrt{(\\mathrm{RMS}_\\parallel^2-\\mathrm{RMS}_\\perp^2)/2}$; a bin with "
             "$\\mathrm{RMS}_\\perp\\geq\\mathrm{RMS}_\\parallel$ has no JER and counts as a failed bin. Each distribution is "
             "first mirror-symmetrised ($h(x)+h(-x)$: the sign of $\\hat{n}$ follows the order of the two jets, and a $\\pt$-ordered "
             "pair folds the parallel projection at zero); the width is then the core $\\sigma$ of a Gaussian centred at zero, "
             "iterated in $\\pm2\\sigma$ (primary), or the plain in-axis rms. Each needs %.0f dijets per $\\pt^{\\mathrm{avg}}$ bin%s. "
             "The truth is the width of $\\pt^{\\mathrm{corr}}/\\pt^{\\mathrm{linkz}}$ of both jets of the same dijets, like with like: "
             "$\\sigma_{\\mathrm{core}}/\\mu_{\\mathrm{core}}$ under the response protocol for the core-fit JER, rms/mean for the "
             "plain-rms JER. The \\emph{reach} is the lowest $\\pt^{\\mathrm{avg}}$ bin with "
             "$|\\mathrm{JER}/\\mathrm{JER}_{\\mathrm{true}}-1|<%.1f$ whose next two measured bins are not significantly outside that "
             "band ($|\\mathrm{JER}/\\mathrm{JER}_{\\mathrm{true}}-1|<%.1f+2\\sigma$, and no failed bin). Besides the 18 "
             "$|\\eta|$ bins, the same analysis is made in four regions of whole bins (barrel, endcap, transition, HF), whose "
             "same-bin samples are the sums of their bins' (@fig:jer_db_regions@ and @fig:jer_mpfx_regions@).", kMinN,
             jerMerge > 1 ? Form(", on a $\\pt^{\\mathrm{avg}}$ axis of %d merged fine bins from 3\\,GeV (3--5, 5--8, 8--12, "
                                 "12--18\\,GeV, \\ldots)", jerMerge) : "", kJerTol, kJerTol));
  {
    // Per method, class, JER bin: the truth (core; rms/mean), the JER (2
    // estimators x 2 width variants), the Q >= P flags, the ratios.  The core
    // JER is compared with the core truth, the plain-rms JER with the rms
    // truth: on a Gaussian toy (drawrev/toy) the plain-rms bisector JER equals
    // the truth rms/mean to 1% above 12 GeV, while the truth core is 5-10%
    // narrower (a pT_avg bin mixes jets of different true pT), so comparing the
    // rms JER with the core truth reads an estimator mismatch as a method bias.
    static TH1D *jt[kNM][kNDC][kNJB], *jtr[kNM][kNDC][kNJB], *jj[kNM][kNDC][kNJB][kNJ][2], *jq[kNM][kNDC][kNJB][kNJ][2], *jf[kNM][kNDC][kNJB][kNJ][2];
    memset(jt, 0, sizeof(jt)); memset(jtr, 0, sizeof(jtr)); memset(jj, 0, sizeof(jj)); memset(jq, 0, sizeof(jq)); memset(jf, 0, sizeof(jf));
    static double reachRms[kNM][kNDC][kNJB][kNJ];
    for (int m = 0; m < kNM; ++m) for (int c = 0; c < kNDC; ++c) for (int k = 0; k < kNJB; ++k) for (int j = 0; j < kNJ; ++j) reachRms[m][c][k][j] = kReachNA;
    bool haveE[kNJB] = {false}, haveMJ[kNM] = {false};
    int nFail[kNJ][2] = {{0}}, nOvf[kNJ] = {0}, nUsed[kNJ] = {0}; double maxOvf[kNJ] = {0};
    for (int m = 0; m < kNM; ++m) for (int c = 0; c < kNDC; ++c) for (int k = 0; k < kNJB; ++k) {
      const TString base = Form("dijet/%s/%s/", MTag(m), CTag(c));
      TH2 *hl = MergeX(SumEta2(base, "hrlink", JBLo(k), JBHi(k), Form("jl%d%d%d", m, c, k)), jerMerge, Form("jlm%d%d%d", m, c, k));
      jt[m][c][k] = SigCoreVsPt(hl, Form("jt%d%d%d", m, c, k), &jtr[m][c][k]);
      for (int j = 0; j < kNJ; ++j) {
        TH2 *hp = MergeX(SumEta2(base, JHist(j, true), JBLo(k), JBHi(k), Form("jp%d%d%d%d", m, c, k, j)), jerMerge, Form("jpm%d%d%d%d", m, c, k, j));
        TH2 *hq = MergeX(SumEta2(base, JHist(j, false), JBLo(k), JBHi(k), Form("jn%d%d%d%d", m, c, k, j)), jerMerge, Form("jnm%d%d%d%d", m, c, k, j));
        if (!hp || !hq) continue;
        haveE[k] = true; haveMJ[m] = true;
        TH1D *core = Like(hp, Form("jc%d%d%d%d", m, c, k, j)), *rms = Like(hp, Form("jr%d%d%d%d", m, c, k, j));
        TH1D *fc = Like(hp, Form("jfc%d%d%d%d", m, c, k, j)), *fr = Like(hp, Form("jfr%d%d%d%d", m, c, k, j));
        for (int b = 1; b <= hp->GetNbinsX(); ++b) {
          const SymShape P = SymBin(hp, b), Q = SymBin(hq, b);
          if (!P.ok || !Q.ok) continue;
          if (P.core && Q.core && !FillJer(core, fc, b, P.sigC, P.sigCerr, Q.sigC, Q.sigCerr) && c == 0 && k < kNP) ++nFail[j][0];
          if (!FillJer(rms, fr, b, P.rms, P.rms/sqrt(2*std::max(P.neff, 1.)), Q.rms, Q.rms/sqrt(2*std::max(Q.neff, 1.))) && c == 0 && k < kNP) ++nFail[j][1];
          // the plain rms is the in-axis rms: count the distributions with more than 5% outside the axis
          if (c == 0 && k < kNP) { ++nUsed[j]; const double o = std::max(P.ovf, Q.ovf); maxOvf[j] = std::max(maxOvf[j], o); if (o > 0.05) ++nOvf[j]; }
        }
        jj[m][c][k][j][0] = core; jj[m][c][k][j][1] = rms; jf[m][c][k][j][0] = fc; jf[m][c][k][j][1] = fr;
        for (int w = 0; w < 2; ++w) {
          const TH1D *tr = w == 0 ? jt[m][c][k] : jtr[m][c][k];
          if (!tr) continue;
          TH1D *q = Like(hp, Form("jq%d%d%d%d%d", m, c, k, j, w));
          for (int b = 1; b <= q->GetNbinsX(); ++b) if (Has(jj[m][c][k][j][w], b) && Has(tr, b)) {
            const double x = jj[m][c][k][j][w]->GetBinContent(b), y = tr->GetBinContent(b);
            const double r = x/y, er = r*sqrt(pow(jj[m][c][k][j][w]->GetBinError(b)/x, 2) + pow(tr->GetBinError(b)/y, 2));
            q->SetBinContent(b, r); q->SetBinError(b, std::max(er, 1e-9));
          }
          jq[m][c][k][j][w] = q;
        }
        reachRms[m][c][k][j] = Reach(jq[m][c][k][j][1], jf[m][c][k][j][1], jtr[m][c][k]);
      }
      reachDb[m][c][k]   = Reach(jq[m][c][k][kJDB][0], jf[m][c][k][kJDB][0], jt[m][c][k]);
      reachMpfx[m][c][k] = Reach(jq[m][c][k][kJMPFX][0], jf[m][c][k][kJMPFX][0], jt[m][c][k]);
    }
    printf("--- bins with RMS_perp >= RMS_par (no JER; failed in the reach), all vertices, the 18 |eta| bins, every method: "
           "DB core %d, DB rms %d, MPFX core %d, MPFX rms %d\n", nFail[kJDB][0], nFail[kJDB][1], nFail[kJMPFX][0], nFail[kJMPFX][1]);
    for (int j = 0; j < kNJ; ++j) if (nUsed[j])
      printf("--- %s: %d of %d (pT_avg, |eta| bin, method) points have more than 5%% of par or perp outside the -1..1 axis (at most %.0f%%): "
             "their plain rms is truncated, the core fit is clipped to the axis\n", JLabel(j), nOvf[j], nUsed[j], 100*maxOvf[j]);
    // the printout: per JER bin, all vertices, the two truths and the four
    // ratios per method (x: RMS_perp >= RMS_par); the bins without a number
    // are listed on one line
    TString empty;
    for (int k = 0; k < kNJB; ++k) {
      if (!haveE[k]) continue;
      TH1D *any = 0; for (int m = 0; m < kNM && !any; ++m) for (int j = 0; j < kNJ && !any; ++j) if (jj[m][0][k][j][0]) any = jj[m][0][k][j][0];
      if (!any) continue;
      std::vector<int> rowsB;
      for (int b = 1; b <= any->GetNbinsX(); ++b) {
        if (any->GetXaxis()->GetBinLowEdge(b) >= kPtDrawHi) continue;
        bool anyb = false;
        for (int m = 0; m < kNM; ++m) {
          anyb = anyb || Has(jt[m][0][k], b) || Has(jtr[m][0][k], b);
          for (int j = 0; j < kNJ; ++j) for (int w = 0; w < 2; ++w)
            anyb = anyb || Has(jj[m][0][k][j][w], b) || (jf[m][0][k][j][w] && jf[m][0][k][j][w]->GetBinContent(b) > 0);
        }
        if (anyb) rowsB.push_back(b);
      }
      if (rowsB.empty()) { empty += Form(" %s", JBTxt(k).Data()); continue; }
      printf("--- all vertices, %s |eta| %s: truth core / rms, and JER/JER_truth for DB core / DB rms / MPFX core / MPFX rms\n  %-8s",
             k < kNP ? "bin" : "region", JBTxt(k).Data(), "pT_avg");
      for (int m = 0; m < kNM; ++m) if (haveMJ[m]) printf(" | %-6s %-34s", MTag(m), "trC    trR    DBc   DBr   MPc   MPr");
      printf("\n");
      auto cell = [&](int m, int j, int w, int b) -> TString {
        if (jf[m][0][k][j][w] && jf[m][0][k][j][w]->GetBinContent(b) > 0) return "x";
        return Val(jq[m][0][k][j][w], b, "%.2f", "-"); };
      for (int b : rowsB) {
        printf("  %-8s", BinLabel(any->GetXaxis(), b).Data());
        for (int m = 0; m < kNM; ++m) if (haveMJ[m])
          printf(" | %6s %6s %5s %5s %5s %5s      ", Val(jt[m][0][k], b, "%.3f", "-").Data(), Val(jtr[m][0][k], b, "%.3f", "-").Data(),
                 cell(m, kJDB, 0, b).Data(), cell(m, kJDB, 1, b).Data(), cell(m, kJMPFX, 0, b).Data(), cell(m, kJMPFX, 1, b).Data());
        printf("\n");
      }
    }
    if (!empty.IsNull()) printf("--- no pT_avg bin with %.0f same-bin dijets (all vertices) at |eta|%s\n", kMinN, empty.Data());
    // one di-canvas per estimator and JER bin (all vertices, core fit); "" if nothing to draw
    auto drawJer = [&](int j, int k) -> TString {
      int np = 0; double r1 = 1e30, r2 = -1e30, q1 = 1e30, q2 = -1e30; TH1D *any = 0;
      for (int m = 0; m < kNM; ++m) {
        if (HasPoints(jj[m][0][k][j][0]) || HasPoints(jt[m][0][k])) { ++np; if (!any) any = jj[m][0][k][j][0] ? jj[m][0][k][j][0] : jt[m][0][k]; }
        Extent(jj[m][0][k][j][0], r1, r2, kPtDrawHi, 0.15); Extent(jt[m][0][k], r1, r2, kPtDrawHi, 0.15);
        Extent(jq[m][0][k][j][0], q1, q2, kPtDrawHi, 0.25);
      }
      if (!np || !any) return "";
      const double xlo = std::max(AxLo(any), 2.);
      double ylo = 0, yhi = 0.5, zlo = 0.5, zhi = 1.5;
      OpenRange(ylo, yhi, r1, r2, 0.55, 0.); if (q1 < q2) { zlo = std::max(0., std::min(zlo, q1 - 0.05)); zhi = std::min(3., std::max(zhi, q2 + 0.05)); }
      TH1D *up = tdrHist(Form("fju%d%d", j, k), "JER", ylo, yhi, "p_{T}^{avg} [GeV]", xlo, kPtDrawHi);
      TH1D *dw = tdrHist(Form("fjd%d%d", j, k), "JER / JER_{truth}", zlo, zhi, "p_{T}^{avg} [GeV]", xlo, kPtDrawHi);
      up->GetXaxis()->SetLabelSize(0.);
      std::unique_ptr<TCanvas> cv(tdrDiCanvas(Form("c_j%d%d", j, k), up, dw, 8, 11));
      cv->cd(1); gPad->SetLogx();
      int nleg = 0; for (int m = 0; m < kNM; ++m) if (HasPoints(jj[m][0][k][j][0])) ++nleg;
      TLegend *leg = tdrLeg(0.42, 0.56, 0.94, 0.82);
      leg->SetTextSize(0.033); leg->SetNColumns(nleg > 3 ? 2 : 1); leg->SetHeader(Form("%s, %s, all vertices", JLabel(j), JBLabel(k).Data()));
      int iw = 0;
      for (int m = 0; m < kNM; ++m) {
        TGraphErrors *gt = Graph(jt[m][0][k], 0, true, kPtDrawHi);
        if (gt->GetN() > 1) tdrDraw(gt, "LX", kNone, MColour[m], MLine[m] == kSolid ? kDashed : MLine[m], -1, kNone, 0, 0, 1);
        TGraphErrors *g = Graph(jj[m][0][k][j][0], 0.025*(iw - kMidM), true, kPtDrawHi);
        if (g->GetN()) { tdrDraw(g, "Pz", MMarker[m], MColour[m], kSolid, -1, kNone, 0, 0.7); leg->AddEntry(g, MShort(m), "PL"); }
        if (jj[m][0][k][j][0] || jt[m][0][k]) ++iw;
      }
      leg->AddEntry((TObject*)0, "markers core fit, lines truth", "");
      ThresholdsDi();
      cv->cd(2); gPad->SetLogx();
      HLine(xlo, kPtDrawHi); HLine(xlo, kPtDrawHi, 1 - kJerTol, kGray, kDotted); HLine(xlo, kPtDrawHi, 1 + kJerTol, kGray, kDotted);
      iw = 0;
      for (int m = 0; m < kNM; ++m) {
        TGraphErrors *g = Graph(jq[m][0][k][j][0], 0.025*(iw - kMidM), true, kPtDrawHi);
        if (g->GetN()) tdrDraw(g, "Pz", MMarker[m], MColour[m], kSolid, -1, kNone, 0, 0.7);
        // a failed bin (RMS_perp >= RMS_par with the truth measured): a cross on the bottom edge
        if (jf[m][0][k][j][0] && jt[m][0][k]) {
          TGraphErrors *gf = new TGraphErrors();
          const TAxis *a = jf[m][0][k][j][0]->GetXaxis();
          for (int b = 1; b <= a->GetNbins(); ++b)
            if (jf[m][0][k][j][0]->GetBinContent(b) > 0 && Has(jt[m][0][k], b) && a->GetBinLowEdge(b) < kPtDrawHi)
              gf->SetPoint(gf->GetN(), GeoCentre(a, b)*(1 + 0.025*(iw - kMidM)), zlo + 0.06*(zhi - zlo));
          if (gf->GetN()) tdrDraw(gf, "P", kMultiply, MColour[m], kSolid, -1, kNone, 0, 1.4);
        }
        if (jj[m][0][k][j][0] || jt[m][0][k]) ++iw;
      }
      fixOverlay();
      return SavePdf(cv.get(), Form("jva_jer_%s_%s", JTag(j), JBTag(k).Data()));
    };
    auto reachCap = [&](int j, int k) -> TString {
      TString cell;
      for (int m : {kDup, kLv, kJva, kOPart, kNobody}) {
        const double r = j == kJDB ? reachDb[m][0][k] : reachMpfx[m][0][k];
        if (r != kReachNA) cell += Form("%s%s %s", cell.IsNull() ? "" : ", ", MTex(m), r > 0 ? Form("%g\\,GeV", r) : "none");
      }
      return cell.IsNull() ? TString("") : TString(Form(" Reach at %s: %s.", JBTex(k).Data(), cell.Data()));
    };
    for (int j = 0; j < kNJ; ++j) {
      // the four regions first (the statistics), then the 18 bins in two figures of nine
      std::vector<TString> reg;
      TString capr;
      for (int k = kNP; k < kNJB; ++k) { if (haveE[k]) drawJer(j, k); reg.push_back(Form("jva_jer_%s_%s", JTag(j), JBTag(k).Data())); capr += reachCap(j, k); }
      TexFigure(reg, Form("JER per jet from the %s method against $\\pt^{\\mathrm{avg}}$ for the same-bin dijets of all vertices "
                "($\\alpha<0.3$) in four regions made of whole $|\\eta|$ bins, %s, %s, %s and %s (left to right, top to bottom; "
                "the regions without %.0f dijets in any $\\pt^{\\mathrm{avg}}$ bin are left out): "
                "core-fit JER per method (markers) and the truth $\\sigma_{\\mathrm{core}}/\\mu_{\\mathrm{core}}$ of "
                "$\\pt^{\\mathrm{corr}}/\\pt^{\\mathrm{linkz}}$ of the same jets (lines, same colour); below, the ratio with "
                "$\\pm10\\%%$ dotted, and a cross on the bottom edge where $\\mathrm{RMS}_\\perp\\geq\\mathrm{RMS}_\\parallel$.%s",
                JLabel(j), JBTex(kNP).Data(), JBTex(kNP+1).Data(), JBTex(kNP+2).Data(), JBTex(kNP+3).Data(), kMinN,
                capr.Data()), Form("jer_%s_regions", JTag(j)), 0.48);
      for (int ie = 0; ie < kNP; ++ie) if (haveE[ie]) drawJer(j, ie);
      for (int half = 0; half < 2; ++half) {
        std::vector<TString> part; TString capj;
        for (int ie = 9*half; ie < 9*(half + 1); ++ie) {
          part.push_back(Form("jva_jer_%s_%s", JTag(j), PTag(ie).Data()));
          if (ie == 9*half || ie == 11 || ie == 15) capj += reachCap(j, ie);
        }
        TexFigure(part, Form("JER per jet from the %s method against $\\pt^{\\mathrm{avg}}$ for the same-bin dijets of all vertices "
                  "($\\alpha<0.3$), per $|\\eta|$ bin from %s to %s (left to right, top to bottom; the bins without %.0f dijets "
                  "in any $\\pt^{\\mathrm{avg}}$ bin are left out); markers, lines and the lower panels as in @fig:jer_%s_regions@.%s",
                  JLabel(j), PTex(9*half).Data(), PTex(9*half + 8).Data(), kMinN, JTag(j), capj.Data()), Form("jer_%s_%d", JTag(j), half), 0.32);
      }
    }
    // the distributions behind one point: jva (else the first method) in the barrel and the transition regions, at 8 GeV
    {
      TString capd, ptd;
      int md = -1; for (int m : {kJva, kOPart, kLv, kDup, kTrk, kOJet, kNobody}) if (haveMJ[m]) { md = m; break; }
      const int kd[2] = {kNP, kNP + 2};
      for (int j = 0; md >= 0 && j < kNJ; ++j) for (int k : kd) {
        const TString base = Form("dijet/%s/all/", MTag(md));
        TH2 *hp = MergeX(SumEta2(base, JHist(j, true), JBLo(k), JBHi(k), Form("dp%d%d", j, k)), jerMerge, Form("dpm%d%d", j, k));
        TH2 *hq = MergeX(SumEta2(base, JHist(j, false), JBLo(k), JBHi(k), Form("dq%d%d", j, k)), jerMerge, Form("dqm%d%d", j, k));
        if (!hp || !hq) continue;
        const int b = BinAt(hp, 8.);
        TH1D *sp = 0, *sq = 0;
        const SymShape P = SymBin(hp, b, &sp), Q = SymBin(hq, b, &sq);
        if (!P.ok || !Q.ok || !sp || !sq) continue;
        ptd = BinTex(hp->GetXaxis(), b);
        const double np = sp->Integral(0, sp->GetNbinsX()+1), nq = sq->Integral(0, sq->GetNbinsX()+1);
        TH1D *up = Unit(sp, Form("udp%d%d", j, k)), *uq = Unit(sq, Form("udq%d%d", j, k));
        const double ymax = std::max(up->GetMaximum(), uq->GetMaximum());
        const double x1 = sp->GetXaxis()->GetXmin(), x2 = sp->GetXaxis()->GetXmax();
        // the legend needs the top 40% of the pad: the curves in the lower 3 of 5.3 decades
        TH1D *frame = tdrHist(Form("fdd%d%d", j, k), "Fraction of dijets", ymax*1e-3, ymax*200,
                              j == kJDB ? "(p_{T,1} + p_{T,2}) #upoint #hat{n} / p_{T}^{avg}" : "p_{T}^{miss,T1} #upoint #hat{n} / p_{T}^{avg}", x1, x2);
        frame->GetXaxis()->SetNdivisions(505);
        std::unique_ptr<TCanvas> cv(tdrCanvas(Form("c_dd%d%d", j, k), frame, 8, 11, kSquare));
        cv->SetLogy();
        tdrDraw(up, "HIST", kNone, kRed+1, kSolid, -1, kNone, 0, 0, 2);
        tdrDraw(uq, "HIST", kNone, kBlue+1, kSolid, -1, kNone, 0, 0, 2);
        TLegend *leg = tdrLeg(0.20, 0.58, 0.94, 0.79);
        leg->SetTextSize(0.028);
        leg->SetHeader(Form("%s, %s, %s GeV, mirror-symmetrised", MShort(md), JBLabel(k).Data(), BinLabel(hp->GetXaxis(), b).Data()));
        leg->AddEntry(up, Form("parallel: rms %.3f, #sigma_{c} %s", P.rms, P.core ? Form("%.3f", P.sigC) : "fail"), "L");
        leg->AddEntry(uq, Form("perpendicular: rms %.3f, #sigma_{c} %s", Q.rms, Q.core ? Form("%.3f", Q.sigC) : "fail"), "L");
        for (int t = 0; t < 2; ++t) {
          const SymShape &S = t == 0 ? P : Q;
          if (!S.core) continue;
          TF1 *f = new TF1(Form("fdg%d%d%d", j, k, t), "gaus", std::max(x1, -2*S.sigC), std::min(x2, 2*S.sigC));
          f->SetParameters(S.amp/(t == 0 ? np : nq), 0., S.sigC);
          f->SetLineColor(t == 0 ? kRed+1 : kBlue+1); f->SetLineStyle(kDashed); f->SetLineWidth(2);
          f->Draw("SAME");
        }
        TString jer = "not defined";
        if (P.core && Q.core && P.sigC > Q.sigC) jer = Form("%.3f", sqrt((P.sigC*P.sigC - Q.sigC*Q.sigC)/2));
        leg->AddEntry((TObject*)0, Form("JER (core) = #sqrt{(#sigma_{#parallel}^{2} - #sigma_{#perp}^{2}) / 2} = %s", jer.Data()), "");
        Thresholds();
        fixOverlay();
        SavePdf(cv.get(), Form("jva_jerdist_%s_%s", JTag(j), JBTag(k).Data()));
        printf("--- %s %s |eta| %s, %s GeV: par rms %.4f core %s, perp rms %.4f core %s, JER core %s\n", JTag(j), MTag(md), JBTxt(k).Data(),
               BinLabel(hp->GetXaxis(), b).Data(), P.rms, P.core ? Form("%.4f", P.sigC) : "fail", Q.rms, Q.core ? Form("%.4f", Q.sigC) : "fail", jer.Data());
        capd += Form(" %s at %s: JER %s.", JLabel(j), JBTex(k).Data(), jer.Data());
      }
      std::vector<TString> order;
      for (int j = 0; j < kNJ; ++j) for (int k : kd) order.push_back(Form("jva_jerdist_%s_%s", JTag(j), JBTag(k).Data()));
      if (md >= 0) TexFigure(order, Form("The mirror-symmetrised parallel (red) and perpendicular (blue) projections behind one point, "
                             "for %s at $\\pt^{\\mathrm{avg}}$ %s\\,GeV, DB bisector (first row) and MPFX (second row), in the "
                             "regions %s and %s, with the core Gaussians (dashed) that give the widths.%s",
                             MTex(md), ptd.Data(), JBTex(kd[0]).Data(), JBTex(kd[1]).Data(), capd.Data()), "jerdist", 0.48, {2, 2});
    }
    // the reach tables: rows the |eta| bins, then the regions (only the rows
    // with something measured); columns method x estimator (core; then rms)
    int nmj = 0; TString mh; std::vector<TString> mtxt;
    for (int m = 0; m < kNM; ++m) if (haveMJ[m]) { ++nmj; mh += Form(" & \\multicolumn{2}{c}{%s}", MTex(m));
                                                    mtxt.push_back(Form("%s:DB", MTag(m))); mtxt.push_back(Form("%s:MPFX", MTag(m))); }
    mtxt.push_back("PVlv:DB"); mtxt.push_back("PVlv:MPFX");
    if (nmj) for (int w = 0; w < 2; ++w) {
      TString sub = " "; for (int m = 0; m < kNM; ++m) if (haveMJ[m]) sub += " & DB & MPFX";
      TexTableBegin(Form("Reach of the dijet JER with the %s widths (against the truth %s): the lowest $\\pt^{\\mathrm{avg}}$ (lower "
                         "bin edge, GeV) with $|\\mathrm{JER}/\\mathrm{JER}_{\\mathrm{true}}-1|<%.1f$ whose next two measured bins "
                         "are not significantly outside that band, per $|\\eta|$ bin and, below the line, per region of whole bins, "
                         "per method and estimator (DB bisector, MPFX), for the same-bin dijets of all vertices and, in the last two "
                         "columns, of the primary vertex under \\texttt{lv} (the standard reconstruction). \\emph{none}: measured, "
                         "never within the band; --: no $\\pt^{\\mathrm{avg}}$ bin with both the JER and the truth. Only the rows "
                         "with a measurement are listed.",
                         w == 0 ? "core-fit" : "plain-rms", w == 0 ? "$\\sigma_{\\mathrm{core}}/\\mu_{\\mathrm{core}}$" : "rms/mean",
                         kJerTol), w == 0 ? "reach" : "reach_rms",
                    Form("l%srr", TString('r', 2*nmj).Data()),
                    TString("$|\\eta|$") + mh + " & \\multicolumn{2}{c}{PV, \\texttt{lv}} \\\\\n" + sub + " & DB & MPFX",
                    mtxt, 10, "\\scriptsize\\setlength{\\tabcolsep}{3pt}");
      bool lined = false;
      for (int k = 0; k < kNJB; ++k) {
        if (!haveE[k]) continue;
        std::vector<TString> cells; bool anyr = false;
        auto rv = [&](int m, int c, int j) { return w == 0 ? (j == kJDB ? reachDb[m][c][k] : reachMpfx[m][c][k]) : reachRms[m][c][k][j]; };
        for (int m = 0; m < kNM; ++m) if (haveMJ[m]) for (int j = 0; j < kNJ; ++j) { anyr = anyr || rv(m, 0, j) != kReachNA; cells.push_back(ReachStr(rv(m, 0, j))); }
        for (int j = 0; j < kNJ; ++j) { anyr = anyr || rv(kLv, 1, j) != kReachNA; cells.push_back(ReachStr(rv(kLv, 1, j))); }
        if (!anyr) continue;
        if (k >= kNP && !lined && gTabRows > 0) { TexHline(); lined = true; }
        TexRow(JBTex(k), JBTxt(k), cells);
      }
      TexTableEnd();
    }
  }

  // ======================================================================
  // 6. summary
  // ======================================================================
  printf("\n=== 6. summary ===\n");
  FlushSection();
  {
    TexTableBegin("Summary per method, vertex level: rms and core $\\sigma$ of the missing-$\\pt$ components over all vertices, the core "
                  "$\\sigma$ of their residual against the owner interaction's true imbalance (GeV); the reconstructed over generated "
                  "jets per owned vertex at $2.5<|\\eta|<3.0$ above 10\\,GeV, and the mean own-interaction fraction of the linked "
                  "$\\pt$ of the jets there above 5\\,GeV; the correct "
                  "fraction of the forward clusters above 5\\,GeV (right vertex plus combination nulled), against which the last "
                  "row is the reference of dropping every cluster.", "summary", "lrrrrrr",
                  "method & rms & $\\sigma_{\\mathrm{core}}$ & $\\sigma_{\\mathrm{core,res}}$ & reco/gen$_{>10}$ & own & correct",
                  {"rms", "core", "core_res", "reco/gen10", "own", "correct"}, 10);
    for (int m = 0; m < kNM; ++m) {
      if (!haveMet[m] && !haveMetRes[m] && !haveSpike[m] && !haveOwn[m] && !haveCorrect[m]) continue;
      TexRow(MTex(m), MTag(m), {Num(haveMet[m], metRms[m], "%.2f"), haveMet[m] ? Num(haveMetCore[m], metCore[m], "%.2f", "fail") : TString("--"),
                                haveMetRes[m] ? Num(haveResCore[m], resCore[m], "%.2f", "fail") : TString("--"),
                                Num(haveSpike[m], spike10[m], "%.2f"), Num(haveOwn[m], own1[m], "%.3f"), Num(haveCorrect[m], correct[m], "%.3f")});
    }
    if (dropAllRef >= 0 && !haveCorrect[kNobody]) TexRow("drop all (ref.)", "dropall", {"--", "--", "--", "--", "--", Form("%.3f", dropAllRef)});
    TexTableEnd();
    TString eh; std::vector<TString> etxt = {"dij/X all", "dij/X PV", "clos a<0.3", "clos all a"};
    for (int r = 0; r < kNJR; ++r) { eh += Form(" & %g--%g", kProbeEta[JRLo[r]], kProbeEta[JRHi[r]+1]); etxt.push_back(Form("reach %s", JBTxt(kNP + r).Data())); }
    TexTableBegin(Form("Summary per method, dijet level: same-bin dijets per crossing (all vertices / primary vertex); the L2Res "
                       "closure, the mean $|R_{\\mathrm{MPF}}/r_{\\mathrm{true}}-1|$ over the probe bins at $\\pt^{\\mathrm{avg}}$ "
                       "8--10\\,GeV (all vertices, with $\\alpha<0.3$ and without the $\\alpha$ cut; number of probe bins in "
                       "parentheses); and the reach of the core-fit MPFX JER (GeV) in the four regions of whole $|\\eta|$ bins, all "
                       "vertices / primary vertex under the method (\\emph{none}: measured, never within %.0f\\%%; --: not measured).", 100*kJerTol),
                  "summary_dijet", Form("lrrrr%s", TString('r', kNJR).Data()),
                  TString("method & \\multicolumn{2}{c}{dijets per crossing} & \\multicolumn{2}{c}{closure} & \\multicolumn{4}{c}{MPFX JER reach at $|\\eta|$, all / PV} \\\\\n"
                          " & all & PV & $\\alpha<0.3$ & all $\\alpha$") + eh,
                  etxt, 12, "\\footnotesize");
    for (int m = 0; m < kNM; ++m) {
      bool any = haveDij[m][0] || haveDij[m][1] || haveL2[m][0] || haveL2a1[m][0];
      for (int r = 0; r < kNJR; ++r) any = any || reachMpfx[m][0][kNP + r] != kReachNA || reachMpfx[m][1][kNP + r] != kReachNA;
      if (!any) continue;
      std::vector<TString> cells = {Num(haveDij[m][0], dijPerX[m][0], "%.4f"), Num(haveDij[m][1], dijPerX[m][1], "%.4f"),
                                    haveL2[m][0] ? TString(Form("%.3f (%d)", l2mpf[m][0], l2n[m][0])) : TString("--"),
                                    haveL2a1[m][0] ? TString(Form("%.3f (%d)", l2mpfa1[m][0], l2na1[m][0])) : TString("--")};
      for (int r = 0; r < kNJR; ++r)
        cells.push_back(ReachStr(reachMpfx[m][0][kNP + r]) + " / " + ReachStr(reachMpfx[m][1][kNP + r]));
      TexRow(MTex(m), MTag(m), cells);
    }
    TexTableEnd();
  }

  if (gTabDone.count("summary") || gTabDone.count("summary_dijet"))
    TexSection("Summary", "The key number of each section per method is collected in @tab:summary@ (vertex level) and "
               "@tab:summary_dijet@ (dijets).", true);
  FlushSection();
  gTexP.close(); gTexT.close();
  printf("\n=== core-fit protocol, responses: %d distributions with n >= %.0f, %d converged, %d failed as sigma_core > 1.1 rms, "
         "%d as mu_core outside [0.5, 1.3] median\n", gNDesc, kMinN, gNConv, gNWide, gNMu);
  printf("=== core-fit protocol, zero-centred: %d distributions with n >= %.0f, %d converged, %d failed as wider than the distribution, "
         "%d as |mu_core - median| > 0.5 sigma_core\n", gNDescS, kMinN, gNConvS, gNWideS, gNMuS);
  if (gMissing > kMaxReport) {
    printf("=== missing or mistyped histograms per directory:");
    for (const auto &d : gMissingDir) printf(" %s: %d;", d.first.c_str(), d.second);
    printf("\n");
  }
  printf("\ndrawJVA: wrote %d PDFs to plots/ (%d cropped), doc/jva_plots.tex and doc/jva_tables.tex from %s (tag %s); %d histogram(s) missing\n",
         gNPdf, gNCrop, fname, tag, gMissing);
}
