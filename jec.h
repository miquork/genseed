// genseed - Copyright 2026 Mikko Voutilainen - Apache License 2.0
#ifndef genseed_jec_h
#define genseed_jec_h
// jec.h - the jet energy correction this analysis derives for itself.
// genseed port of lowptjets/jec.h (v4): same table format and interpolation,
// TString/Form replaced by std::string/snprintf so that the header, like
// puppi.h and genlink.h, compiles without ROOT.
//
// WHY NOT A STANDARD ONE.  The published AK4 PUPPI corrections start at
// 15 GeV, because that is where CMS stores and calibrates jets; below it they
// are an extrapolation of a fit, over exactly the range where the pileup
// offset, the tracking efficiency and the neutral thresholds all change
// fastest.  This analysis reconstructs jets from 1 GeV with its own PUPPI
// weights, so it has neither a valid correction nor the right one: the jets
// are not CMS's jets.  The correction is therefore measured here, in the
// MC-truth way, against the generated energy the linker routes into each
// reco jet: pT^link, at generated scale (genlink.h, GenSeed.h).
//
// THE LOW-pT BIAS, WHICH IS THE WHOLE DIFFICULTY.  The response is
// R = pT^raw / pT^link in bins of pT^link, and the correction is 1/R.  At high
// pT the distribution of R is near Gaussian and its mean, median and Gaussian
// core all agree.  At low pT they do not, for a reason that has nothing to do
// with the shape: a jet whose response fluctuates low is not reconstructed at
// all - it falls under the clustering threshold (and a dR pairing loses the
// match as well) - so the left tail of R is cut away and every estimator that
// averages over what survives comes out too high.  Correcting with a response
// that is too high makes the corrected pT too low, and the error grows as the
// threshold is approached, which is precisely the region this analysis is
// about.
//
// Three things are done about it, in order of how much they help:
//
//  1. Cluster far below the measurement.  The reco jets are clustered from
//     1 GeV although the spectrum starts at 5, so a 5 GeV jet has to
//     fluctuate down by a factor five before it disappears.  This is the only
//     one of the three that removes the truncation rather than working around
//     it, and it is available here only because the jets are reclustered from
//     the candidates.  It is not available to an analysis that reads jets from
//     NanoAOD, which is why the effect is normally fought rather than avoided.
//  2. Use the median of R, not the mean.  The median is unmoved by a tail
//     being cut as long as less than half the distribution goes, and it does
//     not need a fit to converge.
//  3. Take only what can be measured.  A pT^link bin becomes a node only with
//     at least 200 jets in it (and a median above 0.05), and the correction
//     is never extrapolated past the nodes: it is flat outside them.
//
// The median is the only estimator.  v4 also kept the arithmetic mean, the
// Gaussian core fit and a mean after an artificial 5 GeV reco threshold as
// alternative corrections, drawn together by drawJEC.C; genseed has none of
// those variants and no drawJEC.C.  A Gaussian fit is the textbook answer to
// a truncated distribution, but it was never the default: at 5-10 GeV with a
// few hundred jets per bin the fit range has to be chosen from the histogram
// it is fitting, and the result moves by more than the effect being
// corrected.  drawGenSeed.C still fits the core of R, but only to describe
// its width, not to correct with.
//
// FORM OF THE CORRECTION.  What is measured is R(pT^link); what is needed is a
// factor applied to pT^raw.  drawGenSeed.C builds the table by walking the
// pT^link bins of resp/hresp_link_<e> that pass the n >= 200 cut: each one
// contributes the node
//     ( pT^raw = median(R) * pT^link ,  correction = 1/median(R) )
// with pT^link the geometric centre of the bin, i.e. the correction is
// 1/median(pT^raw/pT^link) per pT^link bin, re-indexed to pT^raw.  Here it is
// then linear in log pT^raw between nodes and flat outside the measured range;
// an |eta| bin with fewer than two nodes is left at 1.  One table per |eta|
// bin, written as plain text to text/jec<tag>.txt (lines of "ietabin pT_raw
// correction", # comments) by drawGenSeed.C and read back here.  Because the
// nodes are against pT^link, pass 2 corrects the reco jet to the seed scale.
#include <cstdio>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>

namespace jec {

  // |eta| binning of the correction.  Coarser than the response is measured
  // in would hide the HF transition; finer runs out of jets.
  const int    kNEta = 10;
  const double kEtaW = 0.5;
  inline int   EtaBin(double eta) {
    const int i = int(fabs(eta)/kEtaW);
    return (i >= 0 && i < kNEta) ? i : -1;
  }
  inline std::string EtaTag(int i) {
    char b[16]; snprintf(b, sizeof(b), "e%02d", i); return b; }
  inline std::string EtaLabel(int i) {
    char b[64]; snprintf(b, sizeof(b), "%.1f < |#eta| < %.1f", i*kEtaW, (i+1)*kEtaW);
    return b; }

  class Correction {
  public:
    Correction() : fLoaded(false) {}

    bool Load(const char *fname)
    {
      std::ifstream in(fname);
      if (!in) { printf("jec: no correction at %s, using 1.0\n", fname);
                 return false; }
      for (int i = 0; i < kNEta; ++i) { fPt[i].clear(); fC[i].clear(); }
      std::string line;
      int n = 0;
      while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        int ie; double pt, c;
        if (!(ss >> ie >> pt >> c)) continue;
        if (ie < 0 || ie >= kNEta) continue;
        fPt[ie].push_back(pt); fC[ie].push_back(c); ++n;
      }
      // the nodes have to be in pT order for the interpolation
      for (int i = 0; i < kNEta; ++i) {
        std::vector<std::pair<double,double> > v;
        for (size_t k = 0; k < fPt[i].size(); ++k)
          v.push_back(std::make_pair(fPt[i][k], fC[i][k]));
        std::sort(v.begin(), v.end());
        fPt[i].clear(); fC[i].clear();
        for (size_t k = 0; k < v.size(); ++k) {
          // a node at the same pT as the last one would divide by zero
          if (!fPt[i].empty() && v[k].first <= fPt[i].back()*1.0001) continue;
          fPt[i].push_back(v[k].first); fC[i].push_back(v[k].second);
        }
      }
      fLoaded = n > 0;
      printf("jec: %d nodes from %s\n", n, fname);
      return fLoaded;
    }

    bool IsLoaded() const { return fLoaded; }

    // Correction factor for a raw jet.  Flat outside the measured range: an
    // extrapolated correction below the lowest measured point is exactly the
    // thing this file exists to avoid.
    double operator()(double ptraw, double eta) const
    {
      if (!fLoaded) return 1.;
      const int i = EtaBin(eta);
      if (i < 0 || fPt[i].size() < 2) return 1.;
      const std::vector<double> &x = fPt[i], &y = fC[i];
      if (ptraw <= x.front()) return y.front();
      if (ptraw >= x.back())  return y.back();
      const size_t k = std::upper_bound(x.begin(), x.end(), ptraw) - x.begin();
      const double t = log(ptraw/x[k-1]) / log(x[k]/x[k-1]);
      return y[k-1] + t*(y[k]-y[k-1]);
    }

  private:
    bool   fLoaded;
    std::vector<double> fPt[kNEta], fC[kNEta];
  };

} // namespace jec
#endif
