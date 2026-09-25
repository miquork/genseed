// genseed - Copyright 2026 Mikko Voutilainen - Apache License 2.0
// cmpPuppiFwd.C - the forward spike and the vertex MET under PUPPI's default
// vertex-blind forward treatment and with its vertex association extended into
// the forward pixel coverage, from JVA productions made with different PUPPI
// knobs (JVAOPT, runJVA.C):
//   root -l -b -q 'cmpPuppiFwd.C("rootfiles/JVA_v2.root,rootfiles/JVA_vtx30.root,rootfiles/JVA_vtx30f.root",
//                               "default,|eta|<3 assoc,|eta|<3 assoc + central floor")'
// Per method: MET rms per component over all vertices and at the PV, its
// residual against the owner's true imbalance, reco/gen jets per owned vertex
// at 2.5<|eta|<3.0 above 5/10/20 GeV, the own-interaction fraction of those
// jets, and the reco/gen ratio in fine |eta| bins across the transition.
#include <TFile.h>
#include <TH1.h>
#include <TProfile.h>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <cstdio>
#include <cmath>
#include <vector>

namespace {
  double Rms(TH1 *h) { return h && h->GetEntries() > 0 ? h->GetRMS() : -1; }
  double Sum(TH1 *h, double lo, double hi) {
    if (!h) return -1;
    return h->Integral(h->FindFixBin(lo + 1e-6), h->FindFixBin(hi - 1e-6));
  }
  double ProfMean(TProfile *p, double lo, double hi) {
    if (!p) return -1;
    double s = 0, n = 0;
    for (int b = p->FindFixBin(lo + 1e-6); b <= p->FindFixBin(hi - 1e-6); ++b) {
      s += p->GetBinContent(b)*p->GetBinEntries(b); n += p->GetBinEntries(b); }
    return n > 0 ? s/n : -1;
  }
}

void cmpPuppiFwd(const char *files = "rootfiles/JVA_v2.root,rootfiles/JVA_vtx30.root,rootfiles/JVA_vtx30f.root",
                 const char *labels = "default,|eta|<3 assoc,|eta|<3 assoc + central floor")
{
  TObjArray *fa = TString(files).Tokenize(","), *la = TString(labels).Tokenize(",");
  const char *meth[] = {"dup", "lv", "jva", "ojet", "opart", "none"};
  for (int f = 0; f < fa->GetEntries(); ++f) {
    const TString fn = ((TObjString*)fa->At(f))->GetString();
    const TString lab = f < la->GetEntries() ? ((TObjString*)la->At(f))->GetString() : fn;
    TFile *tf = TFile::Open(fn);
    if (!tf || tf->IsZombie()) { printf("cannot open %s\n", fn.Data()); continue; }
    printf("\n=== %s (%s)\n", lab.Data(), fn.Data());
    printf("%-6s %8s %8s %8s | %8s %8s %8s | %6s %6s %6s\n", "method", "MET all", "MET PV", "res all",
           "r/g >5", "r/g >10", "r/g >20", "own", "own lv", "");
    double g[3];
    const int pts[3] = {5, 10, 20};
    for (int k = 0; k < 3; ++k) g[k] = Sum((TH1*)tf->Get(Form("jets/gen/hjeteta_pt%d", pts[k])), 2.5, 3.0);
    for (const char *m : meth) {
      TH1 *ha = (TH1*)tf->Get(Form("met/%s/hmetxy_all", m));
      if (!ha) continue;
      double r[3];
      for (int k = 0; k < 3; ++k) { const double s = Sum((TH1*)tf->Get(Form("jets/%s/hjeteta_pt%d", m, pts[k])), 2.5, 3.0); r[k] = g[k] > 0 ? s/g[k] : -1; }
      printf("%-6s %8.2f %8.2f %8.2f | %8.2f %8.2f %8.2f | %6.3f\n", m, Rms(ha), Rms((TH1*)tf->Get(Form("met/%s/hmetxy_pv", m))),
             Rms((TH1*)tf->Get(Form("met/%s/hmetres_all", m))), r[0], r[1], r[2],
             ProfMean((TProfile*)tf->Get(Form("jets/%s/pown_eta", m)), 2.5, 3.0));
    }
    // the shape across the transition, pT > 5 GeV, lv (the standard: PV keeps the vertex-blind part)
    TH1 *hg = (TH1*)tf->Get("jets/gen/hjeteta_pt5");
    printf("reco/gen, pT > 5, per |eta| bin 1.8-3.6:  ");
    for (const char *m : {"dup", "lv", "none"}) {
      TH1 *hr = (TH1*)tf->Get(Form("jets/%s/hjeteta_pt5", m));
      if (!hr || !hg) continue;
      printf("\n   %-5s", m);
      for (int b = hg->FindFixBin(1.81); b <= hg->FindFixBin(3.59); ++b)
        printf(" %5.2f", hg->GetBinContent(b) > 0 ? hr->GetBinContent(b)/hg->GetBinContent(b) : -1.);
    }
    printf("\n   |eta| lo:");
    for (int b = hg->FindFixBin(1.81); b <= hg->FindFixBin(3.59); ++b) printf(" %5.1f", hg->GetBinLowEdge(b));
    printf("\n");
  }
}
