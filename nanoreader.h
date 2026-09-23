// genseed - Copyright 2026 Mikko Voutilainen - Apache License 2.0
#ifndef genseed_nanoreader_h
#define genseed_nanoreader_h
// nanoreader.h - the one place in genseed that knows a NanoAOD branch name.
//
// WHY A READER OF ITS OWN.  The v4 analysis read the tree through a MakeClass
// holder (lowptjets.h) with hand-sized arrays: 120 interaction slots x 2048
// particles x five floats, 16384 PF candidates x twelve branches, 3 MB of
// static storage sized from what one production happened to write, with a
// CheckSizes() to refuse an entry that would overflow.  That is the right
// defence when the arrays are fixed, and the wrong design when the same code
// has to run on a production with 200 slots or on JMENANO with other names.
// TTreeReader sizes every array per entry from its counter branch, reads only
// the branches that have a reader attached (the 600 per-interaction branches
// are the whole cost of this tree), and gives a range [first, last) on a
// TChain for free, which is how the jobs are split (with one trap, see
// SetEntriesRange: an empty range must be caught here).  The number of interaction
// slots is not a constant here: the constructor walks nGenPartCandBX0PUEvent<i>
// from i = 0 until the first missing one, so a tree with 51 or 120 or 200
// slots is read whole, and one with none says so.
//
// WHAT THE SAMPLE CONTAINS (MC24NanoV15_PU_IT_OOT / NANOAODSIM_1.root, 25000
// entries).  A NeutrinoGun with the Run 3 pileup library mixed in and the
// generated particles of EVERY pileup interaction kept in a branch of its own:
//   GenPartCandBX0PUEvent<i>_{pt,eta,phi,mass,pdgId}   i = 0 ... 119
//   GenVtxBX0PUEvent<i>_{x,y,z,t0}
// Slot 0 is the hard scatter, i.e. the neutrino itself (one particle); slots
// 1 ... Pileup_nPU are the in-time interactions, one each; slots above
// Pileup_nPU are empty.  <Pileup_nPU> = 45.4, the largest interaction has
// 1121 particles and the largest event 18503.  Despite the directory name
// every per-interaction branch is BX0: there is no out-of-time content to
// read.  GenJet is empty (made from the hard scatter, a neutrino), so the gen
// jets are clustered by the analyzer, which is the point.  PFCandV2 is the
// whole event, 2103 candidates on average and 3397 at the top of 500 events;
// PFCandV2_pt is the raw PF pT, PFCandV2_puppiWeight the weight CMS computed
// against its primary vertex and puppi.h recomputes for every other one.
//
// WHAT FillPF AND FillGen ARE FOR.  The analyzer never touches a branch: it
// asks for a puppi::Event (FillPF: pt eta phi mass, ztrk = PV_z + dz, charged,
// pdg, q, vertexRef, fromPVvertexRef, pvAssocQuality, nPV) and for a
// genlink::GenList (FillGen: all interactions flattened, |eta| < etaMax,
// neutrinos optionally dropped, z of the interaction, the charge sign from
// TDatabasePDG).  Both are what v4's LoadCandidates() and ReadInteractions()
// did, cut by cut, so the v4 numbers carry over: entry 9374 gives nPU = 50,
// 2573 PF candidates, PV_z = 4.488 cm, 51 filled slots of the 120 in the tree
// and 9636 generated particles, of which 7738 are inside |eta| < 5.5 and not
// neutrinos (7 neutrinos, 1891 beyond 5.5).  FillGen keeps the particles
// of one interaction contiguous and in slot order, so the analyzer clusters
// each interaction from the run of equal GenList::inter.  A nucleus (pdgId
// above 1e9) is not in TDatabasePDG and gets charge 0, exactly as in v4.
//
// A MISSING BRANCH IS LOUD.  A TTreeReaderValue on a branch that is not there
// makes every Next() fail with a message that names the wrong branch, but
// only at the first event; the constructor checks each required branch up
// front, prints the missing ones and sets Ok() false, so a wrong production is
// caught before any library is loaded and any job started.
#include "puppi.h"
#include "genlink.h"

#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>

// NAMESPACE.  The spec said nano::Reader, but rootcling compiles every ACLiC
// dictionary with "using namespace std" in effect and std::nano (the <ratio>
// typedef) makes nano:: ambiguous in every macro that names it; nanoaod:: is
// unambiguous everywhere.
namespace nanoaod {

  class Reader {
  public:
    // The tree (or TChain) must outlive the reader.  Prefix names the
    // per-interaction branches: <prefix><i>_pt etc. with counter n<prefix><i>.
    explicit Reader(TTree *tree, const char *prefix = "GenPartCandBX0PUEvent",
                    const char *vtxPrefix = "GenVtxBX0PUEvent")
      : fTree(tree), fReader(tree), fOk(true), fNSlots(0)
    {
      // event
      fRun   = Val<UInt_t>("run");
      fLumi  = Val<UInt_t>("luminosityBlock");
      fEvent = Val<ULong64_t>("event");
      fNPU   = Val<Int_t>("Pileup_nPU");
      fPVz   = Val<Float_t>("PV_z");
      fNPV   = Val<UChar_t>("PV_npvs");
      // particle flow
      fPFpt   = Arr<Float_t>("PFCandV2_pt");
      fPFeta  = Arr<Float_t>("PFCandV2_eta");
      fPFphi  = Arr<Float_t>("PFCandV2_phi");
      fPFmass = Arr<Float_t>("PFCandV2_mass");
      fPFdz   = Arr<Float_t>("PFCandV2_dz");
      fPFw    = Arr<Float_t>("PFCandV2_puppiWeight");
      fPFpdg  = Arr<Int_t>("PFCandV2_pdgId");
      fPFq    = Arr<Int_t>("PFCandV2_charge");
      fPFvref = Arr<Int_t>("PFCandV2_vertexRef");
      fPFfpv  = Arr<Int_t>("PFCandV2_fromPVvertexRef");
      fPFfpv0 = Arr<Int_t>("PFCandV2_fromPV0");
      fPFpaq  = Arr<Int_t>("PFCandV2_pvAssocQuality");
      // generated interactions: walk the slots until the first missing counter
      for (int i = 0; ; ++i) {
        const std::string n = "n" + std::string(prefix) + std::to_string(i);
        if (!tree->GetBranch(n.c_str())) break;
        const std::string p = std::string(prefix) + std::to_string(i) + "_";
        const std::string v = std::string(vtxPrefix) + std::to_string(i) + "_z";
        fGpt  .push_back(Arr<Float_t>((p+"pt").c_str()));
        fGeta .push_back(Arr<Float_t>((p+"eta").c_str()));
        fGphi .push_back(Arr<Float_t>((p+"phi").c_str()));
        fGmass.push_back(Arr<Float_t>((p+"mass").c_str()));
        fGpdg .push_back(Arr<Int_t>((p+"pdgId").c_str()));
        fGvz  .push_back(Val<Float_t>(v.c_str()));
        ++fNSlots;
      }
      if (fNSlots == 0) {
        printf("nanoaod::Reader: no n%s<i> branches in the tree\n", prefix);
        fOk = false;
      }
      if (!fOk)
        printf("nanoaod::Reader: required branches are missing - the tree is not "
               "the production this reader was written for\n");
    }

    bool Ok() const { return fOk; }
    int  NSlots() const { return fNSlots; }
    Long64_t GetEntries() const { return fTree->GetEntries(); }

    // Iteration.  Next() steps to the next entry of the range (the whole tree
    // if none was set) and returns false at its end; Entry() is the GLOBAL
    // chain entry number, which is what the a/b halves are taken from.
    bool Next() { return fReader.Next(); }
    Long64_t Entry() const { return fReader.GetCurrentEntry(); }
    bool SetEntry(Long64_t e) { return fReader.SetEntry(e) == TTreeReader::kEntryValid; }
    // [first, last); last < 0 means to the end, last <= first is an EMPTY
    // range and returns false.  Not ROOT's convention: TTreeReader takes any
    // end <= begin as "to the end" and says kEntryValid, so an empty job
    // [k, k) handed to it straight would read from k to the end of the chain.
    bool SetEntriesRange(Long64_t first, Long64_t last) {
      if (last >= 0 && last <= first) return false;
      return fReader.SetEntriesRange(first, last < 0 ? -1 : last) == TTreeReader::kEntryValid; }

    // ---- event ----
    UInt_t    Run()   const { return **fRun; }
    UInt_t    Lumi()  const { return **fLumi; }
    ULong64_t Event() const { return **fEvent; }
    int   NPU() const { return **fNPU; }
    float PVz() const { return **fPVz; }
    int   NPV() const { return **fNPV; }

    // ---- particle flow, the whole event ----
    unsigned NPF() const { return fPFpt->GetSize(); }
    const TTreeReaderArray<Float_t> &PFpt()   const { return *fPFpt; }
    const TTreeReaderArray<Float_t> &PFeta()  const { return *fPFeta; }
    const TTreeReaderArray<Float_t> &PFphi()  const { return *fPFphi; }
    const TTreeReaderArray<Float_t> &PFmass() const { return *fPFmass; }
    const TTreeReaderArray<Float_t> &PFdz()   const { return *fPFdz; }
    const TTreeReaderArray<Float_t> &PFpuppiWeight() const { return *fPFw; }
    const TTreeReaderArray<Int_t> &PFpdgId()   const { return *fPFpdg; }
    const TTreeReaderArray<Int_t> &PFcharge()  const { return *fPFq; }
    const TTreeReaderArray<Int_t> &PFvertexRef() const { return *fPFvref; }
    const TTreeReaderArray<Int_t> &PFfromPVvertexRef() const { return *fPFfpv; }
    const TTreeReaderArray<Int_t> &PFfromPV0() const { return *fPFfpv0; }
    const TTreeReaderArray<Int_t> &PFpvAssocQuality() const { return *fPFpaq; }

    // ---- generated interactions, per slot ----
    int   NGen(int s)  const { return (int)fGpt[s]->GetSize(); }
    float GenVtxZ(int s) const { return **fGvz[s]; }
    const TTreeReaderArray<Float_t> &GenPt(int s)   const { return *fGpt[s]; }
    const TTreeReaderArray<Float_t> &GenEta(int s)  const { return *fGeta[s]; }
    const TTreeReaderArray<Float_t> &GenPhi(int s)  const { return *fGphi[s]; }
    const TTreeReaderArray<Float_t> &GenMass(int s) const { return *fGmass[s]; }
    const TTreeReaderArray<Int_t>   &GenPdgId(int s) const { return *fGpdg[s]; }
    // Slots that carry an interaction this entry: 0 (hard scatter) .. nPU,
    // clipped to what the tree has.
    int NFilledSlots() const {
      const int n = NPU() + 1; return n < fNSlots ? n : fNSlots; }

    // The puppi::Event of this entry, as v4 LoadCandidates(): ztrk = PV_z + dz,
    // charged = charge != 0, q = its sign, the vertex keys as stored.  Does not
    // call puppi::Prepare (the analyzer owns the Config).
    void FillPF(puppi::Event &pf) const
    {
      pf.clear();
      const unsigned n = NPF();
      const float pvz = PVz();
      pf.pt.reserve(n); pf.eta.reserve(n); pf.phi.reserve(n); pf.mass.reserve(n);
      pf.ztrk.reserve(n); pf.charged.reserve(n); pf.pdg.reserve(n); pf.q.reserve(n);
      pf.vref.reserve(n); pf.fpvref.reserve(n); pf.paq.reserve(n);
      for (unsigned j = 0; j < n; ++j) {
        pf.pt  .push_back((*fPFpt)[j]);
        pf.eta .push_back((*fPFeta)[j]);
        pf.phi .push_back((*fPFphi)[j]);
        pf.mass.push_back((*fPFmass)[j]);
        pf.ztrk.push_back(pvz + (*fPFdz)[j]);
        const int c = (*fPFq)[j];
        pf.charged.push_back(c != 0 ? 1 : 0);
        pf.pdg   .push_back((*fPFpdg)[j]);
        pf.q     .push_back((signed char)(c > 0 ? 1 : (c < 0 ? -1 : 0)));
        pf.vref  .push_back((*fPFvref)[j]);
        pf.fpvref.push_back((char)(*fPFfpv)[j]);
        pf.paq   .push_back((char)(*fPFpaq)[j]);
      }
      pf.nPV = NPV();
    }

    // The generated particles of every filled slot, flattened in slot order
    // (interaction i is the contiguous run with inter == i), |eta| < etaMax,
    // neutrinos dropped if dropNu, z [cm] of the interaction, q the sign of the
    // charge.  vz gets one z per filled slot.  Returns the number of filled
    // slots (nPU + 1, clipped), as v4 ReadInteractions().
    int FillGen(genlink::GenList &gen, std::vector<double> &vz,
                double etaMax = 5.5, bool dropNu = true) const
    {
      gen.clear();
      const int nslot = NFilledSlots();
      vz.assign(nslot, 0.);
      for (int i = 0; i < nslot; ++i) {
        vz[i] = GenVtxZ(i);
        const int n = NGen(i);
        const TTreeReaderArray<Float_t> &pt = *fGpt[i], &eta = *fGeta[i],
                                        &phi = *fGphi[i], &mass = *fGmass[i];
        const TTreeReaderArray<Int_t> &pdg = *fGpdg[i];
        for (int j = 0; j < n; ++j) {
          const int id = pdg[j], a = std::abs(id);
          if (dropNu && (a == 12 || a == 14 || a == 16)) continue;
          const double e = eta[j];
          if (std::fabs(e) > etaMax) continue;
          gen.pt.push_back(pt[j]); gen.eta.push_back(e); gen.phi.push_back(phi[j]);
          gen.mass.push_back(mass[j]); gen.z.push_back(vz[i]);
          gen.pdg.push_back(id); gen.inter.push_back(i);
          gen.q.push_back((signed char)PdgCharge(id));
        }
      }
      return nslot;
    }

    // Sign of the charge of a PDG code from TDatabasePDG, cached: the lookup
    // is a string-keyed map inside ROOT and there are ~9000 particles an event.
    static int PdgCharge(int pdg)
    {
      static std::unordered_map<int,int> cache;
      auto it = cache.find(pdg);
      if (it != cache.end()) return it->second;
      int q = 0;
      if (TParticlePDG *p = TDatabasePDG::Instance()->GetParticle(pdg)) {
        const double c = p->Charge(); q = c > 0 ? 1 : (c < 0 ? -1 : 0);
      }
      cache[pdg] = q;
      return q;
    }

  private:
    template <class T> std::unique_ptr<TTreeReaderValue<T> > Val(const char *b) {
      Check(b); return std::make_unique<TTreeReaderValue<T> >(fReader, b); }
    template <class T> std::unique_ptr<TTreeReaderArray<T> > Arr(const char *b) {
      Check(b); return std::make_unique<TTreeReaderArray<T> >(fReader, b); }
    void Check(const char *b) {
      if (fTree->GetBranch(b)) return;
      printf("nanoaod::Reader: MISSING branch %s\n", b);
      fOk = false;
    }

    TTree      *fTree;
    TTreeReader fReader;
    bool fOk;
    int  fNSlots;
    std::unique_ptr<TTreeReaderValue<UInt_t> >    fRun, fLumi;
    std::unique_ptr<TTreeReaderValue<ULong64_t> > fEvent;
    std::unique_ptr<TTreeReaderValue<Int_t> >     fNPU;
    std::unique_ptr<TTreeReaderValue<Float_t> >   fPVz;
    std::unique_ptr<TTreeReaderValue<UChar_t> >   fNPV;
    std::unique_ptr<TTreeReaderArray<Float_t> > fPFpt, fPFeta, fPFphi, fPFmass, fPFdz, fPFw;
    std::unique_ptr<TTreeReaderArray<Int_t> >   fPFpdg, fPFq, fPFvref, fPFfpv, fPFfpv0, fPFpaq;
    std::vector< std::unique_ptr<TTreeReaderArray<Float_t> > > fGpt, fGeta, fGphi, fGmass;
    std::vector< std::unique_ptr<TTreeReaderArray<Int_t> > >   fGpdg;
    std::vector< std::unique_ptr<TTreeReaderValue<Float_t> > > fGvz;
  };

} // namespace nanoaod
#endif
