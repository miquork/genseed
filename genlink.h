// genseed - Copyright 2026 Mikko Voutilainen - Apache License 2.0
#ifndef genseed_genlink_h
#define genseed_genlink_h
// genlink.h - the generated-particle <-> particle-flow linker: which PF
// candidate(s) carry the energy of each generated particle, once per event and
// independent of the vertex hypothesis.  genseed port of lowptjets/genlink.h
// (v4); the algorithm is unchanged, only the puppi.h it includes lost ROOT.
//
// WHY.  v3 seeded the generated jet from the gen particles of the interaction
// nearest the reconstructed vertex, weighted by a z window.  The pairing became
// complete (93% vs 68% at 5-6 GeV) but the response R = pT_raw/pT_seed grew a
// right tail that a dR < 0.2 pairing does not have: the seed did not know what
// PUPPI had kept in the reconstructed jet.  A neutral from a far interaction
// that PUPPI keeps is in the jet but not in the seed; a neutral of the own
// interaction that PUPPI floors is in the seed but not in the jet; and the
// 1/n_near sharing halved the own seed whenever two interactions sat inside
// +-0.2 cm, which at four interactions per cm is half the time.  v4 gives every
// gen particle a PF partner, so that it can enter the seed of vertex hypothesis
// v with the PUPPI weight w_i(v) OF ITS PARTNER: whatever PUPPI decides about
// candidate i, it decides about the gen energy behind it, and the seed knows.
//
// THE SEED STAYS AT GEN SCALE.  The routed 4-vector G_i of candidate i is the
// sum of f_gi p_g over its gen partners, at the GENERATED pT with the generated
// mass, never pT_i.  The seed of a jet is then sum_i w_i(v) G_i (plus the
// unlinked own ghosts, see LowPtJetsV4Fill), and R = pT_raw/pT_seed is a pure
// detector response: calorimeter non-linearity, thresholds, tracking
// inefficiency, the neutral floor.  Had the seed used pT_i, R would be 1 by
// construction for everything PUPPI kept - a tautology, not a measurement -
// and every PUPPI decision would have been hidden inside the seed instead of
// being measured in the bridge seed -> pure where it belongs.
//
// GEOMETRY BEFORE KERNELS.  A PF track's direction is measured at the track,
// so its target is (eta_g, phi_g) as generated.  A calorimeter deposit is not:
// PF builds the direction of a neutral candidate from the reconstructed PRIMARY
// vertex, while the gen particle left interaction j at z_j, and |z_j - PV_z| is
// ~6 cm rms.  Seen from PV_z the same crystal is at a different eta:
// d(eta) = (z_j - PV_z)/(r cosh eta), 0.045 rms in the barrel against an ECAL
// crystal of 0.0175 and a photon cone of 0.06.  Without the parallax term a
// quarter of the barrel photons from far interactions miss the cone outright -
// and those are exactly the kept-pileup photons the linker exists to find.  So
// every calorimeter target is the point (r_D, z_D) where the gen particle
// crosses surface D, re-referenced to z_ref = PV_z (parallaxRef = 1; 0 uses the
// origin, diagnostic hparallax decides).  A charged particle also bends:
// x = 0.2998 B r/(2 pT) is the half-chord over the helix radius, the crossing
// is at phi_g - bendSign q asin(x) (positive charge to smaller phi in the CMS
// +z field), the arc is longer than the chord so z advances by asin(x)/x, and
// at x >= 1 the particle never reaches the surface: a barrel charged hadron
// below 0.735 GeV never sees the ECAL, below 1.08 GeV never the HCAL (looper).
// In HF the field ends at the solenoid, so the bend accrues only up to
// r = zSolenoid/|sinh eta|.
//
// PASS 1, PROMPT TRACKS, ONE-TO-ONE.  A prompt charged gen particle and a PF
// track are the same object seen twice, so the pairing is a hard assignment,
// not a share.  Pairs inside dR < 0.03 (0.05 past |eta| 1.5) with the same
// charge sign, |ln(pT_i/pT_g)| < 0.7 and |z_i - z_g| inside max(0.1 cm, 3
// sigma_z), sigma_z = sqrt(0.003^2 + (0.02 cosh(eta)^1.5/pT)^2) cm, are scored by
// chi2 = (dR/0.008)^2 + (dz/sigma_z)^2 + (ln ratio/0.15)^2, sorted, and accepted
// greedily while both ends are free.  The charge sign is the cheapest and
// strongest discriminator available (half of every accidental is gone), dz at
// the real track resolution separates interactions 0.25 cm apart, and the
// greedy order settles the jet cores where several tracks share 0.03.  A track
// taken in pass 1 is claimed: no gen particle can share it later.
//
// PASS 2, CALORIMETER AND HF, TIERED SHARING.  A calorimeter cluster genuinely
// merges several particles and a particle genuinely spreads over clusters, so
// here the gen particle is PARTITIONED: kernel a_i = C exp(-dR^2/2 sigma^2) over
// the PF candidates inside a hard cone at the class-appropriate target, with a
// compatibility C per (gen class, PF class) and sigma the granularity of that
// PF type; f_gi = a_i/sum a, so a gen particle with any candidate is fully
// linked.  The TIER RULE (any candidate with C == 1 drops the C < 1 ones) is the
// 'track before calo, photon before hadron' prior in one line: a photon with a
// PF photon nearby flows into photons only, a K_L with a PF neutral hadron
// nearby into that only, and the cross terms (photon -> NH, NH -> photon,
// photon -> conversion track, lost track -> photon) are used only when the
// natural partner is missing.  An untracked charged hadron (no pass-1 link, not
// a looper) goes to the BENT calorimeter target, or to an unclaimed same-sign
// track within kinkDR and kinkDz (decay in flight).  A V0 spreads its decay
// products over min(0.35, 0.05 + 0.45/pT) - the 0.206 GeV decay momentum - onto
// unclaimed tracks of either sign (no dz: it is displaced), photons and neutral
// hadrons, along straight lines.  Beyond |eta| 2.9 everything goes to the HF
// candidates beyond 2.85, hadrons to HFH and photons/electrons to HFE at C = 1,
// crossed at 0.5.  Every sigma and cone is scaled by kernelScale, and by
// regionFactor = 1.4 at 1.5 < |eta| < 3, but a cone is then clamped at 0.4,
// the cell of the coarse grid: already at kernelScale = 1 that stops the V0
// cone below 1.9 GeV past |eta| 1.5 (0.35 x 1.4 = 0.49 at the most) and the HF
// cone at 2.9-3.0 (0.30 x 1.4 = 0.42), and the barrel V0 cone follows above
// kernelScale 1.14, so a stability test at a larger kernelScale widens the
// widest cones by less than the knob says.  A TRACKED charged hadron gets
// nothing further: its calorimeter excess stays an orphan, because the ratio
// of the track pT to the gen pT IS the single-particle response that R has to
// show.
//
// NO RESCUE PASS.  A PF candidate with no gen partner inside a physically sized
// cone is either the calorimeter excess of a tracked hadron, or a deposit the
// gen record does not contain (out-of-time pileup, noise, a floored fragment).
// Widening a cone until the nearest gen particle takes it would move detector
// response into the seed in the first case and pileup at full gen pT into the
// seed in the second - a narrower R manufactured by construction, which the
// pre-mortem flagged as the one failure mode that looks like success.  Orphans
// are counted per PF class and per N_PU (horphan) and stay in the reconstructed
// jet, where R measures them.
//
// The output is a CSR link list both ways (gen -> PF, PF -> gen), the routed
// 4-vector per PF, flags, and one LinkRec per accepted link for diagnostics.
// Cost: ~8000 gen x <= 3 targets x ~25-40 PF in a 3x3 block of the puppi::Grid,
// a few ms per event (39 ms measured in v4 with everything on).  Pure C++:
// puppi.h for Event and Grid, nothing of ROOT.
//
// WHAT THE DIAGNOSTICS MUST SHOW.  The linker has two geometric assumptions
// built in and both are checked by the LinkRec output, not assumed; the v4
// production on 25000 events gave the reference values and any port has to
// reproduce them before its seeds mean anything.  (1) hparallax: eta of the
// linked PF neutral against the eta of the parallax-corrected target, minus the
// generated eta, has slope +0.98 (+0.979 +- 0.001): PF neutrals really point
// from the primary vertex, and without the term a quarter of the barrel photons
// from far interactions would miss their cone.  (2) hbend: phi_PF - phi_gen
// against q/pT for untracked charged hadrons has slope -0.97 rad GeV (positive
// charge to smaller phi); it is steeper than the ECAL-barrel -0.74 because the
// targets are a mixture of the ECAL and HCAL surfaces.  (3) pass-1 track links
// have a dz core of 0.014 cm and a pT_PF/pT_gen peak at 1.002 with sigma 0.016.
// (4) The linked fraction inside the tracker: charged hadrons 90% at 0.5 GeV
// and 95-98% above 1 GeV; photons 56% at 0.5, 68% at 1, 93% at 5, 100% at 10;
// neutral hadrons 49% at 1, 60% at 2, 83% at 5, 94% at 10; HF 93-98% - the
// calorimeter thresholds and nothing else.  (5) The orphan fraction of the
// reconstructed PUPPI-weighted jet pT is 1-2% (by PF class at N_PU = 50: tracks
// 0.3%, photons 11%, neutral hadrons 0.3%).  A parallax slope away from 1, a
// bend slope of the wrong sign or an orphan fraction of 10% means a geometry
// knob (parallaxRef, bendSign) or a unit (z in cm here, m inside CaloTarget)
// is wrong, not that the detector changed.
#include "puppi.h"
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <type_traits>

namespace genlink {

  // Linking classes.  Gen: prompt charged (c tau > 1 m: pi, K, p, e, mu and any
  // other charged), V0 (K_S, Lambda, Sigma, Xi, Omega: they decay before the
  // tracker is done with them), photon, neutral hadron (K_L, n, other neutral).
  enum GenClass { kGTrk = 0, kGV0, kGGam, kGNH, kNGenClass };
  // PF: track (211, 11, 13), photon (22), neutral hadron (130), HF hadron (1),
  // HF em (2), other.
  enum PFClass  { kPTrk = 0, kPGam, kPNH, kPHFH, kPHFE, kPOther, kNPFClass };
  // Diagnostic class, for the linker histograms only.
  enum DiagClass { kDPich = 0, kDKp, kDPh, kDNh, kDV0, kDLep, kNDiagClass };

  inline bool IsV0(int a) {
    return a == 310 || a == 3122 || a == 3112 || a == 3222 || a == 3312 ||
           a == 3322 || a == 3334;
  }
  // -1 for a neutrino: nothing to link.
  inline int GenClassOf(int pdg, int q) {
    const int a = std::abs(pdg);
    if (a == 12 || a == 14 || a == 16) return -1;
    if (IsV0(a)) return kGV0;
    if (a == 22) return kGGam;
    if (q != 0) return kGTrk;
    return kGNH;
  }
  inline int DiagClassOf(int pdg) {
    const int a = std::abs(pdg);
    if (a == 211) return kDPich;
    if (a == 321 || a == 2212 || a > 1000000000) return kDKp;
    if (a == 22) return kDPh;
    if (a == 130 || a == 2112) return kDNh;
    if (IsV0(a)) return kDV0;
    if (a == 11 || a == 13) return kDLep;
    return kDNh;
  }
  inline int PFClassOf(int pdg, int q) {
    const int a = std::abs(pdg);
    if (a == 211 || a == 11 || a == 13) return kPTrk;
    if (a == 22)  return kPGam;
    if (a == 130) return kPNH;
    if (a == 1)   return kPHFH;
    if (a == 2)   return kPHFE;
    return q != 0 ? kPTrk : kPOther;
  }
  inline const char *GenClassTag(int c) {
    static const char *t[kNGenClass] = {"trk","v0","gam","nh"};
    return c >= 0 && c < kNGenClass ? t[c] : "none"; }
  inline const char *PFClassTag(int c) {
    static const char *t[kNPFClass] = {"trk","gam","nh","hfh","hfe","other"};
    return c >= 0 && c < kNPFClass ? t[c] : "none"; }
  inline const char *DiagTag(int c) {
    static const char *t[kNDiagClass] = {"pich","kp","ph","nh","v0","lep"};
    return c >= 0 && c < kNDiagClass ? t[c] : "none"; }

  struct Config {
    // geometry [m]
    double rEcal = 1.29, zEcal = 3.17, etaEcalBarrel = 1.479;
    double rHcal = 1.90, zHcal = 3.90, etaHcalBarrel = 1.30;
    double zHF = 11.15, zSolenoid = 3.4;          // HF bend accumulates to min(r, zSolenoid/|sinh eta|)
    double bField = 3.8;                          // T; x = 0.2998*B*r/(2 pT) = 0.5696 r/pT
    int    bendSign = +1;                         // phi_D = phi_g - bendSign*q*asin(x)
    int    parallaxRef = 1;                       // 1: directions seen from PV_z, 0: from z = 0
    // pass 1 (tracks)
    double dRTrk[2] = {0.03, 0.05};               // |eta|<1.5, >=1.5 ; hard cuts
    double sigRTrk = 0.008, sigLnPt = 0.15, lnPtMax = 0.7;
    double dzFloor[2] = {0.10, 0.30};             // cm, |eta|<2.5, 2.5-3.0
    double dzN = 0.003, dzSlope = 0.02;           // sigma_z = sqrt(dzN^2 + (dzSlope cosh(eta)^1.5 / pT)^2) cm
    int    useDz = 1, useCharge = 1;
    // pass 2 (calo/HF): sigma / hard cut / compatibility; sigma and cut scaled by kernelScale, x regionFactor
    // for 1.5<|eta|<3, and the cut then clamped at detail::kCell = 0.4 (see collect() in Link)
    double kernelScale = 1.0, regionFactor = 1.4;
    double etaHF = 2.9, etaHFpf = 2.85, sigHF = 0.12, cutHF = 0.30;
    double v0Cut0 = 0.05, v0Cut1 = 0.45, v0CutMax = 0.35, v0Sig = 0.08;
    double kinkDR = 0.05, kinkDz = 0.3;
  };

  // The generated particles of one event, all interactions.  Filled by the
  // filler in ReadInteractions; z [cm] is the z of the interaction the
  // particle came from, q the sign of its charge.  q is signed char for the
  // reason given at puppi::Event::q: as plain char on an unsigned-char
  // platform a pi- would carry q = 255 and CaloTarget would bend it by
  // 255 asin(x), so no untracked negative hadron would find its deposit.
  struct GenList {
    std::vector<double> pt, eta, phi, mass, z;
    std::vector<int> pdg, inter; std::vector<signed char> q;
    size_t size() const { return pt.size(); }
    void clear() { pt.clear(); eta.clear(); phi.clear(); mass.clear(); z.clear();
                   pdg.clear(); inter.clear(); q.clear(); }
  };
  static_assert(std::is_signed<decltype(GenList::q)::value_type>::value,
                "genlink::GenList::q holds -1 and must be signed");

  // One accepted link: dr to the target used, dz = z_i - z_g [cm] for pass 1
  // (0 otherwise), pass 1 or 2, and the eta of the target.
  struct LinkRec { int g, i; float f; char pass; float dr, dz; float etaTarget; };

  struct Links {
    // gen side
    std::vector<char>  linked;          // 1 if the particle has links (then sum_i f_gi == 1)
    std::vector<char>  looper;          // charged that never reaches the ECAL surface
    std::vector<int>   gfirst, gcount;  // CSR into gpf/gfrac
    std::vector<int>   gpf; std::vector<float> gfrac;
    // PF side
    std::vector<int>   pfirst, pcount;  // CSR into pgen/pfrac
    std::vector<int>   pgen; std::vector<float> pfrac;
    std::vector<double> Gpx, Gpy, Gpz, GE, Gpt;   // routed gen 4-vector per PF: sum_g f_gi p_g (gen scale, gen mass)
    std::vector<char>  orphan;          // PF with no gen link
    std::vector<LinkRec> recs;
    void clear() {
      linked.clear(); looper.clear(); gfirst.clear(); gcount.clear();
      gpf.clear(); gfrac.clear(); pfirst.clear(); pcount.clear();
      pgen.clear(); pfrac.clear(); Gpx.clear(); Gpy.clear(); Gpz.clear();
      GE.clear(); Gpt.clear(); orphan.clear(); recs.clear();
    }
  };

  namespace detail {
    enum Surface { kECAL = 0, kHCAL, kHF, kNSurface };
    // Two grids: a 3x3 block of the coarse one holds every cone below (largest
    // 0.35, 0.49 with the region factor, clamped to the cell), the fine one
    // serves the cones up to 0.2, which is most of them, at a quarter of the
    // pairs.  The clamp is what makes the coarse block sufficient, and it is
    // applied after kernelScale, so it also caps what that knob can widen.
    const double kCell = 0.4, kCellFine = 0.2;
    const double kEtaTrk = 3.0;   // pass 1 reaches this far

    inline double wrap(double p) {
      if (p > M_PI) p -= 2*M_PI; else if (p < -M_PI) p += 2*M_PI; return p; }

    // Where the particle (pt, eta, phi, q) from z = zg [m] crosses surface D,
    // as PF sees it from zref [m].  ok = 1 target, 0 no crossing (wrong side),
    // -1 looper (never reaches D); signed, or tE.ok < 0 never flags a looper
    // where char is unsigned.
    struct Target { double eta, phi; signed char ok; };
    inline Target CaloTarget(int D, double pt, double eta, double phi,
                             double zg, int q, double zref, const Config &c)
    {
      Target t; t.eta = t.phi = 0; t.ok = 0;
      const double sh = sinh(eta), ash = fabs(sh), sgn = eta < 0 ? -1. : 1.;
      const bool hf = (D == kHF);
      const double R  = (D == kECAL ? c.rEcal : c.rHcal);
      const double Z  = hf ? c.zHF : (D == kECAL ? c.zEcal : c.zHcal);
      const double eb = (D == kECAL ? c.etaEcalBarrel : c.etaHcalBarrel);
      const double k  = q ? 0.2998*c.bField/(2.*pt) : 0.;         // x = k r
      const double rb = hf ? c.zSolenoid/std::max(ash, 1e-9) : 1e9; // field ends at the solenoid
      double rD = 0, zD = 0, x = 0;
      bool barrel = !hf && fabs(eta) < eb;
      if (barrel) {
        rD = R; x = k*std::min(rD, rb);
        if (x >= 1.) { t.ok = -1; return t; }
        zD = zg + rD*sh*(x > 0 ? asin(x)/x : 1.);
        if (fabs(zD) > Z) barrel = false;
      }
      if (!barrel) {
        zD = sgn*Z;
        if (ash < 1e-9 || (zD - zg)*sgn <= 0.) return t;
        rD = (zD - zg)/sh;                          // straight line first
        // helix: the arc is longer than the chord, so the plane is reached at
        // a smaller radius; a few fixed-point steps settle it
        for (int it = 0; it < 3 && k > 0.; ++it) {
          x = k*std::min(rD, rb);
          if (x >= 1.) { t.ok = -1; return t; }
          rD = (zD - zg)/(sh*asin(x)/x);
        }
        x = k*std::min(rD, rb);
        if (x >= 1.) { t.ok = -1; return t; }
      }
      t.eta = asinh((zD - zref)/rD);
      t.phi = wrap(phi - c.bendSign*q*asin(x));
      t.ok = 1;
      return t;
    }

    struct Pair { float chi2; int g, i; float dr, dz; };
    struct Cand { int i; float a, dr, etaT, C; };
  }

  // pf: puppi::Event (needs pt, eta, phi, mass, ztrk, charged, pdg, q). pvz: PV_z [cm].
  inline void Link(const puppi::Event &pf, const GenList &gen, double pvz,
                   const Config &c, Links &out)
  {
    using namespace detail;
    const int ng = gen.size(), np = pf.size();
    out.clear();
    out.linked.assign(ng, 0); out.looper.assign(ng, 0);
    out.gfirst.assign(ng, 0); out.gcount.assign(ng, 0);
    out.pfirst.assign(np, 0); out.pcount.assign(np, 0);
    out.orphan.assign(np, 1);
    out.Gpx.assign(np, 0.); out.Gpy.assign(np, 0.); out.Gpz.assign(np, 0.);
    out.GE.assign(np, 0.);  out.Gpt.assign(np, 0.);
    out.recs.reserve(ng);

    // signed: gcls is -1 for a neutrino
    std::vector<signed char> pcls(np), gcls(ng), claimed(np, 0);
    for (int i = 0; i < np; ++i) pcls[i] = PFClassOf(pf.pdg[i], pf.q[i]);
    for (int g = 0; g < ng; ++g) gcls[g] = GenClassOf(gen.pdg[g], gen.q[g]);

    puppi::Grid grid, fine;
    grid.build(pf.eta, pf.phi, kCell); fine.build(pf.eta, pf.phi, kCellFine);
    const double zref = c.parallaxRef ? 0.01*pvz : 0.;   // m

    // ---- pass 1: prompt tracks, greedy one-to-one on chi2 -----------------
    std::vector<Pair> pairs; pairs.reserve(ng);
    for (int g = 0; g < ng; ++g) {
      if (gcls[g] != kGTrk) continue;
      const double eg = gen.eta[g], pg = gen.phi[g], ptg = gen.pt[g], zg = gen.z[g];
      const double aeta = fabs(eg);
      if (aeta > kEtaTrk) continue;
      const double cut2 = c.dRTrk[aeta < 1.5 ? 0 : 1]*c.dRTrk[aeta < 1.5 ? 0 : 1];
      const double ch = cosh(eg), s1 = c.dzSlope*ch*sqrt(ch)/ptg;
      const double sz = sqrt(c.dzN*c.dzN + s1*s1);
      const double dzMax = std::max(c.dzFloor[aeta < 2.5 ? 0 : 1], 3.*sz);
      fine.forNear(eg, pg, [&](int i){
        if (pcls[i] != kPTrk) return;
        if (c.useCharge && pf.q[i] != gen.q[g]) return;
        const double de = pf.eta[i]-eg, df = puppi::dphi(pf.phi[i], pg);
        const double dr2 = de*de + df*df;
        if (dr2 >= cut2) return;
        const double lr = log(pf.pt[i]/ptg);
        if (fabs(lr) >= c.lnPtMax) return;
        const double dz = pf.ztrk[i] - zg;
        double chi2 = dr2/(c.sigRTrk*c.sigRTrk) + lr*lr/(c.sigLnPt*c.sigLnPt);
        if (c.useDz) { if (fabs(dz) >= dzMax) return; chi2 += dz*dz/(sz*sz); }
        Pair p; p.chi2 = chi2; p.g = g; p.i = i; p.dr = sqrt(dr2); p.dz = dz;
        pairs.push_back(p);
      });
    }
    std::sort(pairs.begin(), pairs.end(),
              [](const Pair &a, const Pair &b){ return a.chi2 < b.chi2; });
    for (size_t k = 0; k < pairs.size(); ++k) {
      const Pair &p = pairs[k];
      if (out.linked[p.g] || claimed[p.i]) continue;
      out.linked[p.g] = 1; claimed[p.i] = 1;
      LinkRec r; r.g = p.g; r.i = p.i; r.f = 1.f; r.pass = 1;
      r.dr = p.dr; r.dz = p.dz; r.etaTarget = gen.eta[p.g];
      out.recs.push_back(r);
    }

    // ---- pass 2: calorimeter and HF, tiered sharing per unlinked gen -------
    std::vector<Cand> cand; cand.reserve(64);
    for (int g = 0; g < ng; ++g) {
      const int cls = gcls[g];
      if (cls < 0) continue;
      const double eg = gen.eta[g], pg = gen.phi[g], ptg = gen.pt[g];
      const double zg = gen.z[g], zgm = 0.01*zg;
      const int q = gen.q[g], apdg = std::abs(gen.pdg[g]);
      const double aeta = fabs(eg);
      // looper flag for every charged particle, tracked or not
      Target tE = {0, 0, 0}; bool haveE = false;
      if (q) { tE = CaloTarget(kECAL, ptg, eg, pg, zgm, q, zref, c); haveE = true;
               if (tE.ok < 0) out.looper[g] = 1; }
      if (out.linked[g]) continue;                 // tracked: nothing further
      const double ks = c.kernelScale*((aeta > 1.5 && aeta < 3.0) ? c.regionFactor : 1.);
      cand.clear();
      // mode 0: any PF of that class; 1: unclaimed track, either sign, no dz;
      // 2: unclaimed same-sign track within kinkDz (decay in flight);
      // 3: unclaimed track that is not a muon (conversion)
      auto collect = [&](const Target &t, int pc, double sig, double cut, double C,
                         int mode, bool hfOnly) {
        if (t.ok != 1 || C <= 0.) return;
        sig *= ks; cut = std::min(cut*ks, kCell);
        const double cut2 = cut*cut, s2 = 2.*sig*sig;
        (cut <= kCellFine ? fine : grid).forNear(t.eta, t.phi, [&](int i){
          if (pcls[i] != pc) return;
          if (hfOnly && fabs(pf.eta[i]) <= c.etaHFpf) return;
          if (mode) {
            if (claimed[i]) return;
            if (mode == 2 && (pf.q[i] != q || fabs(pf.ztrk[i]-zg) >= c.kinkDz)) return;
            if (mode == 3 && std::abs(pf.pdg[i]) == 13) return;
          }
          const double de = pf.eta[i]-t.eta, df = puppi::dphi(pf.phi[i], t.phi);
          const double dr2 = de*de + df*df;
          if (dr2 >= cut2) return;
          Cand k; k.i = i; k.a = C*exp(-dr2/s2); k.dr = sqrt(dr2);
          k.etaT = t.eta; k.C = C;
          cand.push_back(k);
        });
      };
      const Target tTrk = {eg, pg, 1};
      switch (cls) {
      case kGGam: {
        const Target e = CaloTarget(kECAL, ptg, eg, pg, zgm, 0, zref, c);
        const Target h = CaloTarget(kHCAL, ptg, eg, pg, zgm, 0, zref, c);
        collect(e, kPGam, 0.025, 0.06, 1.0, 0, false);
        collect(h, kPNH,  0.06,  0.15, 0.3, 0, false);
        collect(tTrk, kPTrk, 0.05, 0.10, 0.3, 3, false);   // conversion legs
        break; }
      case kGNH: {
        const Target h = CaloTarget(kHCAL, ptg, eg, pg, zgm, 0, zref, c);
        const Target e = CaloTarget(kECAL, ptg, eg, pg, zgm, 0, zref, c);
        collect(h, kPNH,  0.08, 0.20, 1.0, 0, false);
        collect(e, kPGam, 0.04, 0.10, 0.3, 0, false);
        break; }
      case kGTrk: {                                 // untracked: bent targets
        if (!haveE) tE = CaloTarget(kECAL, ptg, eg, pg, zgm, q, zref, c);
        const Target h = CaloTarget(kHCAL, ptg, eg, pg, zgm, q, zref, c);
        collect(tE, kPGam, 0.04, 0.10, 0.5, 0, false);
        collect(h,  kPNH,  0.08, 0.20, 1.0, 0, false);
        collect(tTrk, kPTrk, 0.5*c.kinkDR, c.kinkDR, 0.5, 2, false);
        break; }
      case kGV0: {                                  // straight lines, wide
        const double vcut = std::min(c.v0CutMax, c.v0Cut0 + c.v0Cut1/ptg);
        const Target e = CaloTarget(kECAL, ptg, eg, pg, zgm, 0, zref, c);
        const Target h = CaloTarget(kHCAL, ptg, eg, pg, zgm, 0, zref, c);
        collect(tTrk, kPTrk, c.v0Sig, vcut, 1.0, 1, false);
        collect(e, kPGam, c.v0Sig, vcut, 0.7, 0, false);
        collect(h, kPNH,  c.v0Sig, vcut, 1.0, 0, false);
        break; }
      }
      if (aeta > c.etaHF) {
        const Target f = CaloTarget(kHF, ptg, eg, pg, zgm, q, zref, c);
        const bool em = (apdg == 22 || apdg == 11);
        collect(f, kPHFH, c.sigHF, c.cutHF, em ? 0.5 : 1.0, 0, true);
        collect(f, kPHFE, c.sigHF, c.cutHF, em ? 1.0 : 0.5, 0, true);
      }
      if (cand.empty()) continue;
      // tier rule: a natural partner present -> the cross terms go
      bool tier1 = false;
      for (size_t k = 0; k < cand.size(); ++k) if (cand[k].C >= 1.f) { tier1 = true; break; }
      double sum = 0;
      for (size_t k = 0; k < cand.size(); ++k) {
        if (tier1 && cand[k].C < 1.f) cand[k].a = 0;
        sum += cand[k].a;
      }
      if (sum <= 0) continue;
      out.linked[g] = 1;
      for (size_t k = 0; k < cand.size(); ++k) {
        if (cand[k].a <= 0) continue;
        LinkRec r; r.g = g; r.i = cand[k].i; r.f = cand[k].a/sum; r.pass = 2;
        r.dr = cand[k].dr; r.dz = 0; r.etaTarget = cand[k].etaT;
        out.recs.push_back(r);
      }
    }

    // ---- freeze: CSR both ways and the routed gen 4-vectors ---------------
    const int nl = out.recs.size();
    for (int k = 0; k < nl; ++k) { ++out.gcount[out.recs[k].g]; ++out.pcount[out.recs[k].i]; }
    for (int g = 1; g < ng; ++g) out.gfirst[g] = out.gfirst[g-1] + out.gcount[g-1];
    for (int i = 1; i < np; ++i) out.pfirst[i] = out.pfirst[i-1] + out.pcount[i-1];
    out.gpf.resize(nl); out.gfrac.resize(nl); out.pgen.resize(nl); out.pfrac.resize(nl);
    std::vector<int> gpos(out.gfirst), ppos(out.pfirst);
    for (int k = 0; k < nl; ++k) {
      const LinkRec &r = out.recs[k];
      out.gpf[gpos[r.g]] = r.i;  out.gfrac[gpos[r.g]++] = r.f;
      out.pgen[ppos[r.i]] = r.g; out.pfrac[ppos[r.i]++] = r.f;
      const double pt = gen.pt[r.g], m = gen.mass[r.g];
      const double px = pt*cos(gen.phi[r.g]), py = pt*sin(gen.phi[r.g]);
      const double pz = pt*sinh(gen.eta[r.g]);
      const double E  = sqrt(px*px + py*py + pz*pz + m*m);
      out.Gpx[r.i] += r.f*px; out.Gpy[r.i] += r.f*py;
      out.Gpz[r.i] += r.f*pz; out.GE[r.i]  += r.f*E;
    }
    for (int i = 0; i < np; ++i) {
      out.Gpt[i] = sqrt(out.Gpx[i]*out.Gpx[i] + out.Gpy[i]*out.Gpy[i]);
      out.orphan[i] = out.pcount[i] == 0;
    }
  }

} // namespace genlink
#endif
