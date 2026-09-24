# genseed
Create particle-level jet (GenJet) seeds, i.e. genseed, from gen particles matched to reco particles and weighed with PUPPI weights, and to do the opposite and create reconstruction-level jets (RecoJets) out of reco particles matched to gen particles, i.e. recoseeds. Goal is to improve matching efficiencies, resolution and gaussianity.

**Document:** [`doc/genseed.pdf`](doc/genseed.pdf) (algorithm, every knob, tuple and
histogram schema, all results; generated from the plots in [`plots/`](plots/)).

## What it does

In high pileup a generated jet and the jet the detector reconstructs from the
same interaction are not the same object, and any pairing rule between the two
collections (a ΔR cut, a ghost assignment) leaves a response whose width and
tails describe the pairing as much as the detector. genseed links every
generated particle of every pileup interaction to the particle-flow
candidate(s) that carry its energy, and uses the links in both directions.
Four jet collections come out, all anti-kT R = 0.4 from 1 GeV:

| collection | built from | per |
|---|---|---|
| **pure GenJets** | the generated particles of one interaction | interaction |
| **PUPPI RecoJets** | PF candidates with the PUPPI weights of one vertex hypothesis, plus ghosts of that vertex's own unlinked generated particles; duplicates across vertices removed | reconstructed vertex |
| **GenJetSeed** (`link`) | for a reco jet: Σ over its PF constituents of w_i(v) · G_i, the generated 4-vector routed to candidate i by the linker, plus the own unlinked ghosts | reco jet |
| **RecoJetSeed** (`rs`) | for a pure GenJet: its linked PF candidates with weight s_iJ · w_i(v), clustered again; s_iJ is the share of i's linked generated pT that came from J, w_i(v) the PUPPI weight at the interaction's vertex | pure GenJet |

plus `rs0`, the same RecoJetSeed without PUPPI weights. pT^raw/pT^link is then
a detector response of the reconstructed jet (the pileup it kept is in the
seed at the weight PUPPI gave it), and pT^rs/pT^pure is what the detector and
PUPPI made of that generated jet, with no pairing at all.

The linker ([`genlink.h`](genlink.h)): pass 1 links tracks one to one (charge
sign, ΔR, pT ratio, dz); pass 2 shares calorimeter energy with Gaussian kernels
between compatible generated particles, at the calorimeter impact point
re-referenced to the primary vertex (parallax) and bent in the solenoid for
charged particles without a track. It is checked on data, not assumed: the
parallax slope is +0.979 ± 0.001 (1 expected), track links have a dz core of
0.014 cm and pT_PF/pT_gen 1.002 ± 0.016, and the orphan (unlinked) PF fraction
of a reconstructed jet is 1%.

## Quick start

```
# input: one NanoAOD path per line in files.txt, or LIST=<file.root>
./run.sh p1 8                                                    # pass 1, 8 jobs -> rootfiles/GenSeed_p1.root
root -l -b -q 'drawGenSeed.C+("rootfiles/GenSeed_p1.root","_v1")'  # -> text/jec_v1.txt (+ plots)
./run.sh p2 8                                                    # pass 2 with the correction -> rootfiles/GenSeed_v1.root
root -l -b -q 'drawGenSeed.C+("rootfiles/GenSeed_v1.root","_v1")'  # plots/, doc/plots.tex, doc/tables.tex
cd doc && pdflatex genseed && pdflatex genseed
```

A 200-event test: `root -l -b -q 'runGenSeed.C("files.txt","_test",0,1,200)'`.
Requirements: ROOT 6 (tested with 6.28/06); FastJet is optional. On an Apple
M2 the full 25 000-crossing sample takes 2.5 minutes on 8 cores (47 ms per
crossing, 10 ms of it the linker).

The input is a NanoAOD with the PF candidates (`PFCandV2_*`: pt, eta, phi,
mass, pdgId, charge, dz, vertexRef and its quality flags) and the generated
particles of every pileup interaction (`GenPartCandBX0PUEvent<i>_*`,
`GenVtxBX0PUEvent<i>_z`), as in `MC24NanoV15_PU_IT_OOT`.

## Files

| file | what |
|---|---|
| [`GenSeed.h`](GenSeed.h), [`GenSeed.C`](GenSeed.C) | the analyzer: collections, topology, bridge, tuple, histograms; the header documents the algorithm and every knob |
| [`runGenSeed.C`](runGenSeed.C), [`run.sh`](run.sh) | driver (knobs, job split by entry range, backend, build dir) and the two-pass script |
| [`genlink.h`](genlink.h) | the generated-particle ↔ PF linker (ROOT-free) |
| [`puppi.h`](puppi.h) | per-vertex PUPPI, an independent implementation (ROOT-free) |
| [`tiledjet.h`](tiledjet.h) | clean-room tiled anti-kT/kT/CA after Cacciari–Salam, hep-ph/0512210 (no FastJet code) |
| [`bench_tiledjet.C`](bench_tiledjet.C) | exactness and speed of tiledjet against an N² reference and FastJet |
| [`nanoreader.h`](nanoreader.h), [`jec.h`](jec.h) | TTreeReader input layer; the reco-jet correction read in pass 2 |
| [`drawGenSeed.C`](drawGenSeed.C) | plots, the correction table and the generated LaTeX (`doc/plots.tex`, `doc/tables.tex`) |
| `doc/` | the document; `algorithm.tex`, `howto.tex`, `abstract.tex`, `benchmark.tex` by hand, `plots.tex`, `tables.tex` generated |

## Outputs

* `rootfiles/GenSeed_<tag>.root`, directory `hist/` (event counts, every knob
  in `hopts`, PUPPI version, correction file), `spec/` (spectra per |y| for the
  four collections, halves a and b), `eff/` (pure-jet efficiency for ΔR < 0.2,
  mutual dominance, linked, RecoJetSeed), `fake/`, `resp/` (response TH2 vs pT
  per |η| for reco/link, reco/linkz, reco/pure ΔR < 0.2, rs/pure, rs0/pure,
  reco/rs, per topology and per N_PU slice), `topo/` (gen- and reco-side
  topology 1↔1/split/merge/tangle/lost|none), `bridge/` (the seed → pure matrix
  with its fakes and misses), `map/` (η–φ occupancy and response maps),
  `link/` (linker diagnostics).
* the TTree `genseed` in the same file: per crossing the vertices, the pure
  GenJets (`gj_*`, with their RecoJetSeed `gj_rs_*` and `gj_rs0_pt`, topology,
  partners) and the kept reco jets (`rj_*`, with pT^corr, both seeds, the
  own/far/unlinked/orphan fractions, topology and partner). About 12 kB per
  crossing.

## Results (25 000 crossings, ⟨N_PU⟩ = 45, 14.2 μb⁻¹, pass 2)

Response at |η| < 0.5 (σ_core/μ of the Gaussian core; P(R > 1.5) the right tail):

| pT [GeV] | reco/link σc/μ | ΔR < 0.2 σc/μ | reco/link P(R>1.5) | ΔR < 0.2 P(R>1.5) | rs/pure median | rs0/pure median |
|---|---|---|---|---|---|---|
| 5–6 | **0.273** | 0.433 | **0.028** | 0.095 | 0.589 | 0.704 |
| 8–10 | **0.252** | 0.394 | **0.015** | 0.047 | 0.618 | 0.719 |
| 12–15 | **0.225** | 0.325 | **0.006** | 0.017 | 0.663 | 0.745 |
| 18–21 | **0.207** | 0.257 | **0.003** | 0.012 | 0.717 | 0.781 |

The GenJetSeed core holds the whole distribution (f_core 1.00–1.06), and its
median moves by at most 0.026 between N_PU slices below 30 GeV.

Pure GenJets with a partner, |η| < 2.5:

| pT [GeV] | ΔR < 0.2 | mutual dominance | linked | has a RecoJetSeed |
|---|---|---|---|---|
| 5–6 | 0.60 | 0.70 | 0.76 | 0.93 |
| 8–10 | 0.70 | 0.81 | 0.83 | 0.97 |
| 12–15 | 0.80 | 0.88 | 0.89 | 0.99 |
| 18–21 | 0.89 | 0.93 | 0.94 | 0.99 |

What the RecoJetSeed shows: at 5–6 GeV in the barrel the detector plus PUPPI
keep a median 59% of a generated jet's pT (70% without the PUPPI weights), and
the reconstructed partner carries 26% more than its RecoJetSeed (reco/rs 1.26),
the pileup that PUPPI kept and clustering confusion; reco/pure is the product.
Of the 41% lost, 30 points are the detector (generated particles below the
calorimeter thresholds, which have no PF partner, and the calorimeter response:
rs0/pure is already 0.70) and 11 points are PUPPI: 59% of the linked neutral pT
and 10% of the linked charged pT get weight 0 at the jet's own vertex (the
neutral floor, 0.18 + 0.013 N_PV GeV, and the 0.03 cm track–vertex
association). The response is filled for jets that have a RecoJetSeed and is
read together with that efficiency.

## Known limitations

* **Forward (|η| > 3).** Per-vertex PUPPI has no vertex handle without tracks,
  and its forward neutral floor (1.85 + 0.07 N_PV GeV per candidate) removes
  almost everything at this pileup: 3–6% as many reco jets as pure GenJets above
  5 GeV, linked efficiency ~1%. Outside the tracker the analysis has to be done
  per crossing, not per vertex. 2.5 < |η| < 3 is in between.
* The RecoJetSeed response has its tail on the left (it is built from a subset
  of the jet's own candidates, so R ≲ 1); the v4 core-fit protocol, which
  expects a right tail, often does not converge for it, and the quantile widths
  are the measure to use.
* Vertex ownership is by z proximity (the interaction nearest to a reco vertex
  within 0.05 cm owns it). When two interactions sit within 0.2 cm, the vertex
  can go to the neighbour rather than to the interaction whose jet made it:
  for 3.5% of the pure GenJets above 10 GeV, the vertex they are processed at
  is owned by another interaction. The seed itself is unaffected, but that
  jet's energy is booked as "far" and the own unlinked ghosts come from the
  wrong interaction. Ownership by the linked track pT² of each vertex (the
  PV definition) would fix it.
* The linker needs the generated particles of every pileup interaction. For a
  standard production with the signal generated separately, it runs on the
  signal's particles alone and everything unlinked counts as pileup (see
  *Towards JMENANO* in the document).

## Towards JMENANO

The per-event inputs are the PF candidate table with `dz` and `vertexRef`
(with its quality flags) and the generated particles with their interaction z.
[`genlink.h`](genlink.h), [`puppi.h`](puppi.h) and [`tiledjet.h`](tiledjet.h)
are plain C++17 with no ROOT dependency, so they can be called from a CMSSW
producer that writes the GenJetSeed and RecoJetSeed collections next to the
existing jet tables.

## Licence

Apache License 2.0 ([`LICENSE`](LICENSE)); see [`NOTICE`](NOTICE) for the
provenance of the clusterer and of PUPPI, and for `tdrstyle_mod22.C`, which
contains the CMS publication-style macros under their own terms. This package
grew out of the lowptjets v4 study of the low-pT inclusive jet spectrum in high
pileup (M. Voutilainen, HIP/CMS).
