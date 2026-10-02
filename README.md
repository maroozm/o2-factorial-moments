# Factorial Moments (ALICE O2Physics)

O2Physics task for the **normalized factorial moments** analysis of particle
multiplicity fluctuations (ULM method), `o2-analysis-cf-factorial-moments`.

Reference: R. C. Hwa and C. B. Yang, *Phys. Rev. C* **85**, 044914 (2012),
[doi:10.1103/PhysRevC.85.044914](https://doi.org/10.1103/PhysRevC.85.044914).

```bibtex
@article{PhysRevC.85.044914,
  title = {Local multiplicity fluctuations as a signature of critical hadronization in heavy-ion collisions at TeV energies},
  author = {Hwa, Rudolph C. and Yang, C. B.},
  journal = {Phys. Rev. C},
  volume = {85},
  issue = {4},
  pages = {044914},
  numpages = {11},
  year = {2012},
  month = {Apr},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevC.85.044914},
  url = {https://link.aps.org/doi/10.1103/PhysRevC.85.044914}
}
```

This repository contains only the factorial-moment code, the JSON
configuration, the run/merge/plot scripts and the documentation needed to
produce `F_q` results from AO2Ds on CVMFS.

## Repository layout

| path | purpose |
|---|---|
| `src/FactorialMomentsTask.cxx` | the task (vendored from the CVMFS O2Physics tag) |
| `src/CMakeLists.txt`, `CMakeLists.txt`, `cmake/` | build definition (single workflow) |
| `Common/Core/RecoDecay.h` | small helper header used by the task |
| `enable` | CVMFS environment setup (O2/O2Physics/FairRoot/ROOT) |
| `build_cmd` | cmake + make + install wrapper |
| `run_mc.sh` | run the full chain on every AO2D of a file list |
| `mc_chain_config.json`, `factorial_mc.json` | task + chain configuration |
| `checkFMo2.C` | draw `ln(F_q)` vs `ln(M^2)` and dump `.dat` files |
| `validateFMo2.C` | independent cross-check of those `.dat` files |
| `cvmfs/` | run the chain inside a CVMFS singularity container |
| `data/545312/file.list` | AO2D input list (the `.root` files are **not** in git) |
| `results/` | run outputs (not in git) |

## 1. Environment (CVMFS)

All commands below run on a machine that mounts
`/cvmfs/alice.cern.ch` (lxplus, swrep, ...) and start with:

```
source ./enable
```

`enable` loads the ALICE software stack from CVMFS. Steering variables:

| variable | meaning (default) |
|---|---|
| `O2OPENACCESS_SW_TAG_CVMFS` | O2Physics CVMFS tag (`daily-20260929-0000-1`) |
| `O2OPENACCESS_USE_LOCAL` | set to use a locally installed stack instead of CVMFS |
| `O2OPENACCESS_SW_TAG_LOCAL` | tag of that local stack (`latest-o2physics-o2`) |

Without `O2OPENACCESS_USE_LOCAL` the CVMFS stack is used, i.e.

```
/cvmfs/alice.cern.ch/el9-x86_64/Packages/O2Physics/<tag>
/cvmfs/alice.cern.ch/el9-x86_64/Packages/O2/<tag>
```

`build_cmd`, `run_mc.sh` and `cvmfs/run_cvmfs_analysis` call `enable`
themselves, so a manual `source ./enable` is only needed for ad-hoc commands
(`root`, `hadd`, `bin/o2-analysis-cf-factorial-moments`, ...).

## 2. Building

```
./build_cmd          # cmake + make -j + make install
./build_cmd debug    # additionally enable verbose cmake tracing
```

This produces:

* `bin/o2-analysis-cf-factorial-moments` - the workflow executable
* `share/dpl/o2-analysis-cf-factorial-moments.json` - its DPL options

## 3. Configuration

Two JSON files are provided (both must be passed to *every* device, see below):

* `factorial_mc.json` - the options of `factorial-moments-task` only.
* `mc_chain_config.json` - the same task options plus the process switches of
  the helper devices (`propagation-service`, `mult-cent-table`, ...). This is
  the one used by `run_mc.sh`.

Main task options:

| option | default | meaning |
|---|---|---|
| `numPt` | 5 | number of pT bins (must match `ptCuts`) |
| `ptCuts` | `0.2, 2` | bin edges in GeV/c, `2*numPt` values |
| `samplesize` | 100 | events per subsample (defines `F_q` samples) |
| `centLimits` | `0, 5` | accepted centrality range (FT0C, Run 3) |
| `centralEta` | 0.9 | `|eta|` cut |
| `vertexXYZ` | `0.3, 0.4, 10` | vertex x, y (cm) and z (cm) cuts |
| `dcaXY`, `dcaZ` | 0.1, 1.0 | echoed into `metaConfig` only, **not applied** (see note) |
| `useMC` | false | echoed into `metaConfig`; MC handling follows the process switches |
| `smearPhi` | false | randomise track phi (Gaussian, sigma=2pi) before filling the eta-phi lattices (the shipped JSONs set it to `true`) |

Note on DCA: the cut that is actually applied to reconstructed tracks is the
hard-coded ITS parameterisation `|dca_xy| < 0.0105 + 0.035 / pT^1.1`
(`kDcaXY0/1/2`); there is no DCA_z cut. The following options are accepted (the
shipped JSON files set them) but are currently **not applied by any cut** - they
are kept only so the configuration files keep parsing: `cfgCutTpcChi2NCl`,
`cfgCutItsChi2NCl`, `cfgITScluster`, `cfgTPCcluster`, `cfgTPCnCrossedRows`,
`cfgTPCnCrossedRowsOverFindableCls`, `isApplyVertexTOFmatched`,
`isApplyVertexTRDmatched`, `isApplyExtraCorrCut`, `isApplyExtraPhiCut`,
`includeGlobalTracks`, `includeTPCTracks`, `includeITSTracks`, `useGlobalTrack`,
`reduceOutput`.

## 4. Running on CVMFS (lxplus and friends)

Input AO2Ds: `data/545312/` (`AO2D_001.root`, `AO2D_002.root`, `AO2D_003.root`,
`AO2D_005.root`, listed in `data/545312/file.list`) - ALICE Run 3 MC,
LHC25f3 / 545312. The multi-GB `.root` files are downloaded separately and are
git-ignored.

```
./run_mc.sh                            # every file of data/545312/file.list
./run_mc.sh data/545312/AO2D_002.root  # a single AO2D
```

Per AO2D it writes `results/<AO2D_NAME>/AnalysisResults.root` and `run.log`,
executing:

```
o2-analysis-event-selection-service -b $JSON $AOD $SHM --evselOpts.isMC 1 |
o2-analysis-propagationservice      -b $JSON $AOD $SHM |
o2-analysis-trackselection          -b $JSON $AOD $SHM |
o2-analysis-multcenttable           -b $JSON $AOD $SHM |
o2-analysis-cf-factorial-moments    -b $JSON $AOD $SHM > run.log 2>&1
```

with `JSON=--configuration json://$PWD/mc_chain_config.json`,
`AOD=--aod-file <AO2D>` and `SHM="--shm-segment-size 8589934592"`.

Gotchas learned the hard way:

* `--configuration` **and** `--aod-file` must be passed to *every* device in the
  pipe, otherwise the reader starts with an empty file name.
* Device process switches (`processMonteCarlo`, `processCentralityRun3`, ...) can
  only be set in the JSON; a CLI bool flag is always interpreted as `true`.
* `--shm-segment-size` is in **bytes**; the default is too small for the MC
  tables (the McParticles message is ~55 MB).
* The effective configuration of a run is dumped to `dpl-config.json` in the
  working directory - use it to check that the JSON was applied.
* Opening files through `/eos/...` prints a benign
  `TNetXNGFile::Open ... Permission denied` error; the file is opened locally
  anyway.

### Running inside a CVMFS singularity container

On machines without a CVMFS mount, see `cvmfs/`:

* `cvmfs/run_cvmfs_analysis` - the same chain, driven by `data/545312/file.list`
  (or `@<list>` as first argument)
* `cvmfs/cvmfs_cmd` - runs the script above inside the container
* `cvmfs/cvmfs2go` - the container/`cvmfs2go` helper

## 5. Merging and plotting

```
mkdir -p results/merged
hadd -f results/merged/AnalysisResults.root results/AO2D_*/AnalysisResults.root
root -l -q -e '.x checkFMo2.C+("results/merged/AnalysisResults.root", "./results")'
```

`checkFMo2.C(<rootfile>, <outdir>)` needs nothing but the root file (the output
directory defaults to `./results` and is created when missing) and detects
everything else from the file itself, then draws **every** pT bin it finds:

| what | where it comes from |
|---|---|
| number of pT bins | existing `mFinalFq2Sampled_bin<N>`, `N = 1, 2, ...` |
| pT ranges | `metaConfig` (`ptCuts=`), else the histogram titles |
| number of M bins | `GetNbinsX()` of `mFinalFq2Sampled_bin1` |
| M values | `GetNbinsX()` of `bin1/Reset/mEtaPhi<i>` (= `binningM[iM]`) |
| Fq orders | existing `mFinalFq<q>Sampled_bin1`, `q = 2, 3, ...` |
| centrality range | `metaConfig` (`centLimits=`), else the filled range of `mCentFT0C`, `mCentFT0M`, `mCentFV0A`, `mCentFT0A` |
| eta / z-vertex cuts | `metaConfig` (`centralEta=`, `vertexXYZ=`), else the filled range of `mEta` / `mVertexZ` |
| number of samples | `GetEntries() / (number of M bins)` |
| error bars | `mFqSum`/`mFqSq` when present, else `mFqError` (single job only) |

For every bin it writes `<outdir>/FactorialMoments_bin<N>.dat` and `.pdf`.
All header lines of the `.dat` start with `#` (source file, pT range,
centrality, subsamples, eta, z-vertex, column names), followed by the columns
`M2 lnM2 lnF2 errF2 ... lnF<qmax> errF<qmax>` - one pair per order found in the
file (`q = 2..7` here), where `err` = `sigma(Fq)/mean(Fq)` = the error on
`ln(Fq)`; points with `Fq <= 0` are written as `-999 -999`.

`validateFMo2.C(<outdir>)` (default `./results`) cross-checks those `.dat`
files against an independent pooled recomputation from the four per-file
results (expected agreement: ~1e-7 on the value, ~1e-5 on the error, i.e. the
float storage level):

```
root -l -q -e '.x validateFMo2.C+("./results")'
```

## 6. Output histograms and merging

Written by device `factorial-moments-task`, directory `factorial-moments-task`,
`q` = 2..7, `N` = pT bin (`bin1..bin3`), 52 M bins (M = 4, 6, ..., 106):

* `mFinalFq{q}_bin{N}`, `mFinalFq{q}Sampled_bin{N}`, `mFinalAvBin*` - filled with
  `Fill(iM, value)` and never reset, so `hadd` adds them. Value =
  `content / (GetEntries()/52)`, where `GetEntries()/52` = number of samples
  (samplesize=3 events each) = 6+8+8+7 = 29 for the four AO2Ds.
* `mFqSum{q}_bin{N}` / `mFqSq{q}_bin{N}` (`TH1D`) - sum and sum of squares of the
  same per-sample Fq. Merge-safe, so the mean **and its error** can be rebuilt
  after `hadd` of any number of jobs:

  ```
  N    = entries / 52
  mean = sum / N
  var  = (sumSq - N * mean^2) / (N - 1)
  err  = sqrt(var / N)
  ```

* `mFqError{q}_bin{N}` - standard error of **one** job (`SetBinContent`, NOT
  merge-safe). Corrected to the formula above (the original
  `TH1::GetStdDev()/sqrt(n)` was biased by sqrt((n-1)/n) and returned 0 whenever
  a value fell outside the 0..10 range of its scratch histogram).
* `metaConfig` - a 1-bin `TH1D` whose **title** carries the effective job
  configuration (see the `metaConfig` section below). Histogram titles are taken
  from the first input of `hadd`, so the merged file stays self-describing.
* `mEventSelected` - event cut-flow, 10 labelled bins (`all`, `sel8`, the
  border/pileup/ITS/z-vertex bits, `centrality`, `accepted`). `processRun3`
  fills the border and `kIsGoodITSLayersAll` steps, `processMCRec` fills
  `kIsVertexITSTPC` instead; steps a mode never applies stay empty.

## 7. metaConfig

`init()` registers a 1-bin histogram `metaConfig` whose **title** carries the
effective configuration of the job (no `;` in it, otherwise ROOT would split the
string into axis titles):

```
FMtask centLimits=0,5 numPt=3 centralEta=0.9 samplesize=3 ptMin=0.2 dcaXY=2.4 dcaZ=2 useMC=1 nfqOrder=6 smearPhi=1 vertexXYZ=0.3,0.4,10 ptCuts=0.2,2,0.4,2,0.4,1
```

Histogram titles are not summed by `hadd`, so the merged file keeps this string.
Files produced before `metaConfig` existed fall back to the ranges observed in
the QA histograms - for a truncated distribution that reports the *observed*
range (e.g. `1-5 %` instead of the configured `0-5 %`), so re-running the chain
is recommended once after upgrading the task.

## 8. Implementation notes (kept out of the source, per house style)

* The `M` binning is `binningM[iM] = 2 * (iM + 2)`, i.e. M = 4, 6, ..., 106
  (52 bins). Every sample fills every one of the 52 bins exactly once, so
  `nSamples = GetEntries() / 52`. `checkFMo2.C` reads the same M values back
  from the `bin<N>/Reset/mEtaPhi<i>` binning instead of hardcoding them.
* The factorial numerator uses the falling factorial
  `n!/(n-q)! = n*(n-1)*...*(n-q+1)` evaluated with an explicit product loop:
  `TMath::Factorial` overflows to `inf` for n >= 171.
* `checkpT()` can replace the track azimuth with `gRandom->Gaus(phi, 2*pi)`
  wrapped into `[0, 2*pi)` before filling the eta-phi lattices, i.e. the phi used
  for the factorial moments is effectively randomised. The switch is the
  `smearPhi` option (default `false`, recorded in `metaConfig`); the shipped
  JSON files set it to `true` and the validated `.dat` files were produced with
  it on.
* The eta-phi lattices are cleared lazily at the start of an event: only the
  pT bins that saw tracks in the previous event are `Reset()`, which is
  equivalent to clearing all of them but avoids zeroing ~5 MB of bins per event.
* `fqEvent` and `binConEvent` are accumulated unconditionally (not only when a
  sample is flushed) so both sums always run over exactly the same set of
  events.
* `TH1::GetStdDev()` is the population sigma (divides by n, not n-1) and
  silently ignores values outside the histogram range - it must not be used to
  build a standard error; use `mFqSum`/`mFqSq` instead.
* `TH1::SetTitle("text;x")` splits on the first `;`, so `GetTitle()` returns
  only `"text"`; `TSubString::Data()` is *not* null-terminated, copy into a
  `TString` before using it as a C string.

## License

GPLv3 - see [LICENSE](LICENSE).
