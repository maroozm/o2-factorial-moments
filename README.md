<p align="center">
  <img src="logo.png" alt="O2 CVMFS — Factorial Moments" width="216">
</p>

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
| `slotpool.C` | pool slots into `K`-event subsamples after `hadd`; writes `.dat` + CSV |
| `plotFqvsM.C` | draw `F_q` vs `M^2` straight from the slotpool CSV |
| `findK.C` | scan candidate `K` values (group count, tail, overshoot) for a merged file |
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
| `sampleSize` | 100 | events per subsample (defines `F_q` samples) |
| `centLimits` | `0, 5` | accepted centrality range (FT0C, Run 3) |
| `centralEta` | 0.9 | `|eta|` cut |
| `vertexXYZ` | `0.3, 0.4, 10` | vertex x, y (cm) and z (cm) cuts |
| `dcaXY`, `dcaZ` | 0.1, 1.0 | echoed into `metaConfig` only, **not applied** (see note) |
| `useMC` | false | echoed into `metaConfig`; MC handling follows the process switches |
| `smearPhi` | false | randomise track phi (Gaussian, sigma=2pi) before filling the eta-phi lattices (the shipped JSONs set it to `true`) |
| `flushTailSample` | true | flush a partial group of events as an extra subsample at end of stream (recorded in `metaConfig` as `flushTail=`) |
| `nSlots` | 1250 | hash buckets for merge-safe subsampling, `slot = splitmix64(globalBC) % nSlots` (recorded in `metaConfig` as `nSlots=`). **Must be identical in every job of a campaign**, see below |

Note on DCA: the cut that is actually applied to reconstructed tracks is the
hard-coded ITS parameterisation `|dca_xy| < 0.0105 + 0.035 / pT^1.1`
(`kDcaXY0/1/2`); there is no DCA_z cut. The following options are accepted (the
shipped JSON files set them) but are currently **not applied by any cut** - they
are kept only so the configuration files keep parsing: `maxTPCChi2NCl`,
`maxITSChi2NCl`, `minITSClusters`, `minTPCClusters`, `minTPCCrossedRows`,
`minTPCCrossedRowsOverFindableCls`, `evSelVertexTOFmatched`,
`evSelVertexTRDmatched`, `applyExtraCorrCut`, `applyExtraPhiCut`,
`includeGlobalTracks`, `includeTPCTracks`, `includeITSTracks`, `useGlobalTrack`,
`reduceOutput`.

## 4. Running on CVMFS (lxplus and friends)

Input AO2Ds: `data/545312/` - ALICE Run 3 MC, LHC25f3 / 545312, 49 files
`AO2D_001.root` ... `AO2D_050.root` listed in `data/545312/file.list` (there is
no `AO2D_004`: the production's `AOD/` directory skips `004/` as well). The
multi-GB `.root` files are downloaded separately from
`alien:///alice/sim/2025/LHC25f3/0/545312/AOD/<NNN>/AO2D.root` and are
git-ignored.

```
alien.py cp alien:///alice/sim/2025/LHC25f3/0/545312/AOD/007/AO2D.root \
           file://$PWD/data/545312/AO2D_007.root
```

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
* An AO2D read can fail **transiently**: the reader has been seen to die with
  `Couldn't get TTree "DF_<id>/O2track" from ...` plus
  `websocket_callback: Error while reading from websocket`, and to exit 128 with
  no `AnalysisResults.root` at all - even though the file is byte-identical to
  the grid and opens cleanly (checked DF by DF). It has passed unchanged on
  immediate retry, so treat a job that produced no output as retryable before
  suspecting the file. Note that no AO2D in this production carries an `O2track`
  tree at all, only `O2track_iu`, so that message is not by itself diagnostic of
  a bad file. If one file fails, `run_mc.sh` stops with it: re-run just that
  file (`./run_mc.sh data/545312/AO2D_<N>.root`) and confirm the output by
  checking `sum(slots/nEv) == mEventSelected[accepted]`.

### Running inside a CVMFS singularity container

On machines without a CVMFS mount, see `cvmfs/`:

* `cvmfs/run_cvmfs_analysis` - the same chain, driven by `data/545312/file.list`
  (or `@<list>` as first argument)
* `cvmfs/cvmfs_cmd` - runs the script above inside the container
* `cvmfs/cvmfs2go` - the container/`cvmfs2go` helper

## 5. Merging, pooling and plotting

```
mkdir -p results/merged
hadd -f results/merged/AnalysisResults.root results/AO2D_*/AnalysisResults.root
root -l -b -q 'slotpool.C+("results/merged/AnalysisResults.root", 17, "results/merged/AnalysisResults.root", "slotpool_K17.csv")'
root -l -b -q 'plotFqvsM.C+("slotpool_K17.csv")'
```

Everything after the merge is post-processing of the merged
`AnalysisResults.root` (Hyperloop delivers exactly **one** merged file and
`hadd` cannot split anything back apart, so subsamples of a chosen size `K`
must be reconstructable from the merge-safe slot output, see below):

* `slotpool.C(<file>, <K>, <baseline>, <csv>, <keepTail>, <datDir>)` pools the
  slots into `K`-event subsamples and writes the pooled `F_q` - values only,
  nothing of the old path is ever written to a file: `FactorialMoments_bin<N>.dat`
  into `datDir` (default `dat_slots/`) in the established layout, and a CSV
  when `csv` is given (`""` = no CSV). With a `baseline` file it additionally
  prints a console sanity table against the old `sampleSize` path; the merged
  file itself serves as baseline, since the old-path objects live in it too.
* `plotFqvsM.C(<csv>, <ptSel>, <logMode>, <ymin>, <ymax>, <outDir>)` draws
  `F_q` vs `M^2` directly from that CSV (it never reads a `.dat`): one PDF per
  pT bin with every order as a marker+band series, established ATLAS style -
  details below.
* `findK.C(<file>, <kMin>, <kMax>, <csv>)` scans candidate `K` values so the
  subsample size can be chosen from the data - details below.

The CSV starts with a `#` header block (source path, date, K/keepTail/
subsamples, slot and event counts, centrality, eta and vertex cuts, pT bin
ranges - enough to label a figure without anything else), followed by the
columns `pt,q,M,nGroups,meanNew,errNew`: the pooled mean of `F_q` and its
error over the `nGroups` pooled groups, absolute (not on `ln(F_q)`). The
`.dat` carries the same pooled values as `lnF` + relative error
(`err = sigma(Fq)/mean(Fq)`, i.e. the error on `ln(Fq)`), with `M2` from the
lattice and points `Fq <= 0` written as `-999 -999`, so a pooled `.dat` can be
diffed against any reference in the same layout.

### Subsamples after `hadd` (`slotpool.C`)

Besides the old `sampleSize` path the task therefore fills a second set of
histograms (`slots/...`, see the output section) in which every event is binned
into a hash bucket. `slotpool.C` turns that into subsamples of any `K` by
post-processing alone: it pools consecutive **non-empty**
slots greedily in index order until the accumulated event count reaches `K`, and
recomputes for every pT bin, order and M bin

```
Fq = (Sfq / n) / (Socc / n)^Q        n = sum of slots/nEv over the pooled set
```

where `Q` is the object suffix `o{q}` (2..7), i.e. exactly the exponent
`flushSample()` uses. The mean and its error are taken over the pooled groups;
on the old side the statistics are rebuilt from the merged
`mFqSum`/`mFqSq`/`mSubsampleEvents`, because `mFqError` is per-job and does not
survive `hadd`.

Properties worth knowing:

* **Any `K` works, without rerunning anything.** Splitting is never possible
  after the fact, only merging, so `nSlots` should be chosen generous enough
  that the smallest `K` you will ever want is still reachable; then larger `K`
  are obtained by pooling further.
* `nSlots` only sets the *granularity*. It does not have to relate to `K` or to
  the number of jobs - the pooler never needs to know it.
* **`nSlots` must be identical in every job** of a campaign, otherwise two jobs
  hash the same `globalBC` into different buckets. This is the one failure mode
  that nothing else can see: `hadd` keeps the axis of its first input without
  warning (measured: a 1250-bin and a 1000-bin `slots/nEv` merge with `rc=0` and
  no diagnostic), and since every job still counts each event exactly once, the
  `mEventSelected` invariant passes on the corrupted file too. The task
  therefore writes `slots/params` (bin 1 = number of merged jobs, bin 2 = sum of
  their `nSlots`), and `slotpool.C` prints `jobs merged` / `mean nSlots` and
  raises an `[!! nSlots MISMATCH]` banner whenever that mean disagrees with
  `GetNbinsX()` of `slots/nEv`. `nSlots` is additionally recorded in
  `metaConfig` for humans, but note that titles come from the first `hadd`
  input only, so `metaConfig` alone cannot detect this.
* **Tail**: `K` may not divide the total exactly. With `keepTail=0` (default)
  the trailing partial group is dropped, so every subsample really holds
  `>= K` events; with `keepTail=1` it is kept and may hold `< K`, mirroring
  `flushTailSample`. The macro prints how many events each choice drops.
* **Invariant it asserts**: `mEventSelected[accepted]` must equal the total
  `slots/nEv`, printed as `[OK]`. Any mismatch means events were lost or
  double-counted between the two representations.
* Empty buckets are skipped (they carry `nEv = 0`), so the average occupancy
  `lambda = N_events / nSlots` does not have to be near 1.
* Only `processMCRec` has been exercised end to end; `processRun3` and
  `processMcRun2`/`processRun2` compile and take the same code path but have not
  been run.

**Validation.**

1. *Exactness.* With `sampleSize` set high enough that the old path emits one
   subsample per job, old and new describe the same event set, and with `K = 0`
   (pool everything into a single group) the macro becomes a direct
   comparison of the two computations. On the four-AO2D run (88 accepted
   events) all 936 (pT, q, M) cells agreed to a worst relative difference of
   `1.2e-7` - pure float32-vs-double storage noise.
2. *Post-`hadd` on a real campaign.* 49 AO2Ds run as 49 separate jobs, then
   `hadd`ed, giving 1247 accepted events over 1250 slots (794 non-empty;
   Poisson expects `1250*(1-e^-1) = 790`). Every one of the 49 outputs
   satisfied `sum(slots/nEv) == mEventSelected[accepted]` individually, and
   `slots/params` reports `jobs = 49`, `mean nSlots = 1250`, matching the
   merged axis. At `K = 25` the pooler recovers **49 groups of ~25 events with
   0 events dropped**, i.e. exactly the same subsample count and size as the
   old path, so the two are directly comparable: of the 863 cells where at
   least one side is non-zero, **100 % agree within 1 sigma**, chi2/dof = 0.05,
   rms z = 0.22 and median z = -0.05 at every order q = 2..7 (no bias). The
   chi2 is low because the two paths partition the *same* events and are
   therefore strongly correlated, not because the test is weak. In the
   well-measured core (lowest multiplicity bin, q <= 4) the worst relative
   difference is 0.042 % and the median 0.014 %. At `K = 20` the pooler forms
   60 groups instead (1242 events, 5 dropped to the tail) and the same picture
   holds.

At a finite `K` the two sides are *different* groupings of the same events, so
raw agreement degrades with order and multiplicity even though the estimates
stay statistically consistent: the q=2, M=4 cells match to ~1e-4 while rare
multiplicity bins at high `q` scatter, because a group that happens to contain
none of the rare events records a structural zero. That scatter is the
statistical resolution of the method, not a discrepancy - it is exactly the
error bar the method reports.

### Plotting (`plotFqvsM.C`)

`plotFqvsM.C(<csv>, <ptSel>, <logMode>, <ymin>, <ymax>, <outDir>)` draws the
pooled result straight from the slotpool CSV (it never reads a `.dat`):

* **One PDF per pT bin**, `FactorialMoments_bin<N>.pdf`, with every order
  found in the CSV (`q = 2..7`) drawn as its own series: the per-order
  colours/markers and translucent error bands of the established ATLAS style,
  the cuts in the header (`|eta|`, `|v_z|`, centrality, pT range - all parsed
  from the CSV's `#` header block) and a legend of orders above the curves.
* `x = M^2` is rebuilt from the CSV's `M` index through the lattice
  `bin1/Reset/mEtaPhi<i>` of the root file named in the CSV's `# source` line;
  if that file is not reachable the lattice falls back to `2*(i+2)` - the M
  binning of this analysis, see above - with a warning.
* `ptSel` = `"all"` (default) or a comma list like `"1,3"`. `logMode` = `0`
  automatic (log y as soon as q = 6 or 7 are present, i.e. the combined range
  spans decades), `1` force linear, `2` force log. `ymin`/`ymax` <= 0 choose
  the range automatically (established ±18 % margins); a positive value pins
  that edge. `outDir` defaults to `<csv dir>/plots_FqvsM`.

### Choosing `K` (`findK.C`)

`findK.C(<file>, <kMin>, <kMax>, <csv>)` (defaults `kMin = 5`, `kMax = 200`;
`csv = ""` means no CSV dump) scans candidate `K` values with the exact
`poolSlots()` walk `slotpool.C` uses, so every number it prints is what a real
run at that `K` would produce. Per `K` it reports the groups formed, the
events left in the tail and the worst group-size overshoot (the
indivisible-slot effect), ranks candidates by `score = overshoot% + tail%` and
lists the Pareto options - the `K` for which no other candidate has both less
overshoot and less tail. That table is how `K` was chosen for the campaign
(e.g. `K = 17`: 438 groups, mean size 17.2, 0 events in the tail).

**Provenance.** The full validation record of this workflow - exactness tests,
the post-`hadd` campaign comparison, the `K` scan and all figures - is issue
[#1](https://github.com/maroozm/o2-factorial-moments/issues/1). The validation
scripts of that campaign (`checkFMo2.C`, `validateFMo2.C`, `plotdat.C`,
`plotratio.C`, `scancampaign.sh`, `makeplots_scan.sh`, `verify_ab.sh`, ...) are
archived there as `validation_scripts.zip` and in
`../removed_validation_backup_20261010.tar.gz`; the campaign output data is in
`../validation_data_removed_20261010.tar.gz`. The repository state at the end
of the campaign is preserved as the branch `11/10/2026`.

## 6. Output histograms and merging

Written by device `factorial-moments-task`, directory `factorial-moments-task`,
`q` = 2..7, `N` = pT bin (`bin1..bin3`), 52 M bins (M = 4, 6, ..., 106):

* `mFinalFq{q}_bin{N}`, `mFinalFq{q}Sampled_bin{N}`, `mFinalAvBin*` - filled with
  `Fill(iM, value)` and never reset, so `hadd` adds them. Value =
  `content / (GetEntries()/52)`, where `GetEntries()/52` = number of samples
  (`sampleSize` events each, plus one tail flush per AO2D - e.g. 7+8+8+7 = 30
  for four AO2Ds).
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
* `mSubsampleEvents` - one entry per flushed subsample, x = number of events in
  the group (bin `sampleSize` = complete groups, smaller bins = end-of-stream
  tail flushes).

The merge-safe subsampling representation (see `slotpool.C` above), all under
`slots/` and all written with `SetBinContent`:

* `slots/nEv` (`TH1F`) - `X = slot`, number of accepted events hashed into each
  bucket. This is the `n` of the formulas below and the source of the invariant
  `sum == mEventSelected[accepted]`.
* `slots/params` (`TH1D`, 2 bins) - campaign marker: bin 1 = 1 per job (so it
  accumulates the number of merged jobs), bin 2 = that job's `nSlots` (so it
  accumulates their sum). Lets `slotpool.C` detect jobs that disagreed on
  `nSlots`, which `hadd` and the event-count invariant both miss.
* `slots/sumOcc_pt{N}` (`TH2F`) - `X = slot`, `Y = M`, sum of `binConEvent`
  over the events of that slot and M bin.
* `slots/sumFq_pt{N}_o{q}` (`TH2F`) - `X = slot`, `Y = M`, sum of `fqEvent`,
  `q = 2..7`.

They are plain histograms, so `hadd` adds them cell by cell and the merged file
is indistinguishable from one written by a single job covering the whole sample.
This is the property the TTree alternative could not offer, and the reason a
dense `TH2F` was chosen over `THnSparse`/`StepTHn`: measured on this chain the
slot output is ~9.4 MB in memory and ~7.9 MB on disk, against ~14.0/13.6 MB for
a `StepTHnF` carrying only the sums and ~200 MB for a `THnSparse` of the same
shape. Sum-of-squares is deliberately *not* stored: the spread that matters is
**between** subsamples, and that is produced by the pooler, which knows which
slots it grouped.

## 7. metaConfig

`init()` registers a 1-bin histogram `metaConfig` whose **title** carries the
effective configuration of the job (no `;` in it, otherwise ROOT would split the
string into axis titles):

```
FMtask centLimits=0,5 numPt=3 centralEta=0.9 sampleSize=3 ptMin=0.2 dcaXY=2.4 dcaZ=2 useMC=1 nfqOrder=6 smearPhi=1 flushTail=1 nSlots=1250 vertexXYZ=0.3,0.4,10 ptCuts=0.2,2,0.4,2,0.4,1
```

Histogram titles are not summed by `hadd`, so the merged file keeps this string -
which is how `nSlots` is checked to be consistent across a campaign.
Files produced before `metaConfig` existed fall back to the ranges observed in
the QA histograms - for a truncated distribution that reports the *observed*
range (e.g. `1-5 %` instead of the configured `0-5 %`), so re-running the chain
is recommended once after upgrading the task.

## 8. Implementation notes (kept out of the source, per house style)

* The `M` binning is `binningM[iM] = 2 * (iM + 2)`, i.e. M = 4, 6, ..., 106
  (52 bins). Every sample fills every one of the 52 bins exactly once, so
  `nSamples = GetEntries() / 52`. `slotpool.C` (`.dat` writer) and
  `plotFqvsM.C` read the same M values back from the
  `bin1/Reset/mEtaPhi<i>` binning instead of hardcoding them; when
  `plotFqvsM.C` cannot reach the source root it falls back to `2*(i+2)`,
  which is exactly this binning.
* The factorial numerator uses the falling factorial
  `n!/(n-q)! = n*(n-1)*...*(n-q+1)` evaluated with an explicit product loop:
  `TMath::Factorial` overflows to `inf` for n >= 171.
* `checkpT()` can replace the track azimuth with `gRandom->Gaus(phi, 2*pi)`
  wrapped into `[0, 2*pi)` before filling the eta-phi lattices, i.e. the phi used
  for the factorial moments is effectively randomised. The switch is the
  `smearPhi` option (default `false`, recorded in `metaConfig`); the shipped
  JSON files set it to `true` and the validated `.dat` files were produced with
  it on.
* With `flushTailSample` (default `true`, recorded in `metaConfig` as
  `flushTail=`), any partial group left at the end of a file
  (`countSamples < sampleSize`) is flushed as one extra subsample of
  `countSamples` events, instead of being silently dropped. The flush runs
  inside the framework's end-of-stream hook: a `Service<FlushTailService>`
  element declared **before** `histos` in the task gets its
  `endOfStream()` called before the histogram registry is snapshotted, and it
  calls back into the task (`mSubsampleEvents` records the group size of every
  flush).
  Subsamples are pooled **unweighted** (each sample counts once, as always), so
  a 1-event tail carries the same weight as a complete `sampleSize`-event group
  and its Fq is only resolvable down to one-event statistics - structural zeros
  at high q / large M shift the merged mean by at most ~3 % (measured worst
  case) and inflate the error modestly. Set `flushTailSample: false` to restore
  the old drop-the-tail behaviour; a killed job never reaches the callback, so
  its tail is lost regardless.
* **Why `globalBC` and not a per-file or per-DataFrame index.** The bucket has
  to be a pure function of the event, because Hyperloop is free to chunk the
  input differently from one job to the next and nothing survives the merge
  except the event itself. `collision.globalBC` is the only such identifier with
  the right properties: verified on the production AO2Ds, 868/868 collisions
  carry a unique, non-zero `globalBC` spanning ~1e11, with exactly one distinct
  BC per collision. The raw `fIndexBCs` column is only meaningful *within* one
  DataFrame (38 DFs, 3-38 collisions each) and would silently collapse under
  re-chunking. The mapping is `slot = splitmix64(globalBC) % nSlots`;
  splitmix64 (Steele, Lea, Flood 2014) is a 64-bit integer finaliser, and
  `globalBC` increments linearly, so a plain `% nSlots` would put consecutive
  BCs in consecutive slots and correlate the buckets with time. Measured on
  `AO2D_001`: 638/1250 slots used (626 expected), max occupancy 4,
  chi2/dof = 0.93.
* **Why a dense `TH2F` filled with `SetBinContent`.** The three candidates were
  measured rather than assumed. Filling a `TH2F` with `Fill(x, y, w)` makes
  ROOT allocate `fSumw2` as a `TArrayD` at 8 bytes/cell, inflating the slot
  output to 27.4 MB in memory and 24.6 MB on disk; `SetBinContent` never
  triggers that and needs no sumw2 at all, giving 9.4 MB / 7.9 MB. A
  `StepTHnF` costs ~14.0 / 13.6 MB for the sums alone (plus occupancy
  bookkeeping, ~16.6 / 15.9 MB), emits one INFO line per bucket, and silently
  discards under/overflow. `THnSparse` was rejected outright: the grid is dense
  in practice, so its ~200 MB merge cost is pure overhead. `hadd` of the `TH2F`
  layout was verified on 36 objects including the `nEv` sum.
* The pooler walks slots in **index order**, and that index is hash order, so
  group *g*'s membership depends only on where the boundaries fall and not on
  how the events happen to be distributed over time - which is what lets one
  merged file be reused at any `K`. Pooling can only ever merge, never split,
  so the rule of thumb for choosing `nSlots` is that the smallest `K` you will
  want must stay reachable. Each group stops as soon as it clears `K` and can
  overshoot it by less than the largest bucket it contains, so with `N` events
  the number of groups `G` satisfies `floor(N / (K + b_max)) <= G <= floor(N/K)`,
  where `b_max` is the largest bucket. The two ends coincide while buckets stay
  well below `K` - which is what a generous `nSlots` buys - and only a too-coarse
  `nSlots` costs you groups. Concretely, `[15, 15, 15]` at `K = 20` yields one
  group of 30 plus a 15-event tail, not `floor(45/20) = 2`.
* The recorded `nSlots` in `metaConfig` is useful documentation but is *not* a
  compatibility check: `hadd` takes histogram titles from the first input only,
  so it describes just one job. The real witness is `slots/params`, whose two
  bins accumulate the job count and the sum of the per-job `nSlots`; `hadd` sums
  them, so the mean `bin2/bin1` must equal `GetNbinsX()` of `slots/nEv` or the
  jobs disagreed. This was added after measuring that `hadd` merges histograms
  of differing binning with `rc=0` and **no warning whatsoever** (a 1250-bin and
  a 1000-bin `TH1F` came back as 1250-bin with bin contents summed), and that a
  mismatched campaign still satisfies `sum(slots/nEv) == mEventSelected[accepted]`
  - the totals survive even though the bucket assignment does not.
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
