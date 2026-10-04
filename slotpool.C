// Post-hadd reconstruction of Fq from the merge-safe slot output.
//
//   root -l -b -q 'slotpool.C+("results_v2/merged/AnalysisResults.root", 20, "results_v2/merged/AnalysisResults.root")'
//
// Arguments: file, K (target events per subsample; <= 0 pools everything into
// one group), baseline file ("" to skip the old-path comparison), csv path,
// keepTail (1 = retain the trailing partial group, 0 = drop it so every
// subsample really has >= K events).
//
// Input objects written by FactorialMomentsTask (under the
// factorial-moments-task/ directory O2 creates for the device):
//   slots/nEv              TH1F  events per slot
//   slots/sumOcc_pt{p}     TH2F  X = slot, Y = M, sum of occupancy over events
//   slots/sumFq_pt{p}_o{q} TH2F  X = slot, Y = M, sum of fqEvent over events
//
// Slots are pooled greedily in index order until the accumulated nEv reaches K,
// so every subsample has >= K events and the bucket granularity (nSlots) never
// has to be known here. Empty slots are skipped, a trailing partial group is
// dropped. For a pooled group S, with Q the order suffix of the object
// (Q = 2..7, i.e. the exponent is iq+2 exactly as in flushSample):
//
//   Fq = (Sfq / n) / (Socc / n)^Q,   n = sum of nEv over S
//
// This is the same expression flushSample() uses on the old path. The sums are
// accumulated in a different order (per slot, then pooled), so agreement with
// the old path is expected at roundoff level rather than bit-for-bit.
//
// K <= 0 pools every non-empty slot into ONE group, which reproduces the whole
// event set; against a single-subsample reference that is an exactness test.
//
// With a baseline given, the old mFqSum/mFqSq/mSubsampleEvents statistics are
// recomputed from their merged sums (mFqError is per-job and does not survive
// hadd) and printed alongside.

#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TKey.h"
#include "TString.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace
{
struct Stat {
  int n = 0;
  double mean = 0;
  double err = 0;
};

Stat statOf(const std::vector<double>& v)
{
  Stat s;
  s.n = static_cast<int>(v.size());
  if (s.n == 0) return s;
  double sum = 0, sq = 0;
  for (double x : v) { sum += x; sq += x * x; }
  s.mean = sum / s.n;
  if (s.n > 1) {
    const double var = (sq - s.n * s.mean * s.mean) / (s.n - 1);
    if (var > 0) s.err = std::sqrt(var / s.n);
  }
  return s;
}

// O2 stores HistogramRegistry output in a TDirectoryFile named after the task
// device, not at the top level. Prefer the conventional name, else probe for a
// directory that actually holds the task's histograms.
TDirectory* taskDir(TFile* f)
{
  if (!f) return nullptr;
  auto* named = dynamic_cast<TDirectory*>(f->Get("factorial-moments-task"));
  if (named) return named;
  TIter next(f->GetListOfKeys());
  while (auto* k = (TKey*)next()) {
    if (std::string(k->GetClassName()).find("TDirectory") == std::string::npos) continue;
    auto* d = dynamic_cast<TDirectory*>(f->Get(k->GetName()));
    if (d && (d->Get("slots/nEv") || d->Get("mSubsampleEvents"))) return d;
  }
  return nullptr;
}

// greedy pooling of slots into groups with >= K events.
//   K <= 0        -> one single group holding every non-empty slot
//   keepTail      -> also keep the trailing partial group (which then has < K
//                    events, exactly like flushTailSample on the old path)
// Returns groups; *droppedEvents is the tail discarded when keepTail is false.
std::vector<std::vector<int>> poolSlots(TH1* nEv, int K, bool keepTail, int* pooledEvents, int* droppedEvents)
{
  std::vector<std::vector<int>> groups;
  std::vector<int> cur;
  double acc = 0;
  const int G = nEv->GetNbinsX();
  for (int bin = 1; bin <= G; ++bin) {
    const double ne = nEv->GetBinContent(bin);
    if (ne <= 0) continue; // empty slot contributes nothing
    cur.push_back(bin);
    acc += ne;
    if (K > 0 && acc >= K) {
      groups.push_back(cur);
      cur.clear();
      acc = 0;
    }
  }
  int tail = static_cast<int>(acc);
  if (!cur.empty() && (K <= 0 || keepTail)) {
    groups.push_back(cur);
    tail = 0;
  }
  int pooled = 0;
  for (const auto& g : groups)
    for (int bin : g) pooled += static_cast<int>(nEv->GetBinContent(bin));
  *pooledEvents = pooled;
  *droppedEvents = tail;
  return groups;
}

TH2* get2(TDirectory* d, const TString& name)
{
  auto* h = dynamic_cast<TH2*>(d->Get(name.Data()));
  if (!h) printf("  !! missing %s\n", name.Data());
  return h;
}
} // namespace

void slotpool(const char* fileName,
              int K = 20,
              const char* baselineFile = "",
              const char* csvFile = "slotpool.csv",
              int keepTail = 0)
{
  auto* f = TFile::Open(fileName, "READ");
  if (!f || f->IsZombie()) { printf("cannot open %s\n", fileName); return; }
  TDirectory* td = taskDir(f);
  if (!td) { printf("no task output directory in %s\n", fileName); return; }
  auto* nEv = dynamic_cast<TH1*>(td->Get("slots/nEv"));
  if (!nEv) { printf("no slots/nEv in %s -- not a slot-enabled output\n", fileName); return; }

  int pooled = 0, dropped = 0;
  const auto groups = poolSlots(nEv, K, keepTail != 0, &pooled, &dropped);
  const int G = nEv->GetNbinsX();
  int nonEmpty = 0;
  for (int b = 1; b <= G; ++b)
    if (nEv->GetBinContent(b) > 0) ++nonEmpty;

  printf("=== pooling ===\n");
  printf("  file                : %s\n", fileName);
  printf("  slots               : %d (%d non-empty)\n", G, nonEmpty);
  // nSlots consistency: hadd keeps the first input's axis silently, so the only
  // witness is this accumulated marker (see FactorialMomentsTask::init).
  if (auto* par = dynamic_cast<TH1*>(td->Get("slots/params"))) {
    const double nJobs = par->GetBinContent(1);
    const double sumN = par->GetBinContent(2);
    printf("  jobs merged         : %g   mean nSlots : %g\n", nJobs, nJobs > 0 ? sumN / nJobs : 0.);
    if (nJobs > 0 && std::fabs(sumN / nJobs - static_cast<double>(G)) > 1e-6) {
      printf("\n"
             "  [!! nSlots MISMATCH] jobs were run with DIFFERENT nSlots: the merged\n"
             "  slots/nEv axis has %d bins but the jobs averaged %g. hadd keeps the first\n"
             "  input's axis without warning and the event totals still balance, so this\n"
             "  file's buckets are meaningless -- rerun those jobs with one nSlots.\n\n",
             G, sumN / nJobs);
    }
  } else {
    printf("  jobs merged         : unknown (no slots/params: output predates the campaign marker)\n");
  }
  printf("  target K            : %d%s\n", K, K <= 0 ? "  (K<=0: single group over all slots)" : "");
  printf("  keep trailing group : %s\n", keepTail ? "yes (tail may have < K events)" : "no (tail dropped)");
  printf("  subsamples formed   : %zu\n", groups.size());
  printf("  events pooled       : %d\n", pooled);
  printf("  events dropped (tail): %d\n", dropped);
  auto* sel = dynamic_cast<TH1*>(td->Get("mEventSelected"));
  if (sel) {
    const int accepted = static_cast<int>(sel->GetBinContent(10)); // kEvAccepted = 10
    printf("  mEventSelected accepted: %d   -> slots total: %d   %s\n",
           accepted, pooled + dropped,
           accepted == pooled + dropped ? "[OK]" : "[!! MISMATCH]");
  }
  auto* subs = dynamic_cast<TH1*>(td->Get("mSubsampleEvents"));
  if (subs) printf("  old-path subsamples : %d\n", static_cast<int>(subs->GetEntries()));
  if (groups.empty()) { printf("no complete subsample -- lower K or raise nSlots\n"); return; }

  // discover numPt / nfqOrder from what is actually present
  std::vector<int> ptBins;
  for (int p = 1; p <= 8; ++p)
    if (td->Get(Form("slots/sumOcc_pt%d", p))) ptBins.push_back(p);
  std::vector<int> orders;
  for (int q = 2; q <= 12; ++q)
    if (td->Get(Form("slots/sumFq_pt%d_o%d", ptBins.empty() ? 1 : ptBins.front(), q))) orders.push_back(q);
  printf("  pt bins found       : %zu, orders found: %zu\n\n", ptBins.size(), orders.size());
  if (ptBins.empty() || orders.empty()) return;

  // ---- optional baseline ----
  TFile* fb = nullptr;
  TDirectory* tdb = nullptr;
  int nOld = 0;
  if (baselineFile && baselineFile[0]) {
    fb = TFile::Open(baselineFile, "READ");
    if (!fb || fb->IsZombie()) { printf("cannot open baseline %s\n", baselineFile); fb = nullptr; }
    else {
      tdb = taskDir(fb);
      auto* se = tdb ? dynamic_cast<TH1*>(tdb->Get("mSubsampleEvents")) : nullptr;
      if (se) nOld = static_cast<int>(se->GetEntries());
      printf("=== baseline ===\n  %s\n  old subsamples: %d\n\n", baselineFile, nOld);
    }
  }

  FILE* csv = nullptr;
  if (csvFile && csvFile[0]) {
    csv = fopen(csvFile, "w");
    if (csv) fprintf(csv, "pt,q,M,nGroups,meanNew,errNew,meanOld,errOld\n");
  }

  printf("=== Fq (new slots vs old) ===\n");
  printf("%3s %3s %3s | %10s %10s %-20s | %10s %10s\n",
         "pt", "q", "M", "meanNew", "errNew", "ratio new/old", "meanOld", "errOld");

  double worstRel = 0;
  int worstQ = -1, worstPt = -1, worstM = -1;
  double worstNew = 0, worstOld = 0;

  for (int p : ptBins) {
    auto* occ = get2(td, Form("slots/sumOcc_pt%d", p));
    if (!occ) continue;
    for (int q : orders) {
      auto* fq = get2(td, Form("slots/sumFq_pt%d_o%d", p, q));
      if (!fq) continue;

      auto* oldSum = (tdb && fb) ? dynamic_cast<TH1*>(tdb->Get(Form("mFqSum%d_bin%d", q, p))) : nullptr;
      auto* oldSq = (tdb && fb) ? dynamic_cast<TH1*>(tdb->Get(Form("mFqSq%d_bin%d", q, p))) : nullptr;

      for (int m = 0; m < fq->GetNbinsY(); ++m) { // M index, matches nBins
        std::vector<double> vals;
        vals.reserve(groups.size());
        for (const auto& g : groups) {
          double sfq = 0, socc = 0, n = 0;
          for (int bin : g) {
            sfq += fq->GetBinContent(bin, m + 1);
            socc += occ->GetBinContent(bin, m + 1);
            n += nEv->GetBinContent(bin);
          }
          if (n <= 0) continue;
          const double av = socc / n;
          const double den = std::pow(av, q);
          vals.push_back(den > 0 ? (sfq / n) / den : 0.);
        }
        const Stat sn = statOf(vals);

        Stat so;
        if (oldSum && oldSq && nOld > 1) {
          const double s = oldSum->GetBinContent(m + 1);
          const double s2 = oldSq->GetBinContent(m + 1);
          so.n = nOld;
          so.mean = s / nOld;
          const double var = (s2 - nOld * so.mean * so.mean) / (nOld - 1);
          if (var > 0) so.err = std::sqrt(var / nOld);
        } else if (oldSum && nOld == 1) {
          so.n = 1;
          so.mean = oldSum->GetBinContent(m + 1); // single subsample: no spread
        }

        if (csv) {
          fprintf(csv, "%d,%d,%d,%d,%.10g,%.10g,%.10g,%.10g\n",
                  p, q, m + 1, sn.n, sn.mean, sn.err, so.mean, so.err);
        }
        // only print a few M per (p,q) to keep the table readable
        if (m == 0 || m == fq->GetNbinsY() / 2 || m == fq->GetNbinsY() - 1) {
          printf("%3d %3d %3d | %10.5f %10.5f ", p, q, m + 1, sn.mean, sn.err);
          if (so.n > 0 && std::fabs(so.mean) > 0)
            printf("%-20.4f ", sn.mean / so.mean);
          else
            printf("%-20s ", "-");
          printf("| %10.5f %10.5f\n", so.mean, so.err);
        }
        if (so.n > 0 && std::fabs(so.mean) > 1e-12) {
          const double rel = std::fabs(sn.mean / so.mean - 1.);
          if (rel > worstRel) { worstRel = rel; worstQ = q; worstPt = p; worstM = m + 1; worstNew = sn.mean; worstOld = so.mean; }
        }
      }
    }
  }

  if (csv) { fclose(csv); printf("\nwrote %s\n", csvFile); }

  if (worstQ > 0) {
    printf("\n=== worst relative difference ===\n");
    printf("  pt=%d q=%d M=%d  new=%.6f  old=%.6f  rel=%.4g%%\n",
           worstPt, worstQ, worstM, worstNew, worstOld, 100 * worstRel);
  }
  f->Close();
  if (fb && fb != f) fb->Close();
}
