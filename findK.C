// findK.C -- pick good subsample sizes K for the slot-pooled path.
//
//   root -l -b -q 'findK.C+("results_scan/S17/merged/AnalysisResults.root")'
//   root -l -b -q 'findK.C+("file.root", 5, 200)'                // custom K range
//   root -l -b -q 'findK.C+("file.root", 5, 200, "findK.csv")'   // + CSV dump
//
// Scans candidate K values and reports, for each, the two things that go wrong
// when pooling whole slots into subsamples of >= K events:
//   overshoot -- the last slot added to a group pushes the size past K, so the
//                subsample is larger than requested (slots are indivisible);
//   tail      -- leftover events at the end that never reach K (dropped when
//                keepTail=0, kept as one undersized group when keepTail=1).
// The walk is slotpool.C's poolSlots() walk verbatim (bins in index order,
// empty slots skipped, group closed when the accumulated nEv reaches K), so
// every number here is what a slotpool run with that K would produce.
//
// Ranking: score = OS% + tail% (both are % of events, so comparable).
// Also lists the Pareto options -- K's for which no other candidate has BOTH
// less overshoot and less tail; those are the honest trade-off choices.
//
// Reads, like slotpool.C: slots/nEv (+ slots/params, mEventSelected, metaConfig
// for the header/consistency lines).

#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "TKey.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace
{
// same directory probe as slotpool.C
TDirectory* taskDir(TFile* f)
{
  if (!f) return nullptr;
  if (auto* named = dynamic_cast<TDirectory*>(f->Get("factorial-moments-task"))) return named;
  TIter next(f->GetListOfKeys());
  while (auto* k = (TKey*)next()) {
    if (std::string(k->GetClassName()).find("TDirectory") == std::string::npos) continue;
    auto* d = dynamic_cast<TDirectory*>(f->Get(k->GetName()));
    if (d && (d->Get("slots/nEv") || d->Get("mSubsampleEvents"))) return d;
  }
  return nullptr;
}

// same walk as slotpool.C's poolSlots() with keepTail=false; returns the group
// sizes (each >= K). *tail = leftover events (< K), *pooled = sum of the groups.
std::vector<int> poolSizes(TH1* nEv, int K, int* pooled, int* tail)
{
  std::vector<int> sizes;
  double acc = 0;
  int cur = 0;
  const int G = nEv->GetNbinsX();
  for (int bin = 1; bin <= G; ++bin) {
    const double ne = nEv->GetBinContent(bin);
    if (ne <= 0) continue;
    cur += static_cast<int>(ne);
    acc += ne;
    if (acc >= K) {
      sizes.push_back(cur);
      *pooled += cur;
      cur = 0;
      acc = 0;
    }
  }
  *tail = static_cast<int>(acc);
  return sizes;
}

struct Row {
  int K = 0;
  int groups = 0;
  int pooled = 0;
  int tail = 0;
  int minSz = 0;
  int maxSz = 0;
  double avg = 0;      // mean group size
  double osPct = 0;    // mean relative overshoot, % of K
  double maxOv = 0;    // worst single overshoot, events
  double exactPct = 0; // % of groups landing exactly on K
  double tailPct = 0;  // tail as % of all events
  double score = 0;    // osPct + tailPct
};
} // namespace

void findK(const char* fileName, int kmin = 3, int kmax = 100, const char* csvFile = "")
{
  auto* f = TFile::Open(fileName, "READ");
  if (!f || f->IsZombie()) { printf("cannot open %s\n", fileName); return; }
  TDirectory* td = taskDir(f);
  if (!td) { printf("no task output directory in %s\n", fileName); return; }
  auto* nEv = dynamic_cast<TH1*>(td->Get("slots/nEv"));
  if (!nEv) { printf("no slots/nEv in %s -- not a slot-enabled output\n", fileName); return; }

  const int G = nEv->GetNbinsX();
  int occupied = 0;
  long events = 0;
  for (int b = 1; b <= G; ++b) {
    const double ne = nEv->GetBinContent(b);
    if (ne > 0) { ++occupied; events += static_cast<long>(ne); }
  }

  printf("=== findK: candidate K scan ===\n");
  printf("  file    : %s\n", fileName);
  printf("  slots   : nSlots=%d  occupied=%d  events=%ld\n", G, occupied, events);
  if (auto* par = dynamic_cast<TH1*>(td->Get("slots/params"))) {
    const double nJobs = par->GetBinContent(1);
    const double sumN = par->GetBinContent(2);
    printf("  jobs merged : %g   mean nSlots : %g%s\n", nJobs, nJobs > 0 ? sumN / nJobs : 0.,
           (nJobs > 0 && std::fabs(sumN / nJobs - static_cast<double>(G)) > 1e-6)
               ? "   [!! nSlots MISMATCH -- buckets meaningless]"
               : "");
  }
  if (auto* sel = dynamic_cast<TH1*>(td->Get("mEventSelected"))) {
    const int accepted = static_cast<int>(sel->GetBinContent(10)); // kEvAccepted = 10
    printf("  accepted events    : %d   -> slots total: %ld   %s\n", accepted, events,
           accepted == events ? "[OK]" : "[!! MISMATCH]");
  }
  if (events == 0) { printf("no events in slots -- nothing to scan\n"); return; }
  if (kmin < 1) kmin = 1;
  if (kmax < kmin) kmax = kmin;
  printf("  scan    : K = %d..%d\n", kmin, kmax);
  printf("  note    : tail = leftover events < K (dropped when keepTail=0);\n"
         "            score = OS%% + tail%%, both in %% of events\n\n");

  std::vector<Row> rows;
  int skipped = 0;
  for (int K = kmin; K <= kmax; ++K) {
    int pooled = 0, tail = 0;
    const auto sizes = poolSizes(nEv, K, &pooled, &tail);
    if (sizes.empty()) { ++skipped; continue; }
    Row r;
    r.K = K;
    r.groups = static_cast<int>(sizes.size());
    r.pooled = pooled;
    r.tail = tail;
    r.minSz = *std::min_element(sizes.begin(), sizes.end());
    r.maxSz = *std::max_element(sizes.begin(), sizes.end());
    r.avg = static_cast<double>(pooled) / r.groups;
    const double sumOv = static_cast<double>(pooled) - static_cast<double>(r.groups) * K;
    r.osPct = 100. * sumOv / (static_cast<double>(r.groups) * K);
    r.maxOv = static_cast<double>(r.maxSz) - K;
    int exact = 0;
    for (int s : sizes)
      if (s == K) ++exact;
    r.exactPct = 100. * exact / r.groups;
    r.tailPct = 100. * tail / static_cast<double>(events);
    r.score = r.osPct + r.tailPct;
    if (pooled + tail != static_cast<int>(events))
      printf("  [!!] K=%d: pooled %d + tail %d != events %ld\n", K, pooled, tail, events);
    rows.push_back(r);
  }
  if (rows.empty()) { printf("no K in %d..%d forms a single subsample -- lower K\n", kmin, kmax); return; }

  printf("%5s %7s %6s %9s %7s %6s %7s %6s %7s %7s  %s\n", "K", "groups", "avg", "min..max",
         "OS%", "maxOS", "exact%", "tail", "tail%", "score", "note");
  for (const auto& r : rows) {
    const char* note = r.groups < 30 ? "few groups" : "";
    printf("%5d %7d %6.1f %4d..%-4d %6.2f %6.0f %6.0f%% %6d %6.2f%% %7.2f  %s\n", r.K, r.groups, r.avg,
           r.minSz, r.maxSz, r.osPct, r.maxOv, r.exactPct, r.tail, r.tailPct, r.score, note);
  }
  if (skipped) printf("  (%d K in range formed no subsample and were skipped)\n", skipped);

  // best by combined score
  const Row* best = &rows.front();
  for (const auto& r : rows)
    if (r.score < best->score) best = &r;

  // Pareto front over (OS%, tail%): no other K is better in both
  printf("\n--- Pareto options (no other K has BOTH less overshoot and less tail) ---\n");
  for (const auto& r : rows) {
    bool dominated = false;
    for (const auto& o : rows) {
      if (o.K == r.K) continue;
      if (o.osPct <= r.osPct && o.tailPct <= r.tailPct &&
          (o.osPct < r.osPct || o.tailPct < r.tailPct)) {
        dominated = true;
        break;
      }
    }
    if (!dominated)
      printf("  K=%-4d %6d groups, size %.1f avg (%d..%d), overshoot +%.2f%% (max +%g), "
             "tail %d (%.2f%%), score %.2f%s\n",
             r.K, r.groups, r.avg, r.minSz, r.maxSz, r.osPct, r.maxOv, r.tail, r.tailPct, r.score,
             r.groups < 30 ? "   [few groups]" : "");
  }

  // extremes, spelled out
  double minOs = 1e30;
  int minOsK = 0, minTailK = 0;
  double minTail = 1e30;
  for (const auto& r : rows) {
    if (r.osPct < minOs) { minOs = r.osPct; minOsK = r.K; }
    if (r.tailPct < minTail - 1e-12) { minTail = r.tailPct; minTailK = r.K; }
  }
  printf("\n--- extremes ---\n");
  printf("  least overshoot : K=%d  (OS%% %.2f)\n", minOsK, minOs);
  printf("  least tail      : K=%d  (tail%% %.2f)\n", minTailK, minTail);
  printf("  best score      : K=%d  (score %.2f = OS%% %.2f + tail%% %.2f)\n", best->K, best->score,
         best->osPct, best->tailPct);

  if (csvFile && *csvFile) {
    FILE* fp = fopen(csvFile, "w");
    if (fp) {
      fprintf(fp, "K,groups,pooled,tail,avg,minSize,maxSize,osPct,maxOvershoot,exactPct,tailPct,score\n");
      for (const auto& r : rows)
        fprintf(fp, "%d,%d,%d,%d,%.4f,%d,%d,%.4f,%.0f,%.4f,%.4f,%.4f\n", r.K, r.groups, r.pooled, r.tail,
                r.avg, r.minSz, r.maxSz, r.osPct, r.maxOv, r.exactPct, r.tailPct, r.score);
      fclose(fp);
      printf("\nCSV written: %s\n", csvFile);
    } else {
      printf("\n!! cannot write %s\n", csvFile);
    }
  }
}
