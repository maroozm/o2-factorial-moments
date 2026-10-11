// Post-hadd reconstruction of Fq from the merge-safe slot output.
//
//   root -l -b -q 'slotpool.C+("results_v2/merged/AnalysisResults.root", 20)'
//
// Arguments: file, K (target events per subsample; <= 0 pools everything into
// one group), csv path ("" = no CSV; otherwise the pooled columns
// pt,q,M,nGroups,meanNew,errNew behind a '#' header block), keepTail (1 =
// retain the trailing partial group, 0 = drop it so every subsample really has
// >= K events), datDir (default "dat_slots": write FactorialMoments_bin*.dat
// there in the checkFMo2.C layout -- new-path values only, diffable against a
// reference .dat; "" = skip).
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
// (Q = 2..7, i.e. the exponent is iq+2 exactly as in the task):
//
//   Fq = (Sfq / n) / (Socc / n)^Q,   n = sum of nEv over S
//
// K <= 0 pools every non-empty slot into ONE group, which reproduces the whole
// event set; against a single-subsample reference that is an exactness test.

#include "TDatime.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TKey.h"
#include "TString.h"
#include "TSystem.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
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
    if (d && d->Get("slots/nEv")) return d;
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

struct Cell {
  double mean = 0.;
  double err = 0.;
};
// pt bin -> M index (0-based, same indexing as the slot TH2 Y axis) -> order
using PtCells = std::map<int, std::map<int, std::map<int, Cell>>>;

const Cell* findCell(const PtCells& c, int p, int m, int q)
{
  auto ip = c.find(p);
  if (ip == c.end()) return nullptr;
  auto im = ip->second.find(m);
  if (im == ip->second.end()) return nullptr;
  auto iq = im->second.find(q);
  if (iq == im->second.end()) return nullptr;
  return &iq->second;
}

// lattice size per M index: cell count of bin1/Reset/mEtaPhi{i}, fallback 2*(i+2)
std::vector<int> latticeM(TDirectory* td, int nM)
{
  std::vector<int> mv;
  mv.reserve(nM);
  for (int i = 0; i < nM; ++i) {
    auto* h = dynamic_cast<TH1*>(td->Get(Form("bin1/Reset/mEtaPhi%d", i)));
    mv.push_back(h ? h->GetNbinsX() : 2 * (i + 2));
  }
  return mv;
}

// .dat byte layout follows checkFMo2.C's processBin() (archived with the
// validation campaign), with the slot-pooled statistics

// read one space-delimited value out of the metaConfig histogram title
static std::string metaValue(TDirectory* td, const char* key)
{
  auto* meta = dynamic_cast<TH1*>(td->Get("metaConfig"));
  if (!meta)
    return {};
  const std::string cfg = meta->GetTitle();
  const std::string k = std::string(key) + "=";
  auto p = cfg.find(k);
  if (p == std::string::npos)
    return {};
  p += k.size();
  auto e = cfg.find(' ', p);
  return (e == std::string::npos) ? cfg.substr(p) : cfg.substr(p, e - p);
}

static std::vector<double> listNums(const std::string& s)
{
  std::vector<double> v;
  std::string tok;
  for (char c : s) {
    if (c == ',' || c == ' ' || c == '\t' || c == '\n') {
      if (!tok.empty()) {
        v.push_back(std::atof(tok.c_str()));
        tok.clear();
      }
    } else {
      tok.push_back(c);
    }
  }
  if (!tok.empty())
    v.push_back(std::atof(tok.c_str()));
  return v;
}

// Centrality range for the headers: metaConfig's centLimits when present,
// otherwise the observed span of the first available centrality histogram
// (same candidate order checkFMo2.C uses).
static bool centRange(TDirectory* td, double* lo, double* hi)
{
  const auto cl = listNums(metaValue(td, "centLimits"));
  if (cl.size() >= 2) {
    *lo = cl[0];
    *hi = cl[1];
    return true;
  }
  const char* cand[4] = {"mCentFT0C", "mCentFT0M", "mCentFV0A", "mCentFT0A"};
  for (int i = 0; i < 4; ++i) {
    auto* h = dynamic_cast<TH1*>(td->Get(cand[i]));
    if (!h)
      continue;
    int first = -1, last = -1;
    for (int b = 1; b <= h->GetNbinsX(); ++b)
      if (h->GetBinContent(b) > 0.) {
        if (first < 0)
          first = b;
        last = b;
      }
    if (first > 0) {
      *lo = h->GetXaxis()->GetBinLowEdge(first);
      *hi = h->GetXaxis()->GetBinUpEdge(last);
      return true;
    }
  }
  return false;
}

// Build checkFMo2.C's "# pT bin N: <range> | centrality .. | subsamples .."
// header line from the file's metaConfig, so a slot .dat is as self-describing
// as one written straight from an AnalysisResults.root.
static std::string ptHeader(TDirectory* td, int ptBin, int nGroups)
{
  std::string s = Form("# pT bin %d:", ptBin);
  const auto pt = listNums(metaValue(td, "ptCuts"));
  if ((int)pt.size() >= 2 * ptBin)
    s += Form(" %.2f-%.2f GeV/c", pt[2 * (ptBin - 1)], pt[2 * (ptBin - 1) + 1]);
  else
    s += Form(" bin %d", ptBin);

  double centLo = 0., centHi = 0.;
  if (centRange(td, &centLo, &centHi))
    s += Form(" | centrality %.0f-%.0f %%", centLo, centHi);
  s += Form(" | subsamples %d", nGroups);

  const auto eta = metaValue(td, "centralEta");
  if (!eta.empty())
    s += Form(" | |eta| <= %.2f", std::atof(eta.c_str()));
  const auto vtx = listNums(metaValue(td, "vertexXYZ"));
  if (vtx.size() >= 3)
    s += Form(" | |vz| <= %.1f", vtx[2]);
  return s;
}

void writeDat(const std::string& path, const char* src, const PtCells& cells,
              int ptBin, int nBins, const std::vector<int>& mv,
              const std::vector<int>& orders, int nGroups,
              const std::string& header, int K)
{
  FILE* fp = fopen(path.c_str(), "w");
  if (!fp) { printf("  !! cannot write %s\n", path.c_str()); return; }
  fprintf(fp, "# source: %s\n", src);
  fprintf(fp, "%s\n", header.c_str());
  fprintf(fp, "# method: slot-pooled, K=%d, groups=%d\n", K, nGroups);
  fprintf(fp, "#M2 lnM2");
  for (int q : orders) fprintf(fp, " lnF%d errF%d", q, q);
  fprintf(fp, "\n");
  fprintf(fp, "# err = sigma(Fq)/mean(Fq), i.e. the error on ln(Fq);"
              " entries with Fq <= 0 are written as -999 -999\n");
  for (int m = 0; m < nBins; ++m) {
    // checkFMo2.C writes mval*mval as an int (setprecision has no effect on
    // integers), so print it as one or the first column would not match.
    const int m2 = mv[m] * mv[m];
    fprintf(fp, "%d %.8f", m2, std::log(static_cast<double>(m2)));
    for (int q : orders) {
      const Cell* c = findCell(cells, ptBin, m, q);
      if (c && c->mean > 0. && std::isfinite(c->mean))
        fprintf(fp, " %.8f %.8f", std::log(c->mean), c->err / c->mean);
      else
        fprintf(fp, " -999 -999");
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  printf("  wrote %s\n", path.c_str());
}
} // namespace

void slotpool(const char* fileName,
              int K = 20,
              const char* csvFile = "",
              int keepTail = 0,
              const char* datDir = "dat_slots")
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
  int accepted = -1;
  auto* sel = dynamic_cast<TH1*>(td->Get("mEventSelected"));
  if (sel) {
    accepted = static_cast<int>(sel->GetBinContent(10)); // kEvAccepted = 10
    printf("  mEventSelected accepted: %d   -> slots total: %d   %s\n",
           accepted, pooled + dropped,
           accepted == pooled + dropped ? "[OK]" : "[!! MISMATCH]");
  }
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

  FILE* csv = nullptr;
  if (csvFile && csvFile[0]) {
    csv = fopen(csvFile, "w");
    if (csv) {
      // '#' header block: self-describing provenance, ignored by '#'-skipping readers
      TString full = gSystem->ExpandPathName(fileName);
      if (!full.BeginsWith("/")) full = TString(gSystem->pwd()) + "/" + full;
      fprintf(csv, "# source : %s\n", full.Data());
      fprintf(csv, "# date   : %s\n", TDatime().AsString());
      fprintf(csv, "# K      : %d      keepTail: %d      subsamples: %zu\n", K, keepTail,
              groups.size());
      fprintf(csv, "# nSlots : %d (%d non-empty)      events: %d", G, nonEmpty, pooled + dropped);
      if (accepted >= 0) fprintf(csv, " (accepted: %d)", accepted);
      fprintf(csv, "\n");
      double centLo = 0., centHi = 0.;
      if (centRange(td, &centLo, &centHi))
        fprintf(csv, "# cent   : %.0f-%.0f %%\n", centLo, centHi);
      const auto eta = metaValue(td, "centralEta");
      if (!eta.empty())
        fprintf(csv, "# eta    : |eta| <= %.2f\n", std::atof(eta.c_str()));
      const auto vtx = listNums(metaValue(td, "vertexXYZ"));
      if (vtx.size() >= 3)
        fprintf(csv, "# vertex : |vz| <= %.1f cm\n", vtx[2]);
      const auto pt = listNums(metaValue(td, "ptCuts"));
      std::string ptLine = "# pt bins: ";
      for (size_t i = 0; i < ptBins.size(); ++i) {
        const int p = ptBins[i];
        if (i) ptLine += "  ";
        if ((int)pt.size() >= 2 * p)
          ptLine += Form("%d: %.2f-%.2f GeV/c", p, pt[2 * (p - 1)], pt[2 * (p - 1) + 1]);
        else
          ptLine += Form("%d: bin %d", p, p);
      }
      if (ptLine.find(':') != std::string::npos) fprintf(csv, "%s\n", ptLine.c_str());
      fprintf(csv, "pt,q,M,nGroups,meanNew,errNew\n");
    }
  }

  printf("=== Fq (slot-pooled) ===\n");
  printf("%3s %3s %3s | %10s %10s %6s\n",
         "pt", "q", "M", "meanNew", "errNew", "nGroups");

  PtCells cells; // pooled statistics, buffered so the .dat can be written in M order

  for (int p : ptBins) {
    auto* occ = get2(td, Form("slots/sumOcc_pt%d", p));
    if (!occ) continue;
    for (int q : orders) {
      auto* fq = get2(td, Form("slots/sumFq_pt%d_o%d", p, q));
      if (!fq) continue;

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

        { Cell c; c.mean = sn.mean; c.err = sn.err; cells[p][m][q] = c; }

        if (csv) {
          fprintf(csv, "%d,%d,%d,%d,%.10g,%.10g\n",
                  p, q, m + 1, sn.n, sn.mean, sn.err);
        }
        // only print a few M per (p,q) to keep the table readable
        if (m == 0 || m == fq->GetNbinsY() / 2 || m == fq->GetNbinsY() - 1) {
          printf("%3d %3d %3d | %10.5f %10.5f %6d\n", p, q, m + 1, sn.mean, sn.err, sn.n);
        }
      }
    }
  }

  if (csv) { fclose(csv); printf("\nwrote %s\n", csvFile); }

  if (datDir && datDir[0]) {
    gSystem->mkdir(datDir, kTRUE);
    printf("\n=== .dat (checkFMo2 layout, from pooled slots) ===\n");
    for (int p : ptBins) {
      int nM = 0;
      if (auto* occY = dynamic_cast<TH2*>(td->Get(Form("slots/sumOcc_pt%d", p))))
        nM = occY->GetNbinsY();
      if (nM <= 0) { printf("  !! pt bin %d: cannot determine M bins\n", p); continue; }
      const std::vector<int> mv = latticeM(td, nM);
      const std::string path =
        std::string(datDir) + Form("/FactorialMoments_bin%d.dat", p);
      const int nGrp = static_cast<int>(groups.size());
      writeDat(path, fileName, cells, p, nM, mv, orders, nGrp,
               ptHeader(td, p, nGrp), K);
    }
  }

  f->Close();
}
