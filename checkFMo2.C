/// \file checkFMo2.C
/// \brief Draw ln(F_q) vs ln(M^2) for every pT bin of an AnalysisResults.root (single job or hadd-ed) and dump the values to text files
/// \author Salman Malik
/// \author Balwan Singh

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

struct FMSetup {
  int nBins = 0;
  int nPtBins = 0;
  int nOrders = 0;
  int nSubsamples = 0;
  bool haveMoments = false;
  bool haveConfig = false;
  std::vector<int> mValues;
  std::vector<double> ptLow;
  std::vector<double> ptHigh;
  double centLow = 0.;
  double centHigh = 0.;
  double etaMax = 0.;
  double vzMax = 0.;
  TString centSource;
};

static TString joinPath(const char *dir, const TString &name)
{
  TString p(dir);
  if (!p.EndsWith("/"))
    p += "/";
  p += name;
  return p;
}

static bool configValue(const TString &cfg, const char *key, TString &val)
{
  TString k(key);
  k += "=";
  Ssiz_t p = cfg.Index(k);
  if (p == kNPOS)
    return false;
  p += k.Length();
  Ssiz_t e = cfg.Index(' ', p);
  val = (e == kNPOS) ? cfg(p, cfg.Length() - p) : cfg(p, e - p);
  return val.Length() > 0;
}

static bool parseList(const TString &s, std::vector<double> &out)
{
  out.clear();
  TString rest = s;
  while (rest.Length() > 0) {
    Ssiz_t c = rest.Index(',');
    TString tok = (c == kNPOS) ? rest : rest(0, c);
    const char *str = tok.Data();
    char *end = nullptr;
    double v = std::strtod(str, &end);
    if (end == str)
      return false;
    out.push_back(v);
    if (c == kNPOS)
      break;
    rest = rest(c + 1, rest.Length() - c - 1);
  }
  return !out.empty();
}

static bool parseTitleRange(const TString &title, double &lo, double &hi)
{
  Ssiz_t p = title.Index("bin ");
  if (p == kNPOS)
    return false;
  TString rng = title(p + 4, title.Length() - (p + 4));
  Ssiz_t semi = rng.Index(';');
  if (semi != kNPOS)
    rng = rng(0, semi);
  return sscanf(rng.Data(), "%lf-%lf", &lo, &hi) == 2;
}

static bool filledRange(TH1 *h, double &lo, double &hi)
{
  if (!h)
    return false;
  int first = -1, last = -1;
  for (int b = 1; b <= h->GetNbinsX(); b++) {
    if (h->GetBinContent(b) > 0.) {
      if (first < 0)
        first = b;
      last = b;
    }
  }
  if (first < 0)
    return false;
  lo = h->GetXaxis()->GetBinLowEdge(first);
  hi = h->GetXaxis()->GetBinUpEdge(last);
  return true;
}

static FMSetup detectSetup(TDirectory *tdir)
{
  FMSetup st;

  while (st.nPtBins < 32 && tdir->Get(Form("mFinalFq2Sampled_bin%i", st.nPtBins + 1)))
    st.nPtBins++;

  while (st.nOrders < 20 && tdir->Get(Form("mFinalFq%iSampled_bin1", st.nOrders + 2)))
    st.nOrders++;

  TH1 *hRef = (TH1 *)tdir->Get("mFinalFq2Sampled_bin1");
  st.nBins = hRef ? hRef->GetNbinsX() : 0;
  st.haveMoments = (tdir->Get("mFqSum2_bin1") != nullptr);

  for (int i = 0; i < st.nBins; i++) {
    TH1 *h = (TH1 *)tdir->Get(Form("bin1/Reset/mEtaPhi%i", i));
    st.mValues.push_back(h ? h->GetNbinsX() : 2 * (i + 2));
  }

  for (int n = 1; n <= st.nPtBins; n++) {
    double lo = 0., hi = 0.;
    TH1 *h = (TH1 *)tdir->Get(Form("mFinalFq2Sampled_bin%i", n));
    if (!h || !parseTitleRange(h->GetTitle(), lo, hi)) {
      lo = 0.;
      hi = 0.;
    }
    st.ptLow.push_back(lo);
    st.ptHigh.push_back(hi);
  }

  TString cfg;
  if (tdir->Get("metaConfig")) {
    TH1 *meta = (TH1 *)tdir->Get("metaConfig");
    cfg = meta->GetTitle();
    st.haveConfig = true;
  }

  TString s;
  std::vector<double> vals;
  if (configValue(cfg, "ptCuts", s) && parseList(s, vals) &&
      (int)vals.size() >= 2 * st.nPtBins && st.nPtBins > 0) {
    for (int n = 0; n < st.nPtBins; n++) {
      st.ptLow[n] = vals[2 * n];
      st.ptHigh[n] = vals[2 * n + 1];
    }
  }

  if (configValue(cfg, "centralEta", s))
    st.etaMax = std::atof(s.Data());
  if (configValue(cfg, "vertexXYZ", s)) {
    std::vector<double> vtx;
    if (parseList(s, vtx) && vtx.size() >= 3)
      st.vzMax = vtx[2];
  }
  if (hRef && st.nBins > 0)
    st.nSubsamples = (int)(hRef->GetEntries() / st.nBins);

  if (configValue(cfg, "centLimits", s)) {
    std::vector<double> cent;
    if (parseList(s, cent) && cent.size() >= 2) {
      st.centLow = cent[0];
      st.centHigh = cent[1];
      st.centSource = "metaConfig";
    }
  }

  if (st.centSource.Length() == 0) {
    const char *cand[4] = {"mCentFT0C", "mCentFT0M", "mCentFV0A", "mCentFT0A"};
    for (int i = 0; i < 4 && st.centSource.Length() == 0; i++) {
      double lo = 0., hi = 0.;
      if (filledRange((TH1 *)tdir->Get(cand[i]), lo, hi)) {
        st.centLow = lo;
        st.centHigh = hi;
        st.centSource = Form("observed in %s", cand[i]);
      }
    }
  }

  if (st.etaMax == 0.) {
    double lo = 0., hi = 0.;
    if (filledRange((TH1 *)tdir->Get("mEta"), lo, hi))
      st.etaMax = std::max(std::fabs(lo), std::fabs(hi));
  }

  if (st.vzMax == 0.) {
    double lo = 0., hi = 0.;
    if (filledRange((TH1 *)tdir->Get("mVertexZ"), lo, hi))
      st.vzMax = std::max(std::fabs(lo), std::fabs(hi));
  }

  if (hRef && st.nBins > 0)
    st.nSubsamples = (int)(hRef->GetEntries() / st.nBins);

  return st;
}

static TString ptLabel(const FMSetup &st, int ptBin)
{
  int n = ptBin - 1;
  if (n >= 0 && n < (int)st.ptLow.size() && st.ptHigh[n] > st.ptLow[n])
    return Form("%.2f #leq p_{T} #leq %.2f GeV/c", st.ptLow[n], st.ptHigh[n]);
  return Form("p_{T} bin %i", ptBin);
}

static TString ptText(const FMSetup &st, int ptBin)
{
  int n = ptBin - 1;
  if (n >= 0 && n < (int)st.ptLow.size() && st.ptHigh[n] > st.ptLow[n])
    return Form("%.2f-%.2f", st.ptLow[n], st.ptHigh[n]);
  return Form("bin %i", ptBin);
}

static void processBin(TDirectory *tdir, const FMSetup &st, int ptBin, const char *infile, const char *outdir)
{
  TH1 *hRef = (TH1 *)tdir->Get(Form("mFinalFq2Sampled_bin%i", ptBin));
  if (!hRef) {
    std::cout << "mFinalFq2Sampled_bin" << ptBin << " not found" << std::endl;
    return;
  }

  TString datPath = joinPath(outdir, Form("FactorialMoments_bin%i.dat", ptBin));
  std::ofstream out(datPath.Data());

  if (!out) {
    std::cout << "Cannot write " << datPath << std::endl;
    return;
  }

  TString label = ptLabel(st, ptBin);

  out << "# source: " << infile << std::endl;
  out << "# pT bin " << ptBin << ": " << ptText(st, ptBin) << " GeV/c"
      << " | centrality " << Form("%.0f-%.0f %%", st.centLow, st.centHigh)
      << " | subsamples " << st.nSubsamples;
  if (st.etaMax > 0.)
    out << " | |eta| <= " << Form("%.2f", st.etaMax);
  if (st.vzMax > 0.)
    out << " | |vz| <= " << Form("%.1f", st.vzMax);
  out << std::endl;
  out << "#M2 lnM2";
  for (int j = 0; j < st.nOrders; j++)
    out << " lnF" << (j + 2) << " errF" << (j + 2);
  out << std::endl;
  out << "# err = sigma(Fq)/mean(Fq), i.e. the error on ln(Fq);"
      << " entries with Fq <= 0 are written as -999 -999" << std::endl;
  out << std::fixed << std::setprecision(8);

  std::vector<std::vector<double>> xq(st.nOrders), yq(st.nOrders), eyq(st.nOrders);
  double xMin = 1e30, xMax = -1e30, yMin = 1e30, yMax = -1e30;

  for (int i = 0; i < st.nBins; i++) {
    int mval = st.mValues[i];
    double lnM = std::log((double)mval * mval);
    std::vector<double> vals(st.nOrders), errs(st.nOrders);

    for (int j = 0; j < st.nOrders; j++) {
      int q = j + 2;
      double mean = 0., err = 0.;

      TH1 *hSum = (TH1 *)tdir->Get(Form("mFqSum%i_bin%i", q, ptBin));
      TH1 *hSq = (TH1 *)tdir->Get(Form("mFqSq%i_bin%i", q, ptBin));
      TH1 *hFq = (TH1 *)tdir->Get(Form("mFinalFq%iSampled_bin%i", q, ptBin));

      if (st.haveMoments && hSum && hSq && hFq) {
        int n = (int)(hSum->GetEntries() / st.nBins);
        double sum = hSum->GetBinContent(i + 1);
        double sumSq = hSq->GetBinContent(i + 1);
        mean = (n > 0) ? sum / n : 0.;
        double var = (n > 1) ? (sumSq - n * mean * mean) / (n - 1) : 0.;
        if (var < 0.)
          var = 0.;
        err = (n > 1) ? std::sqrt(var / n) : 0.;
      } else if (hFq) {
        int n = (int)(hFq->GetEntries() / st.nBins);
        mean = (n > 0) ? hFq->GetBinContent(i + 1) / n : 0.;
        TH1 *hErr = (TH1 *)tdir->Get(Form("mFqError%i_bin%i", q, ptBin));
        err = hErr ? hErr->GetBinContent(i + 1) : 0.;
      }

      vals[j] = mean;
      errs[j] = err;
    }

    out << mval * mval << " " << lnM;

    for (int j = 0; j < st.nOrders; j++) {
      bool ok = (vals[j] > 0. && std::isfinite(vals[j]));
      if (ok) {
        out << " " << std::log(vals[j]) << " " << errs[j] / vals[j];
        xq[j].push_back(lnM);
        yq[j].push_back(std::log(vals[j]));
        eyq[j].push_back(errs[j] / vals[j]);
        yMin = std::min(yMin, std::log(vals[j]));
        yMax = std::max(yMax, std::log(vals[j]));
        xMin = std::min(xMin, lnM);
        xMax = std::max(xMax, lnM);
      } else {
        out << " -999 -999";
      }
    }
    out << std::endl;
  }
  out.close();

  if (xMin > xMax) {
    std::cout << "No plottable point for pT bin " << ptBin << " (all Fq <= 0)." << std::endl;
    return;
  }

  double yPad = std::max(0.02, 0.1 * (yMax - yMin));
  double xPad = 0.15;

  TCanvas *c = new TCanvas(Form("c_bin%i", ptBin), "Factorial Moments", 800, 700);

  TH2D *frame = new TH2D(Form("frame_bin%i", ptBin), "", 50, xMin - xPad, xMax + xPad,
                         50, yMin - yPad, yMax + yPad);
  frame->GetXaxis()->SetTitle("ln(M^{2})");
  frame->GetYaxis()->SetTitle("ln(F_{q})");
  frame->Draw();

  int style[6] = {20, 24, 21, 25, 22, 26};
  std::vector<TGraphAsymmErrors *> g(st.nOrders, nullptr);

  for (int j = 0; j < st.nOrders; j++) {
    if (xq[j].empty())
      continue;
    int n = (int)xq[j].size();
    std::vector<double> ex(n, 0.);
    g[j] = new TGraphAsymmErrors(n, &xq[j][0], &yq[j][0], &ex[0], &ex[0],
                                 &eyq[j][0], &eyq[j][0]);
    g[j]->SetMarkerStyle(style[j % 6]);
    g[j]->Draw("P SAME");
  }

  TLegend *leg = new TLegend(0.15, 0.60, 0.35, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  for (int j = 0; j < st.nOrders; j++)
    if (g[j])
      leg->AddEntry(g[j], Form("q = %i", j + 2), "p");
  leg->Draw();

  TPaveText *txt = new TPaveText(0.50, 0.70, 0.88, 0.88, "NDC");
  txt->SetFillStyle(0);
  txt->SetBorderSize(0);
  txt->AddText("Pb-Pb #sqrt{s_{NN}} = 5.36 TeV");
  txt->AddText(label);
  txt->AddText(Form("centrality %.0f-%.0f %%", st.centLow, st.centHigh));
  if (st.etaMax > 0.)
    txt->AddText(Form("|#eta| #leq %.2f", st.etaMax));
  txt->Draw();

  TString pdf = joinPath(outdir, Form("FactorialMoments_bin%i.pdf", ptBin));
  c->SaveAs(pdf.Data());

  std::cout << "pT bin " << ptBin << " (" << label << ")" << std::endl;
  std::cout << "  subsamples : " << st.nSubsamples << std::endl;
  std::cout << "  orders     : q = 2 .. " << (st.nOrders + 1) << std::endl;
  std::cout << "  written    : " << datPath << std::endl;
  std::cout << "               " << pdf << std::endl;
}

void checkFMo2(const char *infile = "/eos/home-s/salman/EPOS_JOBS_BA/BALWAN_O2/results/merged/AnalysisResults.root",
               const char *outdir = "./results")
{
  gStyle->SetOptStat(0);

  TFile *file = TFile::Open(infile);
  if (!file || file->IsZombie()) {
    std::cout << "Cannot open " << infile << std::endl;
    return;
  }

  TDirectory *tdir = (TDirectory *)file->Get("factorial-moments-task");
  if (!tdir) {
    std::cout << "Directory factorial-moments-task not found in " << infile << std::endl;
    return;
  }

  FMSetup st = detectSetup(tdir);
  if (st.nPtBins == 0 || st.nBins == 0) {
    std::cout << "No factorial-moments histograms found in " << infile << std::endl;
    return;
  }

  gSystem->mkdir(outdir, kTRUE);

  std::cout << "file        : " << infile << std::endl;
  std::cout << "pT bins     : " << st.nPtBins;
  for (int n = 0; n < st.nPtBins; n++)
    std::cout << "   [" << ptText(st, n + 1) << "]";
  std::cout << std::endl;
  std::cout << "M bins      : " << st.nBins << "  (M = " << st.mValues.front() << " .. " << st.mValues.back() << ")" << std::endl;
  std::cout << "orders      : q = 2 .. " << (st.nOrders + 1) << std::endl;
  std::cout << "centrality  : " << Form("%.0f-%.0f %%", st.centLow, st.centHigh)
            << "  (from " << st.centSource << ")" << std::endl;
  if (st.etaMax > 0.)
    std::cout << "eta         : |eta| <= " << Form("%.2f", st.etaMax) << std::endl;
  if (st.vzMax > 0.)
    std::cout << "vertex z    : |vz| <= " << Form("%.1f", st.vzMax) << std::endl;
  std::cout << "subsamples  : " << st.nSubsamples << std::endl;
  if (!st.haveMoments)
    std::cout << "[WARN] mFqSum/mFqSq not in this file (older task): errors are correct for ONE job only." << std::endl;
  std::cout << "output dir  : " << outdir << std::endl;

  for (int ptBin = 1; ptBin <= st.nPtBins; ptBin++)
    processBin(tdir, st, ptBin, infile, outdir);
}
