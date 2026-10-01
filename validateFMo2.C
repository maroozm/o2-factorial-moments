/// \file validateFMo2.C
/// \brief Cross-check the checkFMo2.C output against an independent pooled recomputation from the per-file results
/// \author Salman Malik
/// \author Balwan Singh

#include <cstdio>
#include <cstdlib>
#include <vector>

#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TMath.h>
#include <TString.h>
#include <TSystem.h>

constexpr int kNRowsMax = 64;
constexpr int kNFiles = 4;
constexpr int kNPtBins = 3;
constexpr int kNOrdersMax = 20;

static TString dataPath(const char *outdir, int ptBin)
{
  TString p(outdir);
  if (!p.BeginsWith("/")) {
    TString base = gSystem->Getenv("PWD");
    if (base.Length() == 0)
      base = ".";
    p = base + "/" + p;
  }
  if (!p.EndsWith("/"))
    p += "/";
  p += Form("FactorialMoments_bin%i.dat", ptBin);
  return p;
}

static bool readRows(const TString &path, std::vector<std::vector<double>> &rows)
{
  FILE *f = fopen(path.Data(), "r");
  if (!f)
    return false;
  char line[8192];
  rows.clear();
  while (fgets(line, sizeof line, f) && (int)rows.size() < kNRowsMax) {
    const char *p = line;
    if (*p == '#' || *p == '\n' || *p == '\0')
      continue;
    std::vector<double> vals;
    while (*p != '\0') {
      while (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r')
        p++;
      if (*p == '\0')
        break;
      char *end = nullptr;
      double v = strtod(p, &end);
      if (end == p)
        break;
      vals.push_back(v);
      p = end;
    }
    if (vals.size() >= 4)
      rows.push_back(vals);
  }
  fclose(f);
  return !rows.empty();
}

void validateFMo2(const char *outdir = "./results")
{
  const char *tags[kNFiles] = {"AO2D_001", "AO2D_002", "AO2D_003", "AO2D_005"};
  TString base = gSystem->Getenv("PWD");
  if (base.Length() == 0)
    base = ".";
  base += "/";

  TFile *files[kNFiles];
  for (int t = 0; t < kNFiles; t++) {
    files[t] = TFile::Open(base + Form("results/%s/AnalysisResults.root", tags[t]));
    if (!files[t] || files[t]->IsZombie()) {
      printf("cannot open results/%s/AnalysisResults.root\n", tags[t]);
      return;
    }
  }
  TFile *fm = TFile::Open(base + "results/merged/AnalysisResults.root");
  if (!fm || fm->IsZombie()) {
    printf("cannot open results/merged/AnalysisResults.root\n");
    return;
  }
  TDirectory *dm = (TDirectory *)fm->Get("factorial-moments-task");

  double maxRelMean = 0, maxRelErr = 0;

  for (int ptBin = 1; ptBin <= kNPtBins; ptBin++) {
    TString datPath = dataPath(outdir, ptBin);
    std::vector<std::vector<double>> rows;
    if (!readRows(datPath, rows)) {
      printf("%s missing or empty - run checkFMo2 first\n", datPath.Data());
      return;
    }

    int nOrders = ((int)rows[0].size() - 2) / 2;
    if (nOrders < 1 || nOrders > kNOrdersMax) {
      printf("%s: unexpected column count %i\n", datPath.Data(), (int)rows[0].size());
      return;
    }

    TDirectory *d0 = (TDirectory *)files[0]->Get("factorial-moments-task");
    TH1 *hRef = (TH1 *)d0->Get(Form("mFinalFq2Sampled_bin%i", ptBin));
    int nBins = hRef ? hRef->GetNbinsX() : 52;

    std::vector<int> mValues(nBins, 0);
    for (int i = 0; i < nBins; i++) {
      TH1 *h = (TH1 *)d0->Get(Form("bin1/Reset/mEtaPhi%i", i));
      mValues[i] = h ? h->GetNbinsX() : 2 * (i + 2);
    }

    double relMean = 0, relErr = 0;
    int nCmp = 0;

    for (int i = 0; i < (int)rows.size(); i++) {
      int b = i + 1;
      if (b > nBins)
        break;
      double m2 = rows[i][0];
      if (TMath::Abs(m2 - (double)mValues[i] * mValues[i]) > 0.5)
        continue;
      if (TMath::Abs(rows[i][1] - TMath::Log(m2)) > 1e-6)
        continue;

      for (int j = 0; j < nOrders; j++) {
        int q = j + 2;
        double lnF = rows[i][2 + 2 * j];
        double errF = rows[i][3 + 2 * j];
        if (lnF < -900.)
          continue;

        double num = 0, Ntot = 0, var = 0;
        double mf[kNFiles], nf[kNFiles], sf[kNFiles];
        for (int t = 0; t < kNFiles; t++) {
          TDirectory *d = (TDirectory *)files[t]->Get("factorial-moments-task");
          TH1 *h = (TH1 *)d->Get(Form("mFinalFq%iSampled_bin%i", q, ptBin));
          TH1 *e = (TH1 *)d->Get(Form("mFqError%i_bin%i", q, ptBin));
          if (!h) {
            printf("missing mFinalFq%iSampled_bin%i in %s\n", q, ptBin, tags[t]);
            return;
          }
          int n = (int)TMath::Nint(h->GetEntries() / nBins);
          num += h->GetBinContent(b);
          Ntot += n;
          mf[t] = n ? h->GetBinContent(b) / n : 0;
          nf[t] = n;
          sf[t] = (e && n > 0) ? e->GetBinContent(b) * TMath::Sqrt((double)n) : 0;
        }
        if (Ntot <= 1)
          continue;
        double pooledMean = num / Ntot;
        for (int t = 0; t < kNFiles; t++) {
          if (nf[t] > 1)
            var += (nf[t] - 1) * sf[t] * sf[t];
          var += nf[t] * (mf[t] - pooledMean) * (mf[t] - pooledMean);
        }
        double pooledErr = TMath::Sqrt(var / (Ntot - 1) / Ntot);

        double datMean = TMath::Exp(lnF);
        double datErr = errF * datMean;

        relMean = TMath::Max(relMean, TMath::Abs(datMean - pooledMean) / TMath::Abs(pooledMean));
        if (pooledErr > 0)
          relErr = TMath::Max(relErr, TMath::Abs(datErr - pooledErr) / pooledErr);
        nCmp++;

        TH1 *hs = (TH1 *)dm->Get(Form("mFqSum%i_bin%i", q, ptBin));
        TH1 *hq = (TH1 *)dm->Get(Form("mFqSq%i_bin%i", q, ptBin));
        if (hs && hq) {
          int nsm = (int)TMath::Nint(hs->GetEntries() / nBins);
          if (nsm > 1) {
            double sum = hs->GetBinContent(b), sumSq = hq->GetBinContent(b);
            double mean = sum / nsm;
            double rawErr = TMath::Sqrt(((sumSq - nsm * mean * mean) / (nsm - 1)) / nsm);
            if (pooledErr > 0)
              relErr = TMath::Max(relErr, TMath::Abs(rawErr - pooledErr) / pooledErr);
          }
        }
      }
    }

    printf("pt bin %i : rows=%3d  orders=q2..q%i  compared=%4d   max rel. dev.  value=%.3e  error=%.3e\n",
           ptBin, (int)rows.size(), nOrders + 1, nCmp, relMean, relErr);
    maxRelMean = TMath::Max(maxRelMean, relMean);
    maxRelErr = TMath::Max(maxRelErr, relErr);
  }

  for (int t = 0; t < kNFiles; t++)
    files[t]->Close();
  fm->Close();

  printf("\nOVERALL: max rel. dev.  value=%.3e   error=%.3e\n", maxRelMean, maxRelErr);
}
