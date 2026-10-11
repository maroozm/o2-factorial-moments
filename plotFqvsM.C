/// \file plotFqvsM.C
/// \brief Draw F_q vs M^2 from a slotpool.C CSV -- one canvas per pt bin with
///        ALL orders overlaid.
///
/// slotpool.C writes the pooled new-path result as a CSV (pt,q,M,nGroups,
/// meanNew,errNew behind a '#' header block that also carries cent/eta/vertex).
/// This macro reads that CSV directly -- no .dat involved -- and draws the
/// established ATLAS-style figure: every order present (q=2..7) as its own
/// marker+band series on a single canvas per pt bin, identified by a legend.
///
/// Drawing follows plotdat.C exactly: gROOT->SetStyle("ATLAS"), raw M^2 on x
/// and raw F_q on y, the per-q colour/marker table, translucent error bands,
/// and the NDC header block at the approved geometry (four lines: energy, pT
/// range, centrality, eta+vertex cuts parsed from the CSV's '# eta' /
/// '# vertex' lines), with the legend of orders above the curves.
///
///   # all pt bins, auto y range (log y auto since q=6,7 span decades)
///   root -l -b -q 'plotFqvsM.C+("ab_scan/S17/slotpool_K17.csv")'
///
///   # only pt bins 1 and 3; force linear y; y range 0.9..1.4; custom out dir
///   root -l -b -q 'plotFqvsM.C+("file.csv", "1,3", 1, 0.9, 1.4, "mydir")'
///
/// Arguments:
///   csvFile  slotpool CSV (new format; the older 8-column CSV with old-path
///            columns is also accepted -- those columns are simply ignored)
///   ptSel    "all" or a comma list, e.g. "1,3"
///   logMode  0 = auto (log y when q=6,7 are present, i.e. the combined F_q
///            range spans decades), 1 = force linear y, 2 = force log y
///   ymin     y-axis min, <= 0 = choose automatically (established margins)
///   ymax     y-axis max, <= 0 = choose automatically
///   outDir   output directory for the PDFs, "" = <csv dir>/plots_FqvsM
///
/// One PDF per pt bin: FactorialMoments_bin<p>.pdf. x = M^2 is rebuilt from
/// the CSV's M index through the lattice: the cell counts of
/// bin1/Reset/mEtaPhi{i} in the root file named by the CSV's '# source' line
/// (the same object slotpool.C's latticeM() reads). If that file is not
/// reachable the lattice falls back to 2*(i+2) WITH A WARNING -- the source
/// root should normally be next to the CSV.
///
/// \author Salman Malik
/// \author Balwan Singh

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2D.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

namespace
{

// ---------------------------------------------------------------------------
// Style -- identical to plotdat.C (which itself lifted it from the
// reference macro): same ATLAS base, same paper-size match, same margins.
// ---------------------------------------------------------------------------
struct QStyle {
  Color_t col;
  Marker_t mk;
  Size_t sz;
  const char* lab;
  bool logy; // set where the F_q range spans a decade
};

const QStyle kQ[6] = {
  {kGreen + 2,  24, 1.43, "#it{q}=2", false}, // green open circle
  {kBlue,       34, 1.71, "#it{q}=3", false}, // blue open cross
  {kRed + 1,    20, 1.33, "#it{q}=4", false}, // red filled circle
  {kYellow + 1, 33, 2.00, "#it{q}=5", false}, // yellow filled diamond
  {kMagenta + 1,25, 1.50, "#it{q}=6", true }, // magenta open square
  {kOrange + 2, 21, 1.50, "#it{q}=7", true }, // orange filled square
};

// See plotdat.C: SetPaperSize must match the canvas aspect or ROOT
// letterboxes the PDF onto A4.
void applyAtlas(int wPx = 800, int hPx = 600)
{
  gROOT->SetStyle("ATLAS");
  gStyle->SetPaperSize(wPx / 72.0 * 2.54, hPx / 72.0 * 2.54);
  gStyle->SetTitleSize(0.046, "X");
  gStyle->SetTitleSize(0.046, "Y");
  gStyle->SetLabelSize(0.045, "X");
  gStyle->SetLabelSize(0.045, "Y");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetPadLeftMargin(0.145);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetPadTopMargin(0.06);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetOptStat(0);
}

// "0.20-2.00 GeV/c" -> "0.20 #leq p_{T} #leq 2.00 GeV/c"
std::string ptLabel(const std::string& range)
{
  auto dash = range.find('-');
  if (dash == std::string::npos)
    return range;
  return range.substr(0, dash) + " #leq p_{T} #leq " + range.substr(dash + 1);
}

// substring between two markers, "" if either is missing; b == "" -> to end
std::string between(const std::string& s, const std::string& a, const std::string& b)
{
  auto p = s.find(a);
  if (p == std::string::npos)
    return {};
  p += a.size();
  if (b.empty())
    return s.substr(p);
  auto e = s.find(b, p);
  if (e == std::string::npos)
    return s.substr(p);
  return s.substr(p, e - p);
}

// ---------------------------------------------------------------------------
// CSV input
// ---------------------------------------------------------------------------
struct Meta {
  std::string src;     // '# source'  root file the CSV came from
  int K = 0;           // '# K'       target subsample size
  int subsamples = 0;  // '# K'       groups pooled
  std::string cent;    // '# cent'    e.g. "0-10 %"
  std::string eta;     // '# eta'     e.g. "|eta| <= 0.90"
  std::string vertex;  // '# vertex'  e.g. "|vz| <= 10.0 cm"
  std::map<int, std::string> ptRange; // '# pt bins'  pt bin -> "0.20-2.00 GeV/c"

  std::string rangeOf(int p) const
  {
    auto it = ptRange.find(p);
    return it == ptRange.end() ? std::string() : it->second;
  }
};

void parseMeta(const std::string& line, Meta& m)
{
  if (line.rfind("# source", 0) == 0) {
    m.src = between(line, ": ", "");
  } else if (line.rfind("# K", 0) == 0) {
    const std::string rest = between(line, ": ", "");
    m.K = atoi(between(rest, "", " ").c_str());
    m.subsamples = atoi(between(rest, "subsamples: ", "").c_str());
  } else if (line.rfind("# cent", 0) == 0) {
    m.cent = between(line, ": ", "");
  } else if (line.rfind("# eta", 0) == 0) {
    m.eta = between(line, ": ", "");
  } else if (line.rfind("# vertex", 0) == 0) {
    m.vertex = between(line, ": ", "");
  } else if (line.rfind("# pt bins", 0) == 0) {
    // "1: 0.20-2.00 GeV/c  2: 0.40-2.00 GeV/c  3: ..."
    const std::string s = between(line, ": ", "");
    size_t pos = 0;
    while (pos < s.size()) {
      while (pos < s.size() && isspace(s[pos])) ++pos;
      int p = 0;
      bool any = false;
      while (pos < s.size() && isdigit(s[pos])) {
        p = 10 * p + (s[pos] - '0');
        ++pos;
        any = true;
      }
      if (!any) break;
      if (pos < s.size() && s[pos] == ':') ++pos;
      while (pos < s.size() && s[pos] == ' ') ++pos;
      const size_t e = s.find("  ", pos);
      std::string range = (e == std::string::npos) ? s.substr(pos) : s.substr(pos, e - pos);
      while (!range.empty() && isspace(range.back())) range.pop_back();
      if (!range.empty()) m.ptRange[p] = range;
      pos = (e == std::string::npos) ? s.size() : e;
    }
  }
}

struct Row {
  int pt = 0;
  int q = 0;
  int M = 0; // 1-based M index
  double mean = 0.;
  double err = 0.;
};

// Accepts both the new 6-column CSV and the older 8-column one (old-path
// columns are ignored): the column header line names the fields, so we map by
// name instead of by position.
bool readCsv(const char* path, Meta& meta, std::vector<Row>& rows)
{
  std::ifstream in(path);
  if (!in)
    return false;
  std::string line;
  int cPt = -1, cQ = -1, cM = -1, cNew = -1, cErr = -1;
  while (std::getline(in, line)) {
    if (line.empty())
      continue;
    if (line[0] == '#') {
      parseMeta(line, meta);
      continue;
    }
    if (line.rfind("pt,", 0) == 0) {
      std::stringstream ss(line);
      std::string t;
      int i = 0;
      while (std::getline(ss, t, ',')) {
        if (t == "pt") cPt = i;
        else if (t == "q") cQ = i;
        else if (t == "M") cM = i;
        else if (t == "meanNew") cNew = i;
        else if (t == "errNew") cErr = i;
        ++i;
      }
      continue;
    }
    if (cNew < 0 || cErr < 0 || cPt < 0 || cQ < 0 || cM < 0)
      continue; // data row before the column header
    std::stringstream ss(line);
    std::string t;
    std::vector<std::string> v;
    while (std::getline(ss, t, ','))
      v.push_back(t);
    const int need = std::max(std::max(cPt, cQ), std::max(cM, std::max(cNew, cErr)));
    if (static_cast<int>(v.size()) <= need)
      continue;
    Row r;
    r.pt = atoi(v[cPt].c_str());
    r.q = atoi(v[cQ].c_str());
    r.M = atoi(v[cM].c_str());
    r.mean = atof(v[cNew].c_str());
    r.err = atof(v[cErr].c_str());
    rows.push_back(r);
  }
  return !rows.empty();
}

// cell counts per M index, from the source root file -- the exact object
// slotpool.C's latticeM() reads; fallback 2*(i+2) with a warning (the same
// default checkFMo2.C/slotpool.C use, and what this production lattice is).
bool latticeFromRoot(const std::string& src, int maxM, std::vector<int>& mv)
{
  if (src.empty() || maxM <= 0)
    return false;
  const TString p = gSystem->ExpandPathName(src.c_str());
  if (gSystem->AccessPathName(p, kFileExists))
    return false;
  TFile* f = TFile::Open(p, "READ");
  if (!f || f->IsZombie()) {
    delete f;
    return false;
  }
  TDirectory* td = dynamic_cast<TDirectory*>(f->Get("factorial-moments-task"));
  if (!td) {
    TIter next(f->GetListOfKeys());
    while (auto* k = (TKey*)next()) {
      if (std::string(k->GetClassName()).find("TDirectory") == std::string::npos)
        continue;
      auto* d = dynamic_cast<TDirectory*>(f->Get(k->GetName()));
      if (d && (d->Get("slots/nEv") || d->Get("mSubsampleEvents"))) {
        td = d;
        break;
      }
    }
  }
  if (!td) {
    f->Close();
    return false;
  }
  mv.clear();
  for (int i = 0; i < maxM; ++i) {
    auto* h = dynamic_cast<TH1*>(td->Get(Form("bin1/Reset/mEtaPhi%d", i)));
    mv.push_back(h ? h->GetNbinsX() : 2 * (i + 2));
  }
  f->Close();
  return true;
}

// ---------------------------------------------------------------------------
// Drawing
// ---------------------------------------------------------------------------

// Five-line NDC header block at plotdat.C's geometry (x=0.185, 0.055 down from
// the original grid). Line 3 is q-free -- the legend identifies the orders --
// and line 4 shows the eta + vertex cuts from the CSV.
void drawHeader(TCanvas* c, const Meta& meta, int ptBin)
{
  c->cd();
  TLatex latex;
  latex.SetNDC();
  latex.SetTextFont(42);

  latex.SetTextColor(kAzure - 1);
  latex.SetTextSize(0.042);
  latex.DrawLatex(0.185, 0.880, "Pb-Pb #sqrt{s_{NN}} = 5.36 TeV");
  const std::string range = meta.rangeOf(ptBin);
  if (!range.empty())
    latex.DrawLatex(0.185, 0.823, ptLabel(range).c_str());
  else if (!meta.cent.empty())
    latex.DrawLatex(0.185, 0.823, Form("centrality %s", meta.cent.c_str()));

  latex.SetTextColor(kBlack);
  latex.SetTextSize(0.030);
  if (!meta.cent.empty())
    latex.DrawLatex(0.185, 0.771, Form("centrality %s", meta.cent.c_str()));

  // "|eta| <= 0.90" / "|vz| <= 10.0 cm" -> "|#eta| #leq 0.90,  |v_{z}| #leq 10.0"
  const std::string eta = between(meta.eta, "|eta| <= ", "");
  const std::string vz = between(meta.vertex, "|vz| <= ", " cm");
  std::string l4;
  if (!eta.empty() || !vz.empty()) {
    if (!eta.empty())
      l4 += Form("|#eta| #leq %s", eta.c_str());
    if (!eta.empty() && !vz.empty())
      l4 += ",  ";
    if (!vz.empty())
      l4 += Form("|v_{z}| #leq %s", vz.c_str());
  }
  if (!l4.empty())
    latex.DrawLatex(0.185, 0.727, l4.c_str());
}

// One canvas per pt bin: every order present drawn as markers + translucent
// band, y range spanning all of them, legend identifying the orders.
void drawPtCanvas(const Meta& meta, int ptBin, const std::vector<int>& orders,
                  const std::map<int, std::vector<Row>>& rows, const std::vector<int>& mv,
                  bool forceLinY, bool forceLogY, double yminUser, double ymaxUser,
                  const TString& out)
{
  struct Series {
    std::vector<double> x, y, ey;
  };
  std::map<int, Series> data;
  std::vector<int> drawn;
  double yLo = 1e30, yHi = -1e30, yLoPos = 1e30, xHi = 0.0;
  int skippedM = 0;

  for (int q : orders) {
    const int iq = q - 2;
    if (iq < 0 || iq >= 6) {
      std::cout << "  pT bin " << ptBin << ": order q=" << q << " outside 2..7, skipped" << std::endl;
      continue;
    }
    auto it = rows.find(q);
    if (it == rows.end())
      continue;
    Series& s = data[q];
    for (const auto& r : it->second) {
      if (r.M < 1 || r.M > static_cast<int>(mv.size())) {
        ++skippedM;
        continue;
      }
      if (!(r.mean > 0.) || !std::isfinite(r.mean))
        continue;
      const double m2 = static_cast<double>(mv[r.M - 1]) * mv[r.M - 1];
      const double efq = (r.err > 0. && std::isfinite(r.err)) ? r.err : 0.;
      s.x.push_back(m2);
      s.y.push_back(r.mean);
      s.ey.push_back(efq);
      const double lo = r.mean - efq, hi = r.mean + efq;
      yLo = std::min(yLo, lo);
      yHi = std::max(yHi, hi);
      // separate strictly-positive lower edge for the log case (see plotdat.C)
      yLoPos = std::min(yLoPos, (lo > 0.) ? lo : r.mean);
      xHi = std::max(xHi, m2);
    }
    if (s.y.empty())
      data.erase(q); // e.g. an order written entirely as zeros (Fq <= 0)
    else
      drawn.push_back(q);
  }
  if (skippedM)
    std::cout << "  !! pT bin " << ptBin << ": " << skippedM
              << " rows with M outside the lattice (1.." << mv.size() << ")" << std::endl;
  if (drawn.empty()) {
    std::cout << "  pT bin " << ptBin << ": nothing plottable" << std::endl;
    return;
  }

  // ---- y range across all drawn orders ----------------------------------
  bool anyLogy = false;
  for (int q : drawn)
    anyLogy = anyLogy || kQ[q - 2].logy;
  const bool spanWide = (yLo > 0. && yHi / yLo > 100.);
  const bool logY = forceLogY || (!forceLinY && (anyLogy || spanWide));
  if (logY) {
    yLo = yLoPos;
    if (!(yLo > 0.))
      yLo = 0.5 * yHi;
    if (!(yHi > yLo))
      yHi = 2. * yLo;
    double span = std::log(yHi / yLo);
    span = std::max(span, 0.05);
    const double mid = 0.5 * (std::log(yLo) + std::log(yHi));
    yLo = std::exp(mid - 0.5 * span - 0.18 * span);
    yHi = std::exp(mid + 0.5 * span + 0.18 * span);
  } else {
    if (!(yHi > yLo)) {
      yLo -= 0.05;
      yHi += 0.05;
    }
    double span = yHi - yLo;
    span = std::max(span, 0.05 * std::max(std::fabs(yHi), 0.01));
    const double mid = 0.5 * (yLo + yHi);
    yLo = mid - 0.5 * span - 0.18 * span;
    yHi = mid + 0.5 * span + 0.18 * span;
  }
  // user overrides (either one or both; <= 0 keeps the automatic edge)
  if (yminUser > 0.)
    yLo = yminUser;
  if (ymaxUser > 0.)
    yHi = ymaxUser;
  if (!(yHi > yLo)) {
    std::cout << "  pT bin " << ptBin << ": empty y range after override, skipped" << std::endl;
    return;
  }

  const double xHiFrame = std::max(11000., xHi * 1.05); // 0..11000 as reference

  // ---- frame ------------------------------------------------------------
  TCanvas* c = new TCanvas(Form("c_bin%i", ptBin), "Factorial Moments", 800, 600);
  if (logY)
    c->SetLogy();
  TH2D* frame = new TH2D(Form("frame_bin%i", ptBin), "", 1000, 0.0, xHiFrame, 200, yLo, yHi);
  frame->SetXTitle("#it{M}^{2}");
  frame->SetYTitle("#it{F}_{q}(#it{M})");
  frame->GetYaxis()->SetMoreLogLabels();
  frame->Draw();

  // ---- series (markers + translucent band, as plotdat's single series) --
  // legend above the curves, top-right: box bottom sits just above the
  // highest band edge in both the flat (linear) and the all-orders (log)
  // canvases; the header block owns the top-left.
  auto* leg = new TLegend(0.61, 0.845, 0.97, 0.94);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.026);
  leg->SetTextFont(42);
  leg->SetNColumns(2);

  for (int q : drawn) {
    const QStyle& st = kQ[q - 2];
    Series& s = data[q];
    const int n = static_cast<int>(s.x.size());
    auto* g = new TGraphAsymmErrors(n, &s.x[0], &s.y[0], nullptr, nullptr, &s.ey[0], &s.ey[0]);
    g->SetMarkerColorAlpha(st.col, 0.95);
    g->SetMarkerStyle(st.mk);
    g->SetMarkerSize(st.sz);
    g->SetLineColorAlpha(st.col, 0.95);
    g->SetFillColorAlpha(st.col, 0.20);
    g->SetFillStyle(1001);
    g->SetLineWidth(1);
    g->Draw("PE3 SAME");
    leg->AddEntry(g, st.lab, "pf");
  }
  leg->Draw();

  drawHeader(c, meta, ptBin);

  const TString pdf = out + Form("/FactorialMoments_bin%d.pdf", ptBin);
  c->SaveAs(pdf.Data());
  const std::string range = meta.rangeOf(ptBin);
  std::cout << "  pT bin " << ptBin
            << (range.empty() ? "" : " (" + range + ")") << "  q=";
  for (size_t i = 0; i < drawn.size(); ++i)
    std::cout << (i ? "," : "") << drawn[i];
  std::cout << (logY ? "  [log y]" : "  [lin y]") << std::endl;
  std::cout << "             " << pdf << std::endl;

  delete c;
}

} // namespace

void plotFqvsM(const char* csvFile, const char* ptSel = "all", int logMode = 0,
               double ymin = 0., double ymax = 0., const char* outDir = "")
{
  applyAtlas();

  Meta meta;
  std::vector<Row> rows;
  if (!readCsv(csvFile, meta, rows)) {
    std::cout << "cannot read a slotpool CSV from " << csvFile << std::endl;
    return;
  }

  // pt selection: "all" or comma list
  std::vector<int> want;
  if (ptSel && ptSel[0] && std::string(ptSel) != "all") {
    std::stringstream ss(ptSel);
    std::string t;
    while (std::getline(ss, t, ','))
      if (!t.empty())
        want.push_back(atoi(t.c_str()));
  }

  // pt bins present, in ascending order; rows grouped per (pt bin, order)
  std::vector<int> ptBins;
  std::map<int, std::map<int, std::vector<Row>>> byPtQ;
  int maxM = 0;
  for (const auto& r : rows) {
    if (std::find(ptBins.begin(), ptBins.end(), r.pt) == ptBins.end())
      ptBins.push_back(r.pt);
    byPtQ[r.pt][r.q].push_back(r);
    maxM = std::max(maxM, r.M);
  }
  std::sort(ptBins.begin(), ptBins.end());

  // lattice: cell counts per M index (same object slotpool.C uses)
  std::vector<int> mv;
  if (!latticeFromRoot(meta.src, maxM, mv)) {
    std::cout << "!! source root not reachable (" << meta.src
              << ") -- lattice FALLBACK 2*(i+2): x = M^2 uses the default cell\n"
              << "   counts (matches this production lattice, but the source root is\n"
              << "   authoritative -- run next to the CSV's AnalysisResults.root)\n";
    mv.clear();
    for (int i = 0; i < maxM; ++i)
      mv.push_back(2 * (i + 2));
  }

  TString out = (outDir && outDir[0]) ? TString(outDir) : TString(csvFile);
  if (!(outDir && outDir[0])) {
    const Ssiz_t sl = out.Last('/');
    out = (sl == kNPOS) ? TString(".") : out(0, sl);
    out += "/plots_FqvsM";
  }
  gSystem->mkdir(out, kTRUE);

  std::cout << "=== plotFqvsM ===\n";
  std::cout << "  csv     : " << csvFile << std::endl;
  if (meta.K > 0 || meta.subsamples > 0)
    std::cout << "  pooling : K=" << meta.K << "   subsamples=" << meta.subsamples << std::endl;
  if (!meta.cent.empty())
    std::cout << "  cent    : " << meta.cent << std::endl;
  if (!meta.eta.empty())
    std::cout << "  eta     : " << meta.eta << std::endl;
  if (!meta.vertex.empty())
    std::cout << "  vertex  : " << meta.vertex << std::endl;
  std::cout << "  pt bins : ";
  for (size_t i = 0; i < ptBins.size(); ++i)
    std::cout << (i ? ", " : "") << ptBins[i];
  if (!want.empty()) {
    std::cout << "   -> plotting only ";
    for (size_t i = 0; i < want.size(); ++i)
      std::cout << (i ? ", " : "") << want[i];
  }
  std::cout << std::endl;
  std::cout << "  output  : " << out << "\n" << std::endl;

  const bool forceLinY = (logMode == 1);
  const bool forceLogY = (logMode == 2);

  for (int p : ptBins) {
    if (!want.empty() && std::find(want.begin(), want.end(), p) == want.end())
      continue;
    std::map<int, std::vector<Row>>& byQ = byPtQ[p];
    std::vector<int> orders;
    for (const auto& kv : byQ)
      orders.push_back(kv.first);
    std::sort(orders.begin(), orders.end());
    drawPtCanvas(meta, p, orders, byQ, mv, forceLinY, forceLogY, ymin, ymax, out);
  }
}
