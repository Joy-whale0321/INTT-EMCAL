#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>

void checkfile(const char* filename = "electron_pt0-10GeV_100k.root", 
                      const char* treename = "tree", 
                      const char* branchname = "trk_layer", 
                      int nbins = 50, 
                      double xmin = -0.5, 
                      double xmax = 10)
{
  // 打开 ROOT 文件
  TFile* file = TFile::Open(filename);
  if (!file || file->IsZombie())
  {
    std::cerr << "Error: Cannot open file " << filename << std::endl;
    return;
  }

  // 获取 TTree
  TTree* tree = dynamic_cast<TTree*>(file->Get(treename));
  if (!tree)
  {
    std::cerr << "Error: Cannot find tree " << treename << std::endl;
    return;
  }

  // 设置 branch 地址（以 vector<float> 为例）
  std::vector<float>* vec = nullptr;
  tree->SetBranchAddress(branchname, &vec);

  // 创建直方图
  TH1D* h = new TH1D("h_vec_branch", Form("%s distribution", branchname), nbins, xmin, xmax);

  // 遍历所有 entries
  Long64_t nentries = tree->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i)
  {
    tree->GetEntry(i);
    for (float val : *vec)
    {
      h->Fill(val);
    }
  }

  // 画图
  TCanvas* c1 = new TCanvas("c1", "Vector Branch Histogram", 800, 600);
  h->GetXaxis()->SetTitle(branchname);
  h->GetYaxis()->SetTitle("Counts");
  h->Draw();

  c1->SaveAs("branchcheck.pdf");
}
