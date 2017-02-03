#include "DQMFlashAlg.h"

#include <functional>
#include <unordered_map>

dqm::DQMFlashAlg::DQMFlashAlg() {
}

void dqm::DQMFlashAlg::AnalyzeFlashes(std::vector<recob::OpFlash> flashVector, TH1F *NFlashes, TH1F *MeanFlashLight, TH1F *FlashLightVar, TH1F *FlashLightSkew)
{
  std::cout << "number of flashes: " << flashVector.size() << std::endl;
  NFlashes->Fill(flashVector.size());
  int num_pmts = 32;
  TH1F *flashlighthist = new TH1F("flashlighthist", "flashlighthist", 1E5, 0, 1E5);
  for (unsigned int i = 0; i < flashVector.size(); i++)
  {
    for (int pmt_index = 0; pmt_index < num_pmts; pmt_index++)
    {
      if (flashVector[i].PE(pmt_index) != 0)
        flashlighthist->Fill(flashVector[i].PE(pmt_index));
    }
  }
  MeanFlashLight->Fill(flashlighthist->GetMean());
  FlashLightVar->Fill((flashlighthist->GetRMS())*(flashlighthist->GetRMS()));
  FlashLightSkew->Fill(flashlighthist->GetSkewness());
  std::cout << "Maximum bin: " << flashlighthist->GetMaximumBin() << std::endl;
}

