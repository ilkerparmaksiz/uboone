#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Services/Optional/TFileService.h"

#include <vector>
#include <string>

#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/MCBase/MCHitCollection.h"
#include "uboone/DQMTools/DQMFlashAlg.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "lardataobj/RecoBase/OpFlash.h"

#include "TH1F.h"

namespace dqm {
  class DQMFlashModule;
}

class dqm::DQMFlashModule : public art::EDAnalyzer {
public:
  explicit DQMFlashModule(fhicl::ParameterSet const & p);
  virtual ~DQMFlashModule();

  void analyze(art::Event const & e) override;

  void beginJob() override;

private:

  TH1F *NFlashes;

  TH1F *MeanFlashLight;
  TH1F *FlashLightVar;
  TH1F *FlashLightSkew;

  DQMFlashAlg analysisAlg;

};


dqm::DQMFlashModule::DQMFlashModule(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{ 
  //this->reconfigure(p); 
}

dqm::DQMFlashModule::~DQMFlashModule()
{
  // Clean up dynamic memory and other resources here.
}

void dqm::DQMFlashModule::analyze(art::Event const & e)
{
  //get the flash data
  art::Handle< std::vector<recob::OpFlash> > flashHandle;
  std::string _flash_producer_name = "opflashSat";
  e.getByLabel(_flash_producer_name,flashHandle);
  if(!flashHandle.isValid()) {
      std::cerr << "\033[93m[ERROR]\033[00m Could not retrieve recob::OpFlash from "
      << _flash_producer_name << std::endl;
      throw std::exception();
  }
  std::vector<recob::OpFlash> const& flashVector(*flashHandle);
 
  analysisAlg.AnalyzeFlashes(flashVector, NFlashes, MeanFlashLight, FlashLightVar, FlashLightSkew);
}

void dqm::DQMFlashModule::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;

  NFlashes       = tfs->make<TH1F>("NFlashes"      , "# of flashes"        , 300, 0, 300);
  MeanFlashLight = tfs->make<TH1F>("MeanFlashLight", "Mean light per flash", 100, 0, 500);
  FlashLightVar  = tfs->make<TH1F>("FlashLightVar" , "Variance on light per flash", 1E3, 0, 1E6);
  FlashLightSkew = tfs->make<TH1F>("FlashLightSkew", "Skewness on light per flash", 100, -500, 500);
 
}

DEFINE_ART_MODULE(dqm::DQMFlashModule)
