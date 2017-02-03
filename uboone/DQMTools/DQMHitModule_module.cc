////////////////////////////////////////////////////////////////////////
// Class:       DQMHitModule
// Module Type: analyzer
// File:        DQMHitModule_module.cc
//
// Generated at Sun Jun 22 08:40:08 2014 by Wesley Ketchum using artmod
// from cetpkgsupport v1_05_04.
////////////////////////////////////////////////////////////////////////

/*!
 * Title:   DQMHitModule
 * Author:  wketchum@lanl.gov
 * Inputs:  recob::Wire (calibrated), recob::Hit, Assns<recob::Wire, recob::Hit>
 * Outputs: validation histograms
 *
 * Description:
 * This module is intended to be yet another hit analyzer module. Its intention is
 *   (1) to compare hit-finding modules against each other, and eventually
 *   (2) to compare those to truth
 */

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
#include "uboone/DQMTools/DQMHitAlg.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "lardataobj/RecoBase/OpFlash.h"

#include "TH1F.h"

namespace dqm {
  class DQMHitModule;
}

class dqm::DQMHitModule : public art::EDAnalyzer {
public:
  explicit DQMHitModule(fhicl::ParameterSet const & p);
  virtual ~DQMHitModule();

  void analyze(art::Event const & e) override;

  void beginJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;

private:
  using HitWireAssns_t = art::Assns<recob::Hit, recob::Wire>;

  void createAssocVector(
    HitWireAssns_t const&,
    std::vector< std::vector<int> > &
    );

  void createAssocVector(std::vector<recob::Wire> const&,
			 std::vector<recob::Hit> const&,
			 std::vector< std::vector<int> > &);

  void createMCAssocVector( std::vector<recob::Wire> const&,
			    std::vector<sim::MCHitCollection> const&,
			    std::vector< std::vector<int> >&);

  // Declare member data here.
  std::vector<std::string> fHitModuleLabels;
  std::string fMCHitModuleLabel;
  std::string fWireModuleLabel;

  //TTree *wireDataTree;
  //std::vector<TTree *> hitDataTree;

  TH1F *NHits_Plane0;
  TH1F *NHits_Plane1;
  TH1F *NHits_Plane2;

  TH1F *HitTimeMean_Plane0;
  TH1F *HitTimeMean_Plane1;
  TH1F *HitTimeMean_Plane2;

  TH1F *HitTimeVar_Plane0;
  TH1F *HitTimeVar_Plane1;
  TH1F *HitTimeVar_Plane2;

  TH1F *HitTimeSkew_Plane0;
  TH1F *HitTimeSkew_Plane1;
  TH1F *HitTimeSkew_Plane2;
  
  TH1F *HitChargeMean_Plane0;
  TH1F *HitChargeMean_Plane1;
  TH1F *HitChargeMean_Plane2;

  TH1F *HitChargeVar_Plane0;
  TH1F *HitChargeVar_Plane1;
  TH1F *HitChargeVar_Plane2;

  TH1F *HitChargeSkew_Plane0;
  TH1F *HitChargeSkew_Plane1;
  TH1F *HitChargeSkew_Plane2;


  TH1F *NFlashes;

  TH1F *MeanFlashLight;

  DQMHitAlg analysisAlg;

};


dqm::DQMHitModule::DQMHitModule(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{ this->reconfigure(p); }

dqm::DQMHitModule::~DQMHitModule()
{
  // Clean up dynamic memory and other resources here.
}

void dqm::DQMHitModule::createAssocVector(
  HitWireAssns_t const& HitToWire,
  std::vector< std::vector<int> > & WireHitAssocVector
) {
  // WireHitAssocVector: for each wire, indices of all the hits associated to it
  
  // the iteration to art::Assns<Hit, Wire> points to a art::Ptr pair (assn_t)
  // with a hit as first element ("left") and a wire as the second one ("right")
  for (HitWireAssns_t::assn_t const& assn: HitToWire)
    WireHitAssocVector.at(assn.second.key()).push_back(assn.first.key());

}

void dqm::DQMHitModule::createAssocVector(std::vector<recob::Wire> const& wireVector,
					  std::vector<recob::Hit> const& hitVector,
					  std::vector< std::vector<int> > & WireHitAssocVector)
{
  WireHitAssocVector.resize(wireVector.size());
  for(size_t iwire=0; iwire<wireVector.size(); iwire++){
    for(size_t ihit=0; ihit<hitVector.size(); ihit++){
      if(hitVector[ihit].Channel()==wireVector[iwire].Channel())
	WireHitAssocVector[iwire].push_back(ihit);
    }
  }
}

void dqm::DQMHitModule::createMCAssocVector( std::vector<recob::Wire> const& wireVector,
					     std::vector<sim::MCHitCollection> const& mcHitVector,
					     std::vector< std::vector<int> > & WireMCHitAssocVector){
  
  WireMCHitAssocVector.clear(); WireMCHitAssocVector.resize(wireVector.size());
  
  //first, store all the MCHitCollection indices in a map keyed on channel
  //then, loop through wires, and lookup mchitcollections based on the wire's channel
  
  std::map<unsigned int,std::vector<int> > mcHitIndicesByChannel;
  for(unsigned int icol=0; icol<mcHitVector.size(); icol++)
    mcHitIndicesByChannel[mcHitVector[icol].Channel()].push_back(icol);
  
  
  for(unsigned int iwire=0; iwire<wireVector.size(); iwire++)
    WireMCHitAssocVector[iwire].insert(WireMCHitAssocVector[iwire].end(),
				       mcHitIndicesByChannel[wireVector[iwire].Channel()].begin(),
				       mcHitIndicesByChannel[wireVector[iwire].Channel()].end());
  
}

void dqm::DQMHitModule::analyze(art::Event const & e)
{
  //get event and run numbers
  unsigned int eventNumber = e.id().event();
  unsigned int runNumber   = e.run();

  analysisAlg.ClearHitModules();

  //get time service
  const detinfo::DetectorClocks* ts = lar::providerFrom<detinfo::DetectorClocksService>();

  //get the wire data
  art::Handle< std::vector<recob::Wire> > wireHandle;
  e.getByLabel(fWireModuleLabel,wireHandle);
  std::vector<recob::Wire> const& wireVector(*wireHandle);

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
 

  //get the MC hits
  /*
  art::Handle< std::vector<sim::MCHitCollection> > mcHitHandle;
  e.getByLabel(fMCHitModuleLabel,mcHitHandle);
  std::vector<sim::MCHitCollection> const& mcHitVector(*mcHitHandle);
  */

  //make the association vector. First index is wire index, second is mcHitCollection index
  std::vector< std::vector<int> > WireMCHitAssocVector;
  //createMCAssocVector(wireVector,mcHitVector,WireMCHitAssocVector);

  //get the hit data
  size_t nHitModules = fHitModuleLabels.size();
  std::vector< art::Handle< std::vector<recob::Hit> > > hitHandles(nHitModules);
  // for each hit module output (first index), for each wire (second index)
  // the list of hits associated with that wire is stored
  std::vector< std::vector< std::vector<int> > > WireHitAssocVectors(nHitModules);
  for(size_t iter=0; iter < nHitModules; iter++){

    e.getByLabel(fHitModuleLabels[iter],hitHandles[iter]);

    //create association vectors by hand for now
    //art::ValidHandle<HitWireAssns_t> HitToWireAssns
    //= e.getValidHandle<HitWireAssns_t>(fHitModuleLabels[iter]);
    //WireHitAssocVectors[iter].resize(wireVector.size());
    //createAssocVector(*HitToWireAssns,WireHitAssocVectors[iter]);

    WireHitAssocVectors[iter].resize(wireVector.size());
    art::Handle<HitWireAssns_t> HitToWireAssns;
    if( e.getByLabel(fHitModuleLabels[iter],HitToWireAssns) )
    {
      createAssocVector(*HitToWireAssns,WireHitAssocVectors[iter]);
    }
    else
    {
      createAssocVector(wireVector,*(hitHandles[iter]),WireHitAssocVectors[iter]);
    }
    
    //load in this hit/assoc pair
    analysisAlg.LoadHitAssocPair( *(hitHandles[iter]),
				  WireHitAssocVectors[iter],
				  fHitModuleLabels[iter]);
    
  }

  //run the analzyer alg
  analysisAlg.AnalyzeWires(wireVector,/*mcHitVector,WireMCHitAssocVector,*/ts,eventNumber,runNumber, NHits_Plane0, NHits_Plane1, NHits_Plane2, HitTimeMean_Plane0, HitTimeMean_Plane1, HitTimeMean_Plane2, HitTimeVar_Plane0,  HitTimeVar_Plane1, HitTimeVar_Plane2, HitTimeSkew_Plane0, HitTimeSkew_Plane1, HitTimeSkew_Plane2, HitChargeMean_Plane0, HitChargeMean_Plane1, HitChargeMean_Plane2, HitChargeVar_Plane0, HitChargeVar_Plane1, HitChargeVar_Plane2, HitChargeSkew_Plane0, HitChargeSkew_Plane1, HitChargeSkew_Plane2);
  analysisAlg.AnalyzeFlashes(flashVector, NFlashes, MeanFlashLight);
}

void dqm::DQMHitModule::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  //wireDataTree = tfs->make<TTree>("wireDataTree","WireDataTree");
  //analysisAlg.SetWireDataTree(wireDataTree);

  NHits_Plane0 = tfs->make<TH1F>("NHits_Plane0", "# hits in plane 0", 150, 0, 30000);
  NHits_Plane1 = tfs->make<TH1F>("NHits_Plane1", "# hits in plane 1", 150, 0, 30000);
  NHits_Plane2 = tfs->make<TH1F>("NHits_Plane2", "# hits in plane 2", 150, 0, 30000);

  HitTimeMean_Plane0 = tfs->make<TH1F>("HitTimeMean_Plane0", "Mean hit time in plane 0"       , 300, -6000, 6000);
  HitTimeMean_Plane1 = tfs->make<TH1F>("HitTimeMean_Plane1", "Mean hit time in plane 1"       , 300, -6000, 6000);
  HitTimeMean_Plane2 = tfs->make<TH1F>("HitTimeMean_Plane2", "Mean hit time in plane 2"       , 300, -6000, 6000);
  HitTimeVar_Plane0  = tfs->make<TH1F>("HitTimeVar_Plane0" , "Variance on hit time in plane 0", 100, 0, 1E5);
  HitTimeVar_Plane1  = tfs->make<TH1F>("HitTimeVar_Plane1" , "Variance on hit time in plane 1", 100, 0, 1E5);
  HitTimeVar_Plane2  = tfs->make<TH1F>("HitTimeVar_Plane2" , "Variance on hit time in plane 2", 100, 0, 1E5);
  HitTimeSkew_Plane0 = tfs->make<TH1F>("HitTimeSkew_Plane0", "Skewness on hit time in plane 0", 100, -100, 100);
  HitTimeSkew_Plane1 = tfs->make<TH1F>("HitTimeSkew_Plane1", "Skewness on hit time in plane 1", 100, -100, 100);
  HitTimeSkew_Plane2 = tfs->make<TH1F>("HitTimeSkew_Plane2", "Skewness on hit time in plane 2", 100, -100, 100);

  HitChargeMean_Plane0 = tfs->make<TH1F>("HitChargeMean_Plane0", "Mean hit charge in plane 0"       , 100, 0, 100);
  HitChargeMean_Plane1 = tfs->make<TH1F>("HitChargeMean_Plane1", "Mean hit charge in plane 1"       , 100, 0, 100);
  HitChargeMean_Plane2 = tfs->make<TH1F>("HitChargeMean_Plane2", "Mean hit charge in plane 2"       , 100, 0, 100);
  HitChargeVar_Plane0  = tfs->make<TH1F>("HitChargeVar_Plane0" , "Variance on hit charge in plane 0", 100, 0, 1000);
  HitChargeVar_Plane1  = tfs->make<TH1F>("HitChargeVar_Plane1" , "Variance on hit charge in plane 1", 100, 0, 1000);
  HitChargeVar_Plane2  = tfs->make<TH1F>("HitChargeVar_Plane2" , "Variance on hit charge in plane 2", 100, 0, 1000);
  HitChargeSkew_Plane0 = tfs->make<TH1F>("HitChargeSkew_Plane0", "Skewness on hit charge in plane 0", 100, -100, 100);
  HitChargeSkew_Plane1 = tfs->make<TH1F>("HitChargeSkew_Plane1", "Skewness on hit charge in plane 1", 100, -100, 100);
  HitChargeSkew_Plane2 = tfs->make<TH1F>("HitChargeSkew_Plane2", "Skewness on hit charge in plane 2", 100, -100, 100);

  NFlashes       = tfs->make<TH1F>("NFlashes"      , "# of flashes"        , 300, 0, 300);
  MeanFlashLight = tfs->make<TH1F>("MeanFlashLight", "Mean light per flash", 50, 0, 50);
 
  // The below creates a Tree with one branch - a recob::Hit branch - for each 
  // Hit module specified in the fcl file. So, don't run this module once per Hit
  // Finder as one normally would. Just run it once and specify all HitFinders.
  // This was the design for the WireDataTree; we follow it here for the 
  // Hit trees.
  
  /*
  for(auto const& label : fHitModuleLabels) {
    std::string firstArg("hitData_");
    firstArg += label;
    std::string secArg("HitDataTree_");
    secArg += label;
    TTree* intermediateTree = tfs->make<TTree>(firstArg.c_str(),secArg.c_str());
    hitDataTree.push_back(intermediateTree);
  }
  analysisAlg.SetHitDataTree(hitDataTree);*/
}

void dqm::DQMHitModule::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.

  fHitModuleLabels  = p.get< std::vector<std::string> >("HitModuleLabels");
  fWireModuleLabel  = p.get< std::string              >("WireModuleLabel");
  fMCHitModuleLabel = p.get< std::string              >("MCHitModuleLabel");

}

DEFINE_ART_MODULE(dqm::DQMHitModule)
