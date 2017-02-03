#include "DQMHitAlg.h"

#include <functional>
#include <unordered_map>

dqm::DQMHitAlg::DQMHitAlg() {
  wireData.NHitModules = 0;
}

void dqm::DQMHitAlg::SetWireDataTree(TTree *wdt){
  wireDataTree = wdt;
  SetupWireDataTree();
}

void dqm::DQMHitAlg::SetHitDataTree(std::vector<TTree*>& trees){

  hitDataTree.clear();
  hitData.clear();

  hitDataTree.reserve (trees.size());
  // This is particularly important: to establish the hitData memory -- all before 
  // individually making the tree->Branch() calls, specifying those addresses.
  hitData.reserve     (trees.size());

  // Construct the local attribute data container
  for(auto const& t : trees) {
    hitDataTree.push_back(t);
    hitData.push_back(new recob::Hit);
  }

  for(size_t i=0; i<hitData.size(); ++i)
    hitDataTree[i]->Branch(hitDataTree[i]->GetName(),"recob::Hit",&(hitData[i]));
}

void dqm::DQMHitAlg::SetupWireDataTree(){
  wireDataTree->Branch("event", &wireData.event, "event/i");
  wireDataTree->Branch("run", &wireData.run, "run/i");
  wireDataTree->Branch("channel", &wireData.channel, "channel/i");
  wireDataTree->Branch("plane", &wireData.plane, "plane/i");
  wireDataTree->Branch("roi_index", &wireData.range_index, "roi_index/i");
  wireDataTree->Branch("roi_start", &wireData.range_start, "roi_start/i");
  wireDataTree->Branch("roi_size", &wireData.range_size, "roi_size/i");
  wireDataTree->Branch("roi_charge", &wireData.integrated_charge, "roi_charge/F");
  wireDataTree->Branch("roi_peak_charge", &wireData.peak_charge, "roi_peak_charge/F");
  wireDataTree->Branch("roi_peak_time", &wireData.peak_time, "roi_peak_time/F");
  wireDataTree->Branch("nHitModules", &wireData.NHitModules, "nHitModules/I");
  wireDataTree->Branch("HitModuleLabels",&wireData.HitModuleLabels);
  wireDataTree->Branch("NHits",&wireData.NHits);
  wireDataTree->Branch("Hits_IntegratedCharge",&wireData.Hits_IntegratedCharge);
  wireDataTree->Branch("Hits_AverageCharge",&wireData.Hits_AverageCharge);
  wireDataTree->Branch("Hits_PeakCharge",&wireData.Hits_PeakCharge);
  wireDataTree->Branch("Hits_Peak",&wireData.Hits_PeakTime);
  wireDataTree->Branch("Hits_wAverageCharge",&wireData.Hits_wAverageCharge);
  wireDataTree->Branch("Hits_wAverageTime",&wireData.Hits_wAverageTime);
  wireDataTree->Branch("Hits_MeanMultiplicity",&wireData.Hits_MeanMultiplicity);

  //---------------------------LEAVE BLANK------------------------------------------
  wireDataTree->Branch("NMCHits",&wireData.NMCHits);
  wireDataTree->Branch("MCHits_IntegratedCharge",&wireData.MCHits_IntegratedCharge);
  wireDataTree->Branch("MCHits_AverageCharge",&wireData.MCHits_AverageCharge);
  wireDataTree->Branch("MCHits_PeakCharge",&wireData.MCHits_PeakCharge);
  wireDataTree->Branch("MCHits_Peak",&wireData.MCHits_PeakTime);
  wireDataTree->Branch("MCHits_wAverageCharge",&wireData.MCHits_wAverageCharge);
  wireDataTree->Branch("MCHits_wAverageTime",&wireData.MCHits_wAverageTime);
  //--------------------------------------------------------------------------------

}


void dqm::DQMHitAlg::ClearHitModules(){
  HitModuleLabels.clear();
  HitProcessingQueue.clear();
  wireData.NHitModules = 0;
}

void dqm::DQMHitAlg::LoadHitAssocPair( std::vector<recob::Hit> const& HitVector,
				       std::vector< std::vector<int> > const& AssocVector,
				       std::string const& HitModuleLabel){
  
  HitProcessingQueue.push_back( std::make_pair( std::cref(HitVector), std::cref(AssocVector)) );
  HitModuleLabels.push_back(HitModuleLabel);

  if(HitProcessingQueue.size()!=HitModuleLabels.size())
    throw dqmhitalgexception;

}

void dqm::DQMHitAlg::AnalyzeWires(std::vector<recob::Wire> const& WireVector,
				  const detinfo::DetectorClocks *ts,
				  unsigned int event, unsigned int run,
      TH1F *nhits_plane0, TH1F *nhits_plane1, TH1F *nhits_plane2,
      TH1F *meant_plane0, TH1F *meant_plane1, TH1F *meant_plane2,
      TH1F *vart_plane0 , TH1F *vart_plane1 , TH1F *vart_plane2 ,
      TH1F *skewt_plane0, TH1F *skewt_plane1, TH1F *skewt_plane2,
      TH1F *meanq_plane0, TH1F *meanq_plane1, TH1F *meanq_plane2,
      TH1F *varq_plane0 , TH1F *varq_plane1 , TH1F *varq_plane2 ,
      TH1F *skewq_plane0, TH1F *skewq_plane1, TH1F *skewq_plane2){

  int plane0_hitno = 0; 
  int plane1_hitno = 0; 
  int plane2_hitno = 0; 

  TH1F *plane0_time = new TH1F ("plane0_time", "plane0_time", 300, -6000, 6000);
  TH1F *plane1_time = new TH1F ("plane1_time", "plane1_time", 300, -6000, 6000);
  TH1F *plane2_time = new TH1F ("plane2_time", "plane2_time", 300, -6000, 6000);

  TH1F *plane0_charge = new TH1F ("plane0_charge", "plane0_charge", 100, 0, 100);
  TH1F *plane1_charge = new TH1F ("plane1_charge", "plane1_charge", 100, 0, 100);
  TH1F *plane2_charge = new TH1F ("plane2_charge", "plane2_charge", 100, 0, 100);

  InitWireData(event,run);
  for(size_t iwire=0 ; iwire < WireVector.size(); iwire++)
  {
    if (WireVector[iwire].View() == 0)
      FillWireInfo(WireVector[iwire], iwire, ts, plane0_hitno, plane0_charge, plane0_time);
    if (WireVector[iwire].View() == 1)
      FillWireInfo(WireVector[iwire], iwire, ts, plane1_hitno, plane1_charge, plane1_time);
    if (WireVector[iwire].View() == 2)
      FillWireInfo(WireVector[iwire], iwire, ts, plane2_hitno, plane2_charge, plane2_time);
  }

  nhits_plane0->Fill(plane0_hitno);
  nhits_plane1->Fill(plane1_hitno);
  nhits_plane2->Fill(plane2_hitno);

  meant_plane0->Fill(plane0_time->GetMean());
  meant_plane1->Fill(plane1_time->GetMean());
  meant_plane2->Fill(plane2_time->GetMean());
  vart_plane0->Fill((plane0_time->GetRMS())*(plane0_charge->GetRMS()));
  vart_plane1->Fill((plane1_time->GetRMS())*(plane1_charge->GetRMS()));
  vart_plane2->Fill((plane2_time->GetRMS())*(plane2_charge->GetRMS()));
  skewt_plane0->Fill(plane0_time->GetSkewness());
  skewt_plane1->Fill(plane1_time->GetSkewness());
  skewt_plane2->Fill(plane2_time->GetSkewness());

  meanq_plane0->Fill(plane0_charge->GetMean());
  meanq_plane1->Fill(plane1_charge->GetMean());
  meanq_plane2->Fill(plane2_charge->GetMean());
  varq_plane0->Fill((plane0_charge->GetRMS())*(plane0_charge->GetRMS()));
  varq_plane1->Fill((plane1_charge->GetRMS())*(plane1_charge->GetRMS()));
  varq_plane2->Fill((plane2_charge->GetRMS())*(plane2_charge->GetRMS()));
  skewq_plane0->Fill(plane0_charge->GetSkewness());
  skewq_plane1->Fill(plane1_charge->GetSkewness());
  skewq_plane2->Fill(plane2_charge->GetSkewness());

}

void dqm::DQMHitAlg::AnalyzeFlashes(std::vector<recob::OpFlash> flashVector, TH1F *NFlashes, TH1F *MeanFlashLight)
{
  std::cout << "number of flashes: " << flashVector.size() << std::endl;
  NFlashes->Fill(flashVector.size());
  double totflashlight = 0;
  int num_pmts = 32;
  for (unsigned int i = 0; i < flashVector.size(); i++)
  {
    for (int pmt_index = 0; pmt_index < num_pmts; pmt_index++)
    {
      totflashlight += flashVector[i].PE(pmt_index);
    }
  }
  double meanflashlight = totflashlight/(double(flashVector.size()*num_pmts));
  MeanFlashLight->Fill(meanflashlight);
}

void dqm::DQMHitAlg::InitWireData(unsigned int event, unsigned int run){

  wireData.event = event;
  wireData.run = run;
  wireData.NHitModules = HitModuleLabels.size();
  wireData.HitModuleLabels = HitModuleLabels;

}

void dqm::DQMHitAlg::ClearWireDataHitInfo(){
  wireData.NMCHits = 0;
  wireData.MCHits_IntegratedCharge = 0;
  wireData.MCHits_AverageCharge = 0;
  wireData.MCHits_PeakCharge = -999;
  wireData.MCHits_PeakTime = 0;
  wireData.MCHits_wAverageCharge = 0;
  wireData.MCHits_wAverageTime = 0;
  
  wireData.NHits.assign(wireData.NHitModules,0);
  wireData.Hits_IntegratedCharge.assign(wireData.NHitModules,0);
  wireData.Hits_AverageCharge.assign(wireData.NHitModules,0);
  wireData.Hits_PeakCharge.assign(wireData.NHitModules,-999);
  wireData.Hits_PeakTime.assign(wireData.NHitModules,0);
  wireData.Hits_wAverageCharge.assign(wireData.NHitModules,0);
  wireData.Hits_wAverageTime.assign(wireData.NHitModules,0);
  wireData.Hits_MeanMultiplicity.assign(wireData.NHitModules,0);
  wireData.Hits.clear(); wireData.Hits.resize(wireData.NHitModules);
}

void dqm::DQMHitAlg::FillWireInfo(recob::Wire const& wire, 
				  int WireIndex,
				  const detinfo::DetectorClocks *ts,
      int &hitno,
      TH1F* charge_hist,
      TH1F* time_hist){
  
  wireData.channel = wire.Channel();
  wireData.plane = wire.View();
  unsigned int range_index = 0;

  for( auto const& range : wire.SignalROI().get_ranges() ){

    wireData.range_index = range_index;
    wireData.range_start = range.begin_index();
    wireData.range_size = range.size();

    ClearWireDataHitInfo();

    ProcessROI(range, WireIndex, ts, hitno, charge_hist, time_hist);
    range_index++;

  }//end loop over roi ranges

}

void dqm::DQMHitAlg::ROIInfo(lar::sparse_vector<float>::datarange_t const& range,
			     float& charge_sum,
			     float& charge_peak,
			     float& charge_peak_time){

  charge_sum=0;
  charge_peak = -999;
  unsigned int counter=range.begin_index();

  for(auto const& value : range){
    charge_sum += value;
    if(value > charge_peak){
      charge_peak = value;
      charge_peak_time = (float)counter;
    }
    counter++;
  }

}

void dqm::DQMHitAlg::ProcessROI(lar::sparse_vector<float>::datarange_t const& range, 
				int WireIndex,
				const detinfo::DetectorClocks *ts,
    int &hitno,
    TH1F* charge_hist,
    TH1F* time_hist){

  ROIInfo(range,wireData.integrated_charge,wireData.peak_charge,wireData.peak_time);

  for(size_t iter = 0; iter < HitProcessingQueue.size(); iter++)
    FindAndStoreHitsInRange(HitProcessingQueue[iter].first, 
			    HitProcessingQueue[iter].second.at(WireIndex),
			    iter,
			    range.begin_index(),
			    range.begin_index()+range.size(),
       hitno,
       charge_hist,
       time_hist);

}

void dqm::DQMHitAlg::FindAndStoreHitsInRange( std::vector<recob::Hit> const& HitVector,
					      std::vector<int> const& HitsOnWire,
					      size_t hitmodule_iter,
					      size_t begin_wire_tdc,
					      size_t end_wire_tdc,
           int &hitno,
           TH1F *charge_hist,
           TH1F *time_hist){

  for( auto const& hit_index : HitsOnWire){
    recob::Hit const& thishit = HitVector.at(hit_index);

    //check if this hit is on this ROI
    if( thishit.PeakTimeMinusRMS() < begin_wire_tdc ||
	thishit.PeakTimePlusRMS() > end_wire_tdc)
      continue;

    hitno++;
    charge_hist->Fill(thishit.Integral());
    time_hist->Fill(thishit.PeakTime());

  }

}

void dqm::DQMHitAlg::FillHitInfo(recob::Hit const& hit, std::vector<HitInfo>& HitInfoVector){
  HitInfoVector.emplace_back(hit.PeakTime(),
			     hit.SigmaPeakTime(),
			     hit.RMS(),
			     hit.StartTick(),
			     hit.EndTick(),
			     hit.Integral(),
			     hit.SigmaIntegral(),
			     hit.PeakAmplitude(),
			     hit.SigmaPeakAmplitude(),
			     hit.GoodnessOfFit());
}
