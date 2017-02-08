////////////////////////////////////////////////////////////////////////
// Class:       HVBurstFilter
// Plugin Type: filter (art v2_05_00)
// File:        HVBurstFilter_module.cc
//
// Generated at Mon Feb  6 16:13:23 2017 by Wesley Ketchum using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "lardataobj/RawData/RawDigit.h"

namespace dqm {
  class HVBurstFilter;
}


class dqm::HVBurstFilter : public art::EDFilter {
public:
  explicit HVBurstFilter(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  HVBurstFilter(HVBurstFilter const &) = delete;
  HVBurstFilter(HVBurstFilter &&) = delete;
  HVBurstFilter & operator = (HVBurstFilter const &) = delete;
  HVBurstFilter & operator = (HVBurstFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;
  void reconfigure(fhicl::ParameterSet const& p);

private:

  // Declare member data here.
  art::InputTag  fRawDataInputTag;

};


dqm::HVBurstFilter::HVBurstFilter(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  this->reconfigure(p);
}

void dqm::HVBurstFilter::reconfigure(fhicl::ParameterSet const & p)
{
  fRawDataInputTag = p.get<art::InputTag>("RawDataInputTag");
}

bool dqm::HVBurstFilter::filter(art::Event & e)
{

  /**********Now RawWaveform***********/
  auto const& rawtpc_handle = e.getValidHandle<std::vector<raw::RawDigit>>(fRawDataInputTag);
  auto const& rawtpc_vec(*rawtpc_handle);
  
  long int ADCInTime[9595] = {0};
  bool burstFound    = false;
  int firstFireTick = 0;
  
  for (int iADC = 0; iADC < 9595; iADC++ ) 
    {
      for (auto const& rawdigit : rawtpc_vec){
	ADCInTime[iADC] += (rawdigit.ADC(iADC)-rawdigit.GetPedestal());
      }
      if (ADCInTime[iADC]/8256. > 1500. ) {burstFound = true; firstFireTick = iADC; break;}
    }

    if (burstFound)    std::cout<<"BURST FOUND! Run: "<<e.run()<<" Event  "<<e.event()
				<<" First Fire Tick "<< firstFireTick<< std::endl
				<<" Event time Low  "<< e.time().timeLow()
				<<" Event time High "<< e.time().timeHigh()<<std::endl;  

    return burstFound;

}

DEFINE_ART_MODULE(dqm::HVBurstFilter)
