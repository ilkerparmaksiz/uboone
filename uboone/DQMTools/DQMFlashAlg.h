#ifndef DQMFLASHALG_H
#define DQMFLASHALG_H

#include <vector>
#include <string>
#include <exception>

#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "lardataobj/RecoBase/OpFlash.h"

#include "TTree.h"
#include "TH1F.h"

namespace dqm{

  class DQMFlashAlgException : public std::exception{
    virtual const char* what() const throw(){
      return "DQMFlashAlg Exception";
    }
  } dqmflashalgexception;

  class DQMFlashAlg{

  public:

    DQMFlashAlg();
    
    void AnalyzeFlashes(std::vector<recob::OpFlash> flashVector,
        TH1F *NFlashes, TH1F *MeanFlashLight, TH1F *FlashLightVar, TH1F *FlashLightSkew);

  private:
    
    //this is for unit testing...class has no other purpose
    friend class DQMFlashAlgTest;

  };

}//end namespace hit


#endif
