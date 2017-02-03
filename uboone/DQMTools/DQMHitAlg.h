#ifndef DQMHITALG_H
#define DQMHITALG_H

#include <vector>
#include <string>
#include <exception>

//#include "lardata/MCBase/MCHitCollection.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "lardataobj/RecoBase/OpFlash.h"

#include "TTree.h"
#include "TH1F.h"

namespace dqm{

  struct HitInfo{
    
    //need a constructor here
    HitInfo(float pt, float pt_s,
	    float w,
	    int st, int et,
	    float c, float c_s,
	    float mc, float mc_s,
	    float gof)
      : peaktime(pt)
      , peaktime_sigma(pt_s)
      , rms(w)
      , starttick(st)
      , endtick(et)
      , charge(c)
      , charge_sigma(c_s)
      , maxcharge(mc)
      , maxcharge_sigma(mc_s)
      , goodness_of_fit(gof)
    {}

    float peaktime;
    float peaktime_sigma;
    float rms;
    int starttick;
    int endtick;
    float charge;
    float charge_sigma;
    float maxcharge;
    float maxcharge_sigma;
    float goodness_of_fit;
  };

  struct WireROIInfo{
    unsigned int event;
    unsigned int run;
    unsigned int channel;
    unsigned int plane;
    unsigned int range_index;
    unsigned int range_start;
    size_t range_size;
    float integrated_charge;
    float peak_charge;
    float peak_time;
    int NHitModules;
    std::vector<std::string> HitModuleLabels;
    std::vector<int> NHits;
    std::vector<float> Hits_IntegratedCharge;
    std::vector<float> Hits_AverageCharge;
    std::vector<float> Hits_PeakCharge;
    std::vector<float> Hits_PeakTime;
    std::vector<float> Hits_wAverageCharge;
    std::vector<float> Hits_wAverageTime;
    std::vector<float> Hits_MeanMultiplicity;
    std::vector< std::vector<HitInfo> > Hits;
    int NMCHits;
    float MCHits_IntegratedCharge;
    float MCHits_AverageCharge;
    float MCHits_PeakCharge;
    float MCHits_PeakTime;
    float MCHits_wAverageCharge;
    float MCHits_wAverageTime;
  };


  class DQMHitAlgException : public std::exception{
    virtual const char* what() const throw(){
      return "DQMHitAlg Exception";
    }
  } dqmhitalgexception;

  class DQMHitAlg{

    typedef std::pair< const std::vector<recob::Hit>& , const std::vector< std::vector<int> >& > HitAssocPair;

  public:

    DQMHitAlg();
    
    void SetWireDataTree(TTree*);

    void SetHitDataTree(std::vector<TTree*>& trees);

    void AnalyzeWires(std::vector<recob::Wire> const&,
		      //std::vector<sim::MCHitCollection> const&,
		      //std::vector< std::vector<int> > const&,
		      const detinfo::DetectorClocks *,
		      unsigned int,
		      unsigned int,
        TH1F *nhits_plane0, TH1F *nhits_plane1, TH1F *nhits_plane2,
        TH1F *meant_plane0, TH1F *meant_plane1, TH1F *meant_plane2,
        TH1F *vart_plane0 , TH1F *vart_plane1 , TH1F *vart_plane2 ,
        TH1F *skewt_plane0, TH1F *skewt_plane1, TH1F *skewt_plane2,
        TH1F *meanq_plane0, TH1F *meanq_plane1, TH1F *meanq_plane2,
        TH1F *varq_plane0 , TH1F *var1_plane1 , TH1F *varq_plane2,
        TH1F *skewq_plane0, TH1F *skewq_plane1, TH1F *skewq_plane2);

    void AnalyzeFlashes(std::vector<recob::OpFlash> flashVector,
        TH1F *NFlashes, TH1F *MeanFlashLight);

    void LoadHitAssocPair( std::vector<recob::Hit> const&, 
			   std::vector< std::vector<int> > const&,
			   std::string const&);

    void ClearHitModules();
    
  private:
    
    void InitWireData(unsigned int, unsigned int);
    void ClearWireDataHitInfo();

    void FillHitInfo(recob::Hit const&, std::vector<HitInfo>&);

    void FillWireInfo(recob::Wire const&, 
		      int,
		      //std::vector<sim::MCHitCollection> const&,
		      //std::vector<int> const&,
		      const detinfo::DetectorClocks *,
        int &hitno,
        TH1F* charge_hist,
        TH1F *time_hist);

    void ProcessROI(lar::sparse_vector<float>::datarange_t const&, int,
		    //std::vector<sim::MCHitCollection> const&,
		    //std::vector<int> const&,
		    const detinfo::DetectorClocks *,
      int &hitno,
      TH1F* charge_hist,
      TH1F *time_hist);

    void ROIInfo(lar::sparse_vector<float>::datarange_t const&,
		 float&,float&,float&);

    void FindAndStoreHitsInRange(std::vector<recob::Hit> const&,
				 std::vector<int> const&,
				 size_t,size_t,size_t,
     int &hitno,
     TH1F *charge_hist,
     TH1F *time_hist);

    /*
    void FindAndStoreMCHitsInRange(std::vector<sim::MCHitCollection> const&,
				   std::vector<int> const&,
				   size_t,size_t,
				   const detinfo::DetectorClocks *);
       */
    
    WireROIInfo wireData;
    std::vector<recob::Hit*> hitData;

    std::vector<std::string> HitModuleLabels;
    std::vector< HitAssocPair > HitProcessingQueue;


    void SetupWireDataTree();
    TTree* wireDataTree;

    std::vector<TTree*> hitDataTree;

    //this is for unit testing...class has no other purpose
    friend class DQMHitAlgTest;

  };

}//end namespace hit


#endif
