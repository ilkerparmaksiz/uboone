////////////////////////////////////////////////////////////////////////
// Class:       CoherentPionAna
// Module Type: analyzer
// File:        CoherentPionAna_module.cc
//
// Generated at Tue May 16 17:21:02 2017 by Ilker Parmaksiz using artmod
// from cetpkgsupport v1_11_00.
////////////////////////////////////////////////////////////////////////
#include "math.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
///////////////////////////////////////////////////////////////
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCStep.h"
#include "nusimdata/SimulationBase/MCFlux.h"
//////////////////////////////////////////////////
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindMany.h"
///////////////////////////////////////////////////////
#include "larsim/MCCheater/BackTracker.h"
//////////////////////
///// HISTOGRAM /////
/////////////////////
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 

// #####################
// ### ROOT includes ###
// #####################

#include "TComplex.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH3.h"
#include "TF1.h"
#include "TTree.h"
#include "TTimeStamp.h"



class CoherentPionAna;

class CoherentPionAna : public art::EDAnalyzer {
public:
  explicit CoherentPionAna(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CoherentPionAna(CoherentPionAna const &) = delete;
  CoherentPionAna(CoherentPionAna &&) = delete;
  CoherentPionAna & operator = (CoherentPionAna const &) = delete;
  CoherentPionAna & operator = (CoherentPionAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void beginRun(art::Run const & r) override;
  void beginSubRun(art::SubRun const & sr) override;
  void endJob() override;
  void endRun(art::Run const & r) override;
  void endSubRun(art::SubRun const & sr) override;
  void reconfigure(fhicl::ParameterSet const & p) override;
  void respondToCloseInputFile(art::FileBlock const & fb) override;
  void respondToCloseOutputFiles(art::FileBlock const & fb) override;
  void respondToOpenInputFile(art::FileBlock const & fb) override;
  void respondToOpenOutputFiles(art::FileBlock const & fb) override;

private:

  // Declare member data here.
 
  std::string fGenieGenModuleLabel;
  std::string fG4ModuleLabel;
  std::string fLArGeantModuleLabel;
  std::string fMCShowerModuleLabel;
  std::string fMCTrackModuleLabel;
  std::string fPOTModuleLabel;
  
  int Muon=13;
  int Muon_Plus=-13;
  int Pion=211;
  int Pion_Minus=-211;
  int TotalNuInteractions=0;
  int TotalNuInteractionsinVertex=0;
  int CCNuInteractions=0;
  int CCpionNuInteractions=0;
  int Kcoh=0;
  int KCohElastic=0;
  int KCCCOH=0;
  int TotalEvents=0;
  
  //Final State Particles
  std::vector<int> FinalStateParticles;
  std::vector<int> Particle={Pion,Pion_Minus,Muon,Muon_Plus};
  int BothPandMInTPC=0;
  int BothPandMOutTPC=0;
  int OnlyPionInTPC=0;
  int OnlyMuonInTPC=0;
  int InCount=0;
  int OutCount=0;

  
  //Histogram variables
  //TH1F *fNuVtxPosition1D; // 1D Histogram of the x-position of the neutrino vetex
  //TH2F *fNuVtxPosition2D; // 2D Histogram of the x vs y position of the neutrino vertex
  TH1D *hPiPx1D;
  TH1D *hPiPy1D;
  TH1D *hPiPz1D;
  TH1D *hPiPt1D;
  TH1D *hMuPx1D;
  TH1D *hMuPy1D;
  TH1D *hMuPz1D;
  TH1D *hMuPt1D;
  TH1D *hMuPtdotPiPt1D;
  TH1D *hQ_Squared1D;
  TH1D *hAbs_t1D;
};


CoherentPionAna::CoherentPionAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
this->reconfigure(p);
}
template <typename T> 
const bool Check( std::vector<T>& Vec, const T& Element ) 
{
    if (std::find(Vec.begin(), Vec.end(), Element) != Vec.end())
        return true;

    return false;
}
void CoherentPionAna::analyze(art::Event const & evt)
{
  //Number of Events
  TotalEvents++;
  

  
  // MC truth information
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  
    // === BackTracker service ===
  art::ServiceHandle<cheat::BackTracker> bt;
  const sim::ParticleList& plist = bt->ParticleList();
  // for G4 particle double check	
    std::vector<const simb::MCParticle* > geant_part;
    for(size_t p = 0; p < plist.size(); ++p) 
      geant_part.push_back(plist.Particle(p));  
  
  if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
	art::fill_ptr_vector(mclist, mctruthListHandle);

  for(unsigned int iList = 0; (iList < mclist.size()) ; ++iList)
  {
    int InTPCcount=0;
    int OutTPCcount=0;
    std::vector <int> InTPCpdgCode;
    std::vector <int> OutTPCpdgCode;
    TLorentzVector Muon_P;
    TLorentzVector Pion_P;
    TLorentzVector G4_Abst;
    TLorentzVector G4_QSquared;
    
    // For GENIE information
    int  pdgCode=mclist[iList]->GetNeutrino().Nu().PdgCode();
    float  Vx=mclist[iList]->GetNeutrino().Nu().Vx();
    float  Vy=mclist[iList]->GetNeutrino().Nu().Vy();
    float  Vz=mclist[iList]->GetNeutrino().Nu().Vz();
    TLorentzVector Nu_P=mclist[iList]->GetNeutrino().Nu().Momentum();
    int CCNC=mclist[iList]->GetNeutrino().CCNC();
    int InteractionType=mclist[iList]->GetNeutrino().InteractionType();
     
    
    if(pdgCode==14 || pdgCode==-14)
    { 
      // Total Number of Neutrino interactions
      TotalNuInteractions++;
      // Total Number of Neutrino interactions inside the vertex
      if(Vx>0 && Vx< 250  && Vy>-125  && Vy< 125  && Vz>0 && Vz<1030 )
      {
	TotalNuInteractionsinVertex++;
	     
        // Total number of Charged current  neutrino interactions
        if(CCNC==0)
        {
	  CCNuInteractions++;
		
           //filltering to Kcoh =3, kCohElastic=4 and kCCCOH =1097
           if(InteractionType==3 || InteractionType ==4 || InteractionType==1097)
	   {
             //counting coherent pion neutrino interactions
	     if(InteractionType==3)
	       Kcoh++;
	     else if(InteractionType==4)
	       KCohElastic++;
	     else if(InteractionType==1097)
	       KCCCOH++; 
              std::cout<<"Vx = "<<Vx<<std::endl;
 	      std::cout<<"Vy = "<<Vy<<std::endl;
	      std::cout<<"Vz = "<<Vz<<std::endl;
	      std::cout<<"inter = "<<InteractionType<<std::endl;
	    
	     
	     CCpionNuInteractions++;
	     for(unsigned int i = 0; (i < geant_part.size()); ++i )
             { 
               // For G4 information
    	       int G4_pdgCode=geant_part[i]->PdgCode();
    	       float G4_Vx=geant_part[i]->Vx();
    	       float G4_Vy=geant_part[i]->Vy();
   	       float G4_Vz=geant_part[i]->Vz();
   	       double G4_Ex=geant_part[i]->EndX();
   	       double G4_Ey=geant_part[i]->EndY();
   	       double G4_Ez=geant_part[i]->EndZ();
   	       double G4_Px=geant_part[i]->Px();
   	       double G4_Py=geant_part[i]->Py();
   	       double G4_Pz=geant_part[i]->Pz();
   	       double G4_Pt=geant_part[i]->Pt();
	       double G4_P=geant_part[i]->P();
	       double G4_E=geant_part[i]->E();
	       double G4_Mass=geant_part[i]->Mass();
	       TLorentzVector G4_Momentum=geant_part[i]->Momentum();
	       double G4_Q_Squared;
	       double G4_Abs_t;
	       
	       
        

               if(G4_pdgCode==Muon || G4_pdgCode==Muon_Plus || G4_pdgCode==Pion || G4_pdgCode==Pion_Minus)
	       {    
     		
		 if(Vx==G4_Vx && Vy==G4_Vy && Vz==G4_Vz)
		 {
          	   std::cout<<"G4_Pdg = "<<G4_pdgCode<<std::endl;
		   std::cout<<"G4_Vx = "<<G4_Vx<<std::endl;
 		   std::cout<<"G4_Vy = "<<G4_Vy<<std::endl;
	  	   std::cout<<"G4_Vz = "<<G4_Vz<<std::endl;
	   	   std::cout<<"G4_Ex = "<<G4_Ex<<std::endl;
	  	   std::cout<<"G4_Ey = "<<G4_Ey<<std::endl;
	 	   std::cout<<"G4_Ez = "<<G4_Ez<<std::endl;
		   std::cout<<"G4_Px = "<<G4_Px<<std::endl;
		   std::cout<<"G4_Py = "<<G4_Py<<std::endl;
		   std::cout<<"G4_Pz = "<<G4_Pz<<std::endl;
		   std::cout<<"G4_Pt = "<<G4_Pt<<std::endl;
		   std::cout<<"G4_P = "<<G4_P<<std::endl;
		   std::cout<<"G4_E = "<<G4_E<<std::endl;
		   std::cout<<"G4_Mass = "<<G4_Mass<<std::endl;

		   
		 
		   

		   // Collecting Final State Particles
	           FinalStateParticles.push_back (G4_pdgCode);
		     
  	     	   // In TPC
	           if(G4_Ex>0 && G4_Ex<250 && G4_Ey>-125 && G4_Ey<125 && G4_Ez>0 && G4_Ez<1030)
	           {
	  	     InTPCpdgCode.push_back (G4_pdgCode);
	             InTPCcount++;
 	           }else // Out of TPC
	           {
		     OutTPCpdgCode.push_back (G4_pdgCode);
	             OutTPCcount++;
	           }
	           
	           if(G4_pdgCode==Muon) //Histogram for muon
	           {
	             hMuPx1D->Fill(G4_Px);
		     hMuPy1D->Fill(G4_Py);
		     hMuPz1D->Fill(G4_Pz);
		     hMuPt1D->Fill(G4_E);
		     Muon_P=G4_Momentum;
	           }else if(G4_pdgCode==Pion)
		   {
		     hPiPx1D->Fill(G4_Px);
		     hPiPy1D->Fill(G4_Py);
		     hPiPz1D->Fill(G4_Pz);
		     hPiPt1D->Fill(G4_E);
		     Pion_P=G4_Momentum;
		   }
		   //Histogram for Dot product of Muon's and Pion's total momentum
		   if(Pion_P.Mag()!=0 && Muon_P.Mag()!=0)
		   {
		     double MudotPi=Pion_P.Dot(Muon_P);
		     double MudotPiDeg=acos(MudotPi)*(180/3.14);   
		     // -|t|-
		     G4_Abst=Nu_P-Muon_P-Pion_P;
		     G4_Abs_t=fabs(G4_Abst.Dot(G4_Abst));
		     // Q
		     G4_QSquared=Nu_P-Muon_P;
		     G4_Q_Squared=fabs(G4_QSquared.Dot(G4_QSquared));
		     // Filling the histogram
		     std::cout<<"G4_Abs_t = "<<G4_Abs_t<<std::endl;
		     std::cout<<"G4_Abs_t = "<<G4_Q_Squared<<std::endl;
		     hAbs_t1D->Fill(G4_Abs_t);
		     hMuPtdotPiPt1D->Fill(MudotPiDeg);  
		     hQ_Squared1D->Fill(G4_Q_Squared);
		     
		   }

	         }
      	       } 
             }
            
            
	     if(OutTPCcount==0 && InTPCcount>0)
	       BothPandMInTPC++;	                      // CC Coh-Pion events where both the muon and the pion are fully contained inside the TPC	   
	     else if (InTPCcount==0 && OutTPCcount>0)
	       BothPandMOutTPC++;                      //CC Coh-Pion Events where both the muon and the pion exit the TPC
	     else
             {
               if(Check(InTPCpdgCode,Muon) || Check(InTPCpdgCode,Muon_Plus))
	         OnlyMuonInTPC++;
               else if(Check(InTPCpdgCode,Pion) || Check(InTPCpdgCode,Pion_Minus))
		 OnlyPionInTPC++;
	     }    
		   
           }    
	 
         }
      }
   }
 }
 
}	   
  
   


void CoherentPionAna::beginJob()
{
  // Implementation of optional member function here.
  
  // Histogram 
  art::ServiceHandle<art::TFileService> tfs;
  
  
  
  /////////////////////////////////////////////
  /////// #### DEFINING HISTOGRAMS #### ///////
  /////////////////////////////////////////////
  
  //fNuVtxPosition1D = tfs->make<TH1F>("Name of the histogram","Title",250,0,250); //(name,title,Num of bins, start x bin, end x bin
  
  //Pion's Momentum
  //Px
  hPiPx1D=tfs->make<TH1D>("Pxpi","Px (Pion)",250,-5,5);
  hPiPx1D->GetYaxis()->SetTitle("Number of Events");
  hPiPx1D->GetXaxis()->SetTitle("Pion's Px Component of momentum [GeV/c]");
  //Py
  hPiPy1D=tfs->make<TH1D>("Pypi","Py (Pion)",250,-5,5);
  hPiPy1D->GetYaxis()->SetTitle("Number of Events");
  hPiPy1D->GetXaxis()->SetTitle("Pion's Py Component of momentum [GeV/c]");
  //Pz
  hPiPz1D=tfs->make<TH1D>("Pzpi","Pz (Pion)",250,-5,5);
  hPiPz1D->GetYaxis()->SetTitle("Number of Events");
  hPiPz1D->GetXaxis()->SetTitle("Pion's Pz Component of momentum [GeV/c]");
  //E
  hPiPt1D=tfs->make<TH1D>("Ptpi","Pt (Pion)",250,-5,5);
  hPiPt1D->GetYaxis()->SetTitle("Number of Events");
  hPiPt1D->GetXaxis()->SetTitle("Pion's E/c Component of momentum [GeV/c]");
  
  //Muon's Momentum
  
  //Px
  hMuPx1D=tfs->make<TH1D>("PxMu","Px (Muon)",250,-5,5);
  hMuPx1D->GetYaxis()->SetTitle("Number of Events");
  hMuPx1D->GetXaxis()->SetTitle("Muon's Px Component of momentum [GeV/c]");
  //Py
  hMuPy1D=tfs->make<TH1D>("PyMu","Py (Muon)",250,-5,5);
  hMuPy1D->GetYaxis()->SetTitle("Number of Events");
  hMuPy1D->GetXaxis()->SetTitle("Muon's Py Component of momentum [GeV/C]");
  //Pz
  hMuPz1D=tfs->make<TH1D>("PzMu","Pz (Muon)",250,-5,5);
  hMuPz1D->GetYaxis()->SetTitle("Number of Events");
  hMuPz1D->GetXaxis()->SetTitle("Muon's Pz Component of momentum [GeV/c]");
  //Pt
  hMuPt1D=tfs->make<TH1D>("PtMu","Pt (Muon)",250,-5,5);
  hMuPt1D->GetYaxis()->SetTitle("Number of Events");
  hMuPt1D->GetXaxis()->SetTitle("Muon's E/c Component of momentum [GeV/c]");
  
  //Mu_P dot Pi_P 
  hMuPtdotPiPt1D=tfs->make<TH1D>("MuPtdotPiPt","Angle from  Muon's Pt dot pion's Pt",360,-180,180);
  hMuPtdotPiPt1D->GetYaxis()->SetTitle("Number of Events");
  hMuPtdotPiPt1D->GetXaxis()->SetTitle("Angle Between Muon and Pion (Degrees)");
  
  //Q^2
  hQ_Squared1D=tfs->make<TH1D>("Q_Squared","Q_Squared for CC Coherent Pion Interaction",250,0,0.5);
  hQ_Squared1D->GetYaxis()->SetTitle("Number of Events");
  hQ_Squared1D->GetXaxis()->SetTitle("Q^2 [(GeV/c)^2]");
  
  //|t|
  hAbs_t1D=tfs->make<TH1D>("Abs_t","|t| for CC Coherent Pion Interaction",250,0,0.5);
  hAbs_t1D->GetYaxis()->SetTitle("Number of Events");
  hAbs_t1D->GetXaxis()->SetTitle("|t| [(GeV/c)^2]");


}

void CoherentPionAna::beginRun(art::Run const & r)
{
  std::cout<<"Run Starts"<<std::endl;
}

void CoherentPionAna::beginSubRun(art::SubRun const & sr)
{
  // Implementation of optional member function here.
}

void CoherentPionAna::endJob()
{
  // Implementation of optional member function here.
  std::cout<<"Total Events = " << TotalEvents <<std::endl;
  std::cout<<"TotalNuInteractions = " << TotalNuInteractions <<std::endl;
  std::cout<<"TotalNuInteractionsinVertex = " << TotalNuInteractionsinVertex <<std::endl;
  std::cout<<"CCNuInteractions = " << CCNuInteractions <<std::endl;
  std::cout<<"Kcoh = " << Kcoh <<std::endl;
  std::cout<<"KCohElastic = " << KCohElastic <<std::endl;
  std::cout<<"KCCCOH = " << KCCCOH <<std::endl;
  std::cout<<"CCpionNuInteractions = " << CCpionNuInteractions <<std::endl;
  
  // if CCpionNuInteractions Exists
  if(CCpionNuInteractions>0)
  {
    // Final State Particles
    std::sort(FinalStateParticles.begin(),FinalStateParticles.end());
    for(unsigned int i=0;(i<Particle.size());i++)
      std::cout<<"PDG Code ="<< Particle[i] << " = "<<std::count(FinalStateParticles.begin(), FinalStateParticles.end(), Particle[i])<<std::endl;
    
  }
   //CC Coh-Pion Events where only the pion is fully contained (the muon exits the TPC)
   std::cout<<"Only Pion in TPC = " << OnlyPionInTPC <<std::endl;
   //CC Coh-Pion Events where only the muon is fully contained (the pion exits the TPC)
   std::cout<<"Only Muon in TPC = " << OnlyMuonInTPC <<std::endl;
   //CC Coh-Pion Events where both the muon and the pion exit the TPC
   std::cout<<"Both Pion and Muon Out of TPC = " << BothPandMOutTPC <<std::endl;
   // CC Coh-Pion events where both the muon and the pion are fully contained inside the TPC
   std::cout<<"Both Pion and Muon in TPC = " << BothPandMInTPC <<std::endl;
}

void CoherentPionAna::endRun(art::Run const & r)
{
  
}
void CoherentPionAna::endSubRun(art::SubRun const & sr)
{
  // Implementation of optional member function here.
}

void CoherentPionAna::reconfigure(fhicl::ParameterSet const & p)
{
 fGenieGenModuleLabel = p.get < std::string >("GenieGenModuleLabel");
 fG4ModuleLabel = p.get < std::string >("G4ModuleLabel");
 fLArGeantModuleLabel = p.get < std::string >("LArGeantModuleLabel");
 fMCShowerModuleLabel= p.get < std::string >("MCShowerModuleLabel");
 fMCTrackModuleLabel= p.get < std::string >("MCTrackModuleLabel");
 fPOTModuleLabel =p.get < std::string >("POTModuleLabel");
  
}

void CoherentPionAna::respondToCloseInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void CoherentPionAna::respondToCloseOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void CoherentPionAna::respondToOpenInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void CoherentPionAna::respondToOpenOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(CoherentPionAna)
