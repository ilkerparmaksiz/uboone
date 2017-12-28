////////////////////////////////////////////////////////////////////////
// Class:       CoherentPionAna
// Module Type: analyzer
// File:        CoherentPionAna_module.cc
//
// Generated at Tue May 16 17:21:02 2017 by Ilker Parmaksiz using artmod
// from cetpkgsupport v1_11_00.
////////////////////////////////////////////////////////////////////////
#include <map>
#include <string>
#include <algorithm>
#include <iostream>
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
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"

/////////////////////////Geometry//////////////////////
#include "larcore/Geometry/GeometryCore.h"


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
#include "TAxis.h"

class CoherentPionAna;

// Class for each Interaction
struct Interaction_Info{
  
  //General Variables For MCTruth
  int Expected_InteractionType;
  int GENI_InteractionType;
  int BothPandMOutTPC=0;
  int BothPandMInTPC=0;
  Double_t pdgCode;
  
  TLorentzVector Nu_P;
  TLorentzVector G4_Abst;
  TLorentzVector G4_QSquared;
  
  std::string FileName;
  std::string Title;
  std::ofstream file;

  std::vector <double> GENI_Particle;
  std::vector<std::string> G4_Variable={"Vx","Vy","Vz","EndX","EndY","EndZ","Px","Py","Pz","E","Mass"};
  
  std::map<int,std::string> Expected_Particles;
  std::map<int,int> ParticleCount;
  std::map<int,int> FinalStateParticles;
  std::map <int,int> OnlyOneInTPC={{13,0},{-13,0},{211,0},{-211,0},{2212,0},{-2212,0},{2112,0}};
  std::map<int,std::pair<int,int>> Pair;
  
  
  //for Histograms
  std::vector <TH1D*> h1D;
  std::vector <const char *> h1D_name;
  std::vector <const char *> h1D_title;
  //std::vector <const char *> h1D_Xaxis={"Px  [GeV/c]","Py  [GeV/c]","Pz  [GeV/c]","E  [GeV/c]","Px  [GeV/c]","Py  [GeV/c]","Pz  [GeV/c]","E  [GeV/c]","Q^2 [(GeV/c)^2]","|t| [(GeV/c)^2]","Cone Angle [Degrees]","Angle [Degrees]","DOCA [cm]"};
  //std::vector<std::vector<int>> h1D_binds={{250,-5,5},{250,-5,5},{250,-5,5},{250,-5,5},{250,-5,5},{250,-5,5},{250,-5,5},{250,-5,5},{250,0,3},{250,0,3},{720,-360,360},{720,-360,360},{250,0,100}};
  std::vector <const char *> h1D_Xaxis={"Q^2 [(GeV/c)^2]","|t| [(GeV/c)^2]","Cone Angle [Degrees]","Angle [Degrees]","DOCA [cm]","VertexEnergy[MeV]"};
  std::vector<std::vector<int>> h1D_binds={{250,0,5},{250,0,5},{720,-360,360},{720,-360,360},{1000,0,1000},{10000,0,10000}}; 
  
  
  // Variables For MCTrack
  std::map <unsigned int,std::vector <Double_t>> MCTrack_Particles;
  std::map <unsigned int,std::vector<Double_t>> MCT_X;
  std::map <unsigned int,std::vector<Double_t>> MCT_Y;
  std::map <unsigned int,std::vector<Double_t>> MCT_Z;  
  int MCT_TrackID;
  double MCT_pdgCode;
  double MCT_Distance;
  std::string MCT_Process;
  int Mu_TrackID=-1;
  int Pi_TrackID=-1;
  unsigned int VertexTotal=0;  
  
  //for efficiency 
  int Total_Cone=0;
  int Total_DOCA=0;
  
  // Vertex activity
  int NPlanes=3;
  int XPlChs;
  int YPlChs;
  int ZPlChs;//collection plane chanels
  unsigned int NearCh;
  std::vector<double> Vertex;
  std::vector<double>OldVertex={-1,-1,-1};
  
  

  // This function does the printing in the end of each job
  // I specificly use this "<>"  so external programs do the counting form me
  void Print()
  {   
    //Prints How many particles in the TPC or out of TPC and so on
   /* std::cout<<"-------"<< Title <<"-------"<<std::endl;
    for (auto k=Expected_Particles.begin();k!=Expected_Particles.end();++k)
    { 
      for (auto i=OnlyOneInTPC.begin();i!=OnlyOneInTPC.end();++i)
      { 
        if(k->first==i->first)
            std::cout<<Title<<"_Only " <<k->second<< " in TPC<>" << i->second <<std::endl;
   ii   }
    }
  std::cout<<Title<<"_Both  Particles Out of TPC<> "<< BothPandMOutTPC <<std::endl;
  std::cout<<Title<<"_Both  Particles in TPC<> "<< BothPandMInTPC <<std::endl;
     i
  for (auto i=FinalStateParticles.begin();i!=FinalStateParticles.end();++i)
    std::cout<<Title<<"_"<<i->first<< "<>" << i->second <<std::endl;
  for(auto i=ParticleCount.begin();i!=ParticleCount.end();++i)
    std::cout<<Title<<"_Extra_"<<i->first<< "<>" << i->second <<std::endl;
    */
    
  // printing efficency info       
  std::cout<<"A_"<<Title<<"TotalCone<>"<<Total_Cone<<std::endl;
  std::cout<<"B_"<<Title<<"TotalDOCA<>"<<Total_DOCA<<std::endl;
  std::cout<<"C_"<<Title<<"VertexTotal<>"<<VertexTotal<<std::endl;
    
  }  
  
  //it counts how many times an integer presents in a map 
  void Count(std::map<int,int> &Particles,int pdg)
  {
    auto sr=Particles.find(pdg);
    if(sr!=Particles.end())
      sr->second+=1;      
  }
  
  // it checks if a certain pdgcode exists in a map or not
  bool DoExist(std::map<int,int> particle,int PdgCode)
  {
    auto sr=particle.find(PdgCode);
    if(sr!=particle.end())
      return true;
    return false;
  }
  
  // it accepts pdgcode and checks if it is expected particle(not in use currently)
  int IsExpectedParticle(int PdgCode)
  {
    auto sr=Expected_Particles.find(PdgCode);
    if(sr!=Expected_Particles.end())
      return PdgCode;
    return -1;
  }
  
  // This one checks any element present in any type of vector such as int to int or double to double(not in use currently)
  template <typename T> 
  const bool Check(std::vector<T>& Vec, const T& Element) 
  {
    if (std::find(Vec.begin(), Vec.end(), Element) != Vec.end())
        return true;

    return false;
  }
  
  //To open a file to store (must be called at begining of the job)(not in use currently)
  void FileOpen(){
    file.open(FileName,std::ofstream::out | std::ofstream::app);
   }  
  //Header for the file 
  void Header()
  {
    if(Expected_InteractionType==GENI_InteractionType)
    {
      file <<"Interaction Type<>"<<GENI_InteractionType<<"<>";      
      file <<"Neutrino<>"<<Nu_P.Px()<<"<>"<<Nu_P.Py()<<"<>"<<Nu_P.Pz()<<"<>"<<Nu_P.Pt()<<"\n";
    }
  }
  
  // it writes the provided information to file
  void WriteToFile(std::vector<double> G4_Particle,std::vector<std::string>G4_Variable,int G4_pdgCode,TLorentzVector G4_Momentum)
  { 
    file << "G4_pdgCode<>"<<G4_pdgCode<<"<>";
    for(unsigned int i=0;i<G4_Variable.size();++i)
      file <<G4_Variable[i]<<"<>"<<G4_Particle[i]<<" <> ";
    file << "G4_Momentum<>"<<G4_Momentum.Px()<<"<>"<<G4_Momentum.Py()<<"<>"<<G4_Momentum.Pz()<<"<>"<<G4_Momentum.Pt()<<"\n";
  }
 
  // close the file (must be at the end of the job)(not in use currently)
  void FileClose(){
    file.close();
  }
 
   
 

  // This function does all the main job for two particles for each interaction
  void Histogram_Fill(std::multimap<int,std::vector<double>> &G4_Particles,std::map<int,std::vector<double>> &Energy)
  { 
   
   if(Expected_InteractionType==GENI_InteractionType && Expected_InteractionType!=0)
   {
      
      //some needed variables 
      std::map<int,int> InTPC;
      std::map<int,int> OutTPC;
      
      TLorentzVector PMuPi[4];
      
      unsigned int mucount=1;
      unsigned int pioncount=1;
      
      int first;
      int second;
      int Particle;
      
      
    
      auto ip=Pair.find(Expected_InteractionType);
      for(auto i=G4_Particles.begin();i!=G4_Particles.end();++i)
      {
        if (ip!=Pair.end())
        {
        
          first=std::get<0>(ip->second);
          second=std::get<1>(ip->second);
            
          if(i->first==first && first!=0) //for pion
          {		
            if(G4_Particles.count(second)!=0)
            {
              Particle=second;
              Pi_TrackID=(int)(i->second)[11];
            }
            else
              Particle=0; 
        
          
            if(G4_Particles.count(i->first)<=G4_Particles.count(Particle) && Particle!=0)
            {
               // Histograms of momentum and Energy for pion
              /*for(unsigned int k=0;k<4;++k)
                  h1D[k]->Fill((i->second)[k+6]);*/
              
              PMuPi[1].SetPxPyPzE((i->second)[6],(i->second)[7],(i->second)[8],(i->second)[9]);
              InorOut(FinalStateParticles,i->first,i->second,InTPC,OutTPC);

            }else if(G4_Particles.count(i->first)>G4_Particles.count(Particle) && Particle!=0)
            {
              auto it =Energy.find(i->first);
           
              if(it!=Energy.end())
              {
                if((i->second[9]==*std::max_element(it->second.begin(),it->second.end()) && pioncount<=G4_Particles.count(Particle)) ) 
                {               
                    // Histograms of momentum and Energy for pion if there are more than one
                    /*for(unsigned int k=0;k<4;++k)
                        h1D[k]->Fill((i->second)[k+6]);*/
                  
                    PMuPi[1].SetPxPyPzE((i->second)[6],(i->second)[7],(i->second)[8],(i->second)[9]);
                    InorOut(FinalStateParticles,i->first,i->second,InTPC,OutTPC);
                    pioncount+=1;

                    it->second.erase(std::max_element(it->second.begin(),it->second.end()));

                }  
                else
                {  
                    auto it2=ParticleCount.find(first);
                    if(it2!=ParticleCount.end())
                    {
                        it2->second+=1;
                        
                        //WriteToFile(i->second,G4_Variable,i->first,PMuPi[1]);
                    
                    }else
                        ParticleCount.insert(std::pair<int,int>(i->first,1));

                }

              }
            }
           
          }
        
          else if(i->first==second && second!=0) //for muon
          {
            if(G4_Particles.count(first)!=0)
            {
              Particle=first;
              Mu_TrackID=(int)(i->second)[11];
            }
            else
              Particle=0;
            
            
            if(G4_Particles.count(i->first)<=G4_Particles.count(Particle) && Particle!=0)
            {
                // Histogram of momentum and energy for muon
              /*for(unsigned int k=0;k<4;++k)
	              h1D[k+4]->Fill((i->second)[k+6]); */
		          
		          PMuPi[0].SetPxPyPzE((i->second)[6],(i->second)[7],(i->second)[8],(i->second)[9]);
		          InorOut(FinalStateParticles,i->first,i->second,InTPC,OutTPC);
		        }
		        else if(G4_Particles.count(i->first)>G4_Particles.count(Particle) && Particle!=0 )
		        {
              auto it =Energy.find(i->first);
              if(it!=Energy.end())
              {
                if(i->second[9]==*std::max_element(it->second.begin(),it->second.end()) && mucount<=G4_Particles.count(Particle))
                { 
                // Histogram of momentum and energy for muon if there are more than one
		            /* for(unsigned int k=0;k<4;++k)
	                    h1D[k+4]->Fill((i->second)[k+6]);*/    
	                    
		              PMuPi[0].SetPxPyPzE((i->second)[6],(i->second)[7],(i->second)[8],(i->second)[9]);
		              InorOut(FinalStateParticles,i->first,i->second,InTPC,OutTPC);
		              mucount+=1;               

		              it->second.erase(std::max_element(it->second.begin(),it->second.end()));
                }

                else
                {  
                    auto it2=ParticleCount.find(first);
                    if(it2!=ParticleCount.end())
                    {
                        it2->second+=1;
                        
                        //WriteToFile(i->second,G4_Variable,i->first,PMuPi[1]); //if info needed to be written to file
                        
                    }else
                        ParticleCount.insert(std::pair<int,int>(i->first,1)); 
                }

              }	      
		        }
          }
        
        }
         
		    //Histogram for openning angle,DOCA,ConeAngle
       if((PMuPi[1].Mag()!=0 && PMuPi[0].Mag()!=0) && (PMuPi[0]!=PMuPi[2] || PMuPi[1]!=PMuPi[3]))
        {

          double G4_Abs_t=0;
          double G4_Q_Squared=0;
          //PMuPi[0].Print();
          //PMuPi[1].Print();
          
          //Cone 
          double newAngle,newAnglelength;
          double Xcomp,Ycomp,Zcomp;
          double DOCA_Out=9999;
        
          double MudotPi=(PMuPi[0].Px()*PMuPi[1].Px() + PMuPi[0].Py()*PMuPi[1].Py() + PMuPi[0].Pz()*PMuPi[1].Pz());
          double length=(sqrt(pow(PMuPi[0].Px(),2) + pow(PMuPi[0].Py(),2) + pow(PMuPi[0].Pz(),2))*(sqrt(pow(PMuPi[1].Px(),2) + pow(PMuPi[1].Py(),2) + pow(PMuPi[1].Pz(),2))));
                              
          double MudotPiDeg=acos(MudotPi/length)*(180/3.14);
	        if(Mu_TrackID!=-1 && Pi_TrackID!=-1)
	        {
	          DOCA_Out=Cal_DOCA(Mu_TrackID,Pi_TrackID);
	          if(DOCA_Out!=9999)
	          {
              h1D[4]->Fill(DOCA_Out);
              Total_DOCA=Total_DOCA+1;
	          }
	        }
            
         
          
          // -|t|-
          G4_Abst=Nu_P-PMuPi[0]-PMuPi[1];
          G4_Abs_t=fabs(G4_Abst.Dot(G4_Abst));
        
          // Q 
          G4_QSquared=Nu_P-PMuPi[0];
          G4_Q_Squared=fabs(G4_QSquared.Dot(G4_QSquared));
          
          Xcomp=PMuPi[0].Px() + PMuPi[1].Px();
          Ycomp=PMuPi[0].Py() + PMuPi[1].Py();
          Zcomp=PMuPi[0].Pz() + PMuPi[1].Pz();
          
          newAnglelength=(sqrt(pow(Xcomp,2) + pow(Ycomp,2) + pow(Zcomp,2)));
          newAngle=(acos(Zcomp/newAnglelength)* (180/3.14));
          
          std::vector <double> combo ={G4_Q_Squared,G4_Abs_t,newAngle,MudotPiDeg};
          
          for (unsigned int i=0;i<combo.size();i++)
            h1D[i]->Fill(combo[i]);
          
          Total_Cone=Total_Cone+1;
          
          PMuPi[0]=PMuPi[2];
          PMuPi[1]=PMuPi[3];
          
          //VertexAct();
          
          // Set TrackID to -1
          Mu_TrackID=-1;
          Pi_TrackID=-1;
         
        }
        
     }
      
     if(OutTPC.size()==0 && InTPC.size()>0)
	    BothPandMInTPC++;	                      // CC Coh-Pion events where both the muon and the pion are fully contained inside the TPC	   
	   else if (InTPC.size()==0 && OutTPC.size()>0)
	    BothPandMOutTPC++;                      //CC Coh-Pion Events where both the muon and the pion exit the TPC
  	 else
  	 {
  	  for(auto p=InTPC.begin(); p!=InTPC.end();++p)
  	   { 
  	      for (auto i=OnlyOneInTPC.begin();i!=OnlyOneInTPC.end();++i)
          {
            if(p->first==i->first)
              i->second+=p->second;   
          }
       }
     }   
   }
   
 }
 
  // calculates how many particles in TPC or out of TPC
  void InorOut(std::map<int,int> &FinalStateParticles,int i,std::vector <double> b,std::map<int,int> &InTPC,std::map<int,int> &OutTPC)
  {
    if(!DoExist(FinalStateParticles,i))
        FinalStateParticles.insert(std::pair<int,int>((i),1));
    else
        Count(FinalStateParticles,(i));
      
      // In TPC
   
      if((b)[3]>0 && (b)[3]<250 && (b)[4]>-125 && (b)[4]<125 && (b)[5]>0 && (b)[5]<1030)
       {
          if(!DoExist(InTPC,(i)))
            InTPC.insert(std::pair<int,int>((i),1) );
          else
            Count(InTPC,(i));
       }
      else // Out of TPC
       {
          if(!DoExist(OutTPC,(i)))
            OutTPC.insert(std::pair<int,int>((i),1));
          else
            Count(OutTPC,(i));
       }
       
   }
  
  //setting Histograms
  void Histogram_Set()
  {
    art::ServiceHandle<art::TFileService> tfs;
    for(unsigned int k=0;k<h1D_name.size();++k)
    {
      h1D.push_back(tfs->make<TH1D>(h1D_name[k],h1D_title[k],h1D_binds[k][0],h1D_binds[k][1],h1D_binds[k][2]));
      h1D[k]->GetYaxis()->SetTitle("Number of Events");
      h1D[k]->GetXaxis()->SetTitle(h1D_Xaxis[k]);
    }
  }
  
  // info from MCTrack to calculate distance of approach
  void DOCA(std::vector<art::Ptr<sim::SimChannel>>&simList,art::Handle< std::vector<sim::MCTrack> > &mctrackh)
  {
  
  
   if(Expected_InteractionType==GENI_InteractionType && Expected_InteractionType!=0)
   {
  
          MCT_X.erase(MCT_X.begin(),MCT_X.end());
          MCT_Y.erase(MCT_Y.begin(),MCT_Y.end());
          MCT_Z.erase(MCT_Z.begin(),MCT_Z.end());
          Vertex={GENI_Particle[0],GENI_Particle[1],GENI_Particle[2]};
          MCTrack_Particles.erase(MCTrack_Particles.begin(),MCTrack_Particles.end());
          double Etotal=0;
    //DOCA 
    for(std::vector<sim::MCTrack>::const_iterator imctrk = mctrackh->begin();imctrk != mctrackh->end(); ++imctrk) 
     {
      
        const sim::MCTrack& mctrk = *imctrk;
        MCT_pdgCode=mctrk.PdgCode();
        MCT_TrackID=mctrk.TrackID();
        //std::cout<<"MCTRACKID="<<MCT_TrackID<<std::endl;
        //std::cout<<"PDGCode="<<MCT_pdgCode<<std::endl;
        MCT_Process=mctrk.Process();           
     
        MCT_Distance=sqrt((mctrk.End().X()-mctrk.Start().X())*(mctrk.End().X()-mctrk.Start().X()) + (mctrk.End().Y()-mctrk.Start().Y())*(mctrk
     .End().Y()-mctrk.Start().Y()) + (mctrk.End().Z()-mctrk.Start().Z())*(mctrk.End().Z()-mctrk.Start().Z()));
     
      std::vector<Double_t> data={MCT_pdgCode,mctrk.Start().X(),mctrk.Start().Y(),mctrk.Start().Z(),mctrk.End().X(),mctrk.End().Y(),mctrk.End().Z(),mctrk.Start().Px(),mctrk.Start().Py(),mctrk.Start().Pz(),mctrk.Start().E(),MCT_Distance};
     
                   
        MCTrack_Particles.insert(std::pair<unsigned int,std::vector<Double_t>>(MCT_TrackID,data));
        
        std::vector<Double_t> XStep;
        std::vector<Double_t> YStep;
        std::vector<Double_t> ZStep;
       // std::cout<<"----PDGCode="<<MCT_pdgCode<<" ------"<<std::endl;
        //std::cout<<"MCTRACKID="<<MCT_TrackID<<std::endl;
        for(auto step:mctrk)
        { 
           /*std::cout<<"MCSTEPX="<<step.X()<<std::endl;
           std::cout<<"Y="<<step.Y()<<std::endl;
           std::cout<<"Z="<<step.Z()<<std::endl;
          */
                
           XStep.push_back(step.X());
           YStep.push_back(step.Y());
           ZStep.push_back(step.Z());
           //std::cout<<"MCSTEPE="<<step.E()<<std::endl;
          if((step.X()>=0 && step.X()<=250) && (step.Y()>=-125 && step.Y()<=125) && (step.Z()>=0 && step.Z()<=1030)){
           if((MCT_pdgCode!=111 || MCT_pdgCode!=2112))
           {
             std::vector<double> XYZ={step.X(),step.Y(),step.Z()};
             if(abs(XYZ[0]-Vertex[0])<5 && abs(XYZ[1]-Vertex[1])< 5 && abs(XYZ[2]-Vertex[2])<5)
               VertexAct(simList,XYZ,Etotal);
           }
          }
        }
        
        
        if(XStep.empty() && YStep.empty() && ZStep.empty())
        { 
          if((Vertex[0]>=0 && Vertex[0]<=250) && (Vertex[1]>=-125 && Vertex[1]<=125) && (Vertex[2]>=0 && Vertex[2]<=1030)){  
            if((MCT_pdgCode!=111 || MCT_pdgCode!=2112))
              VertexAct(simList,Vertex,Etotal);
           }
        
            XStep={data[1],data[4]};
            YStep={data[2],data[5]};
            ZStep={data[3],data[6]};
        }         


        MCT_X.insert(std::pair<unsigned int,std::vector<Double_t>>(MCT_TrackID,XStep));
        MCT_Y.insert(std::pair<unsigned int,std::vector<Double_t>>(MCT_TrackID,YStep));
        MCT_Z.insert(std::pair<unsigned int,std::vector<Double_t>>(MCT_TrackID,ZStep));
     }
     //std::cout<<"Etotal="<<Etotal<<std::endl;
       h1D[5]->Fill(Etotal);
      VertexTotal=VertexTotal+1;     
    }
  }
  
  // Distance of Closest approach calculation
  double Cal_DOCA(int MCT_TrackID1,int MCT_TrackID2)
  {
      
    unsigned int sizet;
    auto x1=MCT_X.find(MCT_TrackID1);
    auto y1=MCT_Y.find(MCT_TrackID1);
    auto z1=MCT_Z.find(MCT_TrackID1);
            
    auto x2=MCT_X.find(MCT_TrackID2);
    auto y2=MCT_Y.find(MCT_TrackID2);
    auto z2=MCT_Z.find(MCT_TrackID2);
                
    if(x1!=MCT_X.end() && x2!=MCT_X.end() && y1!=MCT_Y.end() && y2!=MCT_Y.end() && z1!=MCT_Z.end() && z2!=MCT_Z.end())
    {   
       (x1->second.size()<x2->second.size()) ? sizet=x1->second.size():sizet=x2->second.size();

       
       Double_t closest=9999;
                                                                          
       for(unsigned int t=sizet;t>0;--t)
       {
         if(x1->second[t]!=x2->second[t] && y1->second[t]!=y2->second[t] && z1->second[t]!=z2->second[t])
         {
            Double_t Result=sqrt(pow(x1->second[t]-x2->second[t],2)+pow(y1->second[t]-y2->second[t],2) + pow(z1->second[t]-z2->second[t],2));
            //std::cout<<"Result="<<Result<<std::endl;
            if(Result<closest)
              closest=Result;       
         }
       }
      //std::cout<<"closest="<<closest<<std::endl;   
      return closest;     
     }
     
     return 9999;
  }
    
    
   void VertexAct(std::vector<art::Ptr<sim::SimChannel>>&simList,std::vector<double> & XYZ,double & Etotal)
   {
   if(Expected_InteractionType==GENI_InteractionType && Expected_InteractionType!=0)
   {
    // === Geometry Service ===
    art::ServiceHandle<geo::Geometry>geo;
    
    
    //unsigned int TrackID1;
     //TrackID2=-1;
      XPlChs=geo->Nwires(0);
      YPlChs=geo->Nwires(1);    
      ZPlChs=geo->Nwires(2);
           
    for(size_t nChan = 0; nChan < simList.size(); nChan++) 
    {
       
     if(simList.at(nChan)->Channel()==geo->NearestChannel(XYZ,2))
     {
       // Get the information of each wire
       const auto & wire = simList.at(nChan)->TDCIDEMap();
       for(auto it=wire.begin();it!=wire.end();it++)
       {
          // Looping over the IDEs in a given time tick
          for(size_t i=0;i<it->second.size();i++)
          {
            //Only non-showering (nonnegative) primary track IDEs
            //if(it->second.at(i).trackID > 13 ) { continue; } 
            
            //TrackID1=it->second.at(i).trackID;
           // std::cout<<TrackID1<<std::endl;
            //std::cout<<"TrackID"<<it->second.at(i).trackID<<std::endl;
            //std::cout<<"Eng: "<<it->second.at(i).energy<<std::endl;
            //std::cout<<"X: "<<it->second.at(i).x<<std::endl;
            //std::cout<<"Y: "<<it->second.at(i).y<<std::endl;
            //std::cout<<"Z: "<<it->second.at(i).z<<std::endl;
            Etotal+=it->second.at(i).energy;
            
           // if(TrackID1!=TrackID2 && i>0)
            
          }
          
         }
         //std::cout<<Etotal<<std::endl;
       }
     }
   }
 
 }
};

// Declaring Interactions

Interaction_Info cccpion;
Interaction_Info ccQE;
Interaction_Info ccres;
Interaction_Info NCRes;
Interaction_Info NCDis;
Interaction_Info Counter;

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
  std::string fSimChanModuleLabel;
  

  int TotalNuInteractions=0;
  int TotalNuInteractionsinVertex=0;
  int CCNuInteractions=0;
  int NCNuInteractions=0;
  //int CCpionNuInteractions=0;
  int TotalEvents=0;
  int CCRes=0;
  int ncres=0;
  
  
  
  // Reference-> http://nusoft.fnal.gov/larsoft/doxsvn/html/namespacesimb.html#a2cce734d1b71408bbc7d98d148ac4360a501129a20d2a93f092d7c8c9ae63fd66
 
  // total Interactions involved
  std::map<int,int> Interactions={{1001,0},{1003,0},{1004,0},{1005,0},{1010,0},{1011,0},{1012,0},{1017,0},{1021,0},{1028,0},{1032,0},{1039,0},{1041,0},{1046,0},{1048,0},{1053,0},{1055,0},{1060,0},{1062,0},{1067,0},{1070,0},{1073,0},{1076,0},{1079,0},{1080,0},{1085,0},{1086,0},{1090,0},{1097,0},{1092,0},{1007,0},{1009,0},{1014,0},{1016,0}};
    
};


CoherentPionAna::CoherentPionAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  this->reconfigure(p);
  
}




// Easy String Print
void SWrite(std::string Variable)
{
 std::cout<<Variable<<std::endl;
}


void CoherentPionAna::analyze(art::Event const & evt)
{

  // MC truth information
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  
  
  
  //McTrack information  	
  art::Handle< std::vector<sim::MCTrack> > mctrackh;
  evt.getByLabel(fMCTrackModuleLabel, mctrackh);
  
  // === BackTracker service ===
  art::ServiceHandle<cheat::BackTracker> bt;
  const sim::ParticleList& plist = bt->ParticleList();
  
  // for G4 particle double check	
  std::vector<const simb::MCParticle* > geant_part;
  for(size_t p = 0; p < plist.size(); ++p) 
    geant_part.push_back(plist.Particle(p));  
  
  if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
    art::fill_ptr_vector(mclist, mctruthListHandle);
    
   //For SimChannel  
    art::Handle< std::vector<sim::SimChannel> > SimListHandle; 
    std::vector<art::Ptr<sim::SimChannel> > simList;    
    
    if(evt.getByLabel(fSimChanModuleLabel, SimListHandle))
    { art::fill_ptr_vector(simList, SimListHandle); }
  
  
  //NChan=geo->NearestChannel();
  // looping to get MCTruth particle information
  for(unsigned int iList = 0; (iList < mclist.size()); ++iList)
  {
    
    // For GENIE information
    std::vector <double>  GENI_Particle={mclist[iList]->GetNeutrino().Nu().Vx(),mclist[iList]->GetNeutrino().Nu().Vy(),mclist[iList]->GetNeutrino().Nu().Vz()};
    
    std::string GENI_Variable[]={"Vx","Vy","Vz"};
    std::map<int,std::vector<double>>Energy;
    
    TLorentzVector Nu_P=mclist[iList]->GetNeutrino().Nu().Momentum();
    
    int GENI_pdgCode=mclist[iList]->GetNeutrino().Nu().PdgCode();
    int CCNC=mclist[iList]->GetNeutrino().CCNC();
    int InteractionType=mclist[iList]->GetNeutrino().InteractionType();
     
    // look for muNu and +MuNu Interactions
    if(GENI_pdgCode==14 ||  GENI_pdgCode==-14)
    { 
      // Total Number of Neutrino interactions
      TotalNuInteractions++;
      
      // Total Number of Neutrino interactions inside the vertex
      if(GENI_Particle[0]>=0 && GENI_Particle[0]<= 250  && GENI_Particle[1]>=-125  && GENI_Particle[1]<= 125  && GENI_Particle[2]>=0 && GENI_Particle[2]<=1030 )
      {
	      TotalNuInteractionsinVertex++;
   
        if(CCNC==0)
        {
	        CCNuInteractions++;
	        
	        //CCCoherent Pion Interaction
	        cccpion.Expected_InteractionType=1097;
	        cccpion.Title="CC-COH";
	        cccpion.Expected_Particles={{211,"pi+"},{13,"mu-"}};
	        cccpion.Pair={{1097,{211,13}}};
	       
	        // Queasi Eleastic Interaction          
	        ccQE.Expected_InteractionType=1001;
	        ccQE.Title="CCQE";
	        ccQE.Expected_Particles={{13,"mu-"},{2212,"p+"}};
	        ccQE.Pair={{1001,{2212,13}}};
	        

  
          // CCres 
          if((InteractionType>=1003 && InteractionType<=1090) && (InteractionType!=1097 || InteractionType!=1001 || InteractionType!=1007 || InteractionType!=1009 || InteractionType!=1014 || InteractionType!=1016 ))
          {    
            ccres.Expected_InteractionType=InteractionType;
            auto it =ccres.Pair.find(InteractionType);
            if(it==ccres.Pair.end())
              ccres.Pair.insert(std::pair<int,std::pair<int,int>>(InteractionType,{211,13}));
          }
	        ccres.Title="CCRes";
	        ccres.Expected_Particles={{13,"mu-"},{211,"pi+"}};
	           
	        
	      }
	      else if(CCNC==1)
	      {
	        NCNuInteractions++;
	        
	        //NCDis
	        NCDis.Expected_InteractionType=1092;
	        NCDis.Title="NCDis";
	        NCDis.Expected_Particles={{211,"pi+"},{-211,"pi-"}};
	        NCDis.Pair={{1092,{211,-211}}};
	           
	         //NCRes  
	        if(InteractionType==1007 || InteractionType==1009 || InteractionType==1014 || InteractionType==1016)
             NCRes.Expected_InteractionType=InteractionType;
                
	        NCRes.Title="NCRes";
	        NCRes.Expected_Particles={{211,"pi+"},{-211,"pi-"}};
	        NCRes.Pair={{1007,{211,-211}},{1009,{211,-211}},{1014,{211,-211}},{1016,{211,-211}}};
	        
	      }    
	      
        if((InteractionType>=1003 && InteractionType<=1090) || InteractionType==1001 || InteractionType==1092 || InteractionType==1097)
	       Counter.Count(Interactions,InteractionType);
	     	      
        std::multimap<int,std::vector <double>> G4_Particles;
        
        
        ccres.GENI_InteractionType=NCDis.GENI_InteractionType=NCRes.GENI_InteractionType=InteractionType;
        ccres.Nu_P=NCDis.Nu_P=NCRes.Nu_P=Nu_P;
        NCDis.GENI_Particle=NCRes.GENI_Particle=ccres.GENI_Particle=GENI_Particle;
           
        cccpion.Nu_P=ccres.Nu_P=ccQE.Nu_P=Nu_P;
	      cccpion.GENI_InteractionType=ccres.GENI_InteractionType=ccQE.GENI_InteractionType=InteractionType;
	      cccpion.GENI_Particle=ccres.GENI_Particle=ccQE.GENI_Particle=GENI_Particle;
      
        //looping over giant info
        for(unsigned int i = 0; (i < geant_part.size()); ++i )
         { 
           // For G4 information
           int G4_pdgCode=geant_part[i]->PdgCode();
           double TrackID=geant_part[i]->TrackId();
   
           
           std::vector <double> G4_Particle={geant_part[i]->Vx(),geant_part[i]->Vy(),geant_part[i]->Vz(),geant_part[i]->EndX(),geant_part[i]->EndY(),geant_part[i]->EndZ(),geant_part[i]->Px(),geant_part[i]->Py(),geant_part[i]->Pz(),geant_part[i]->E(),geant_part[i]->Mass(),TrackID};
           TLorentzVector G4_Momentum=geant_part[i]->Momentum();	
                       
	         if((float)GENI_Particle[0]==(float)G4_Particle[0] && (float)GENI_Particle[1]==(float)G4_Particle[1] && (float)GENI_Particle[2]==(float)G4_Particle[2] && (float)G4_Particle[0]!=(float)G4_Particle[3])
           {
                
             G4_Particles.insert(std::pair<int,std::vector<double>>(G4_pdgCode,G4_Particle));
                
             auto it=Energy.find(G4_pdgCode);
             if(it!=Energy.end())
                it->second.push_back(G4_Particle[9]);
             else
                Energy.insert(std::pair<int,std::vector<double>>(G4_pdgCode,{G4_Particle[9]}));

             //WriteToFile(G4_Particle,G4_Variable,G4_pdgCode,G4_Momentum);
             // Collecting Final State Particles
            } 
    
          }
          // Distance of closest approach
          cccpion.DOCA(simList,mctrackh);
	        ccQE.DOCA(simList,mctrackh);
	        ccres.DOCA(simList,mctrackh);
	        NCDis.DOCA(simList,mctrackh);
	        NCRes.DOCA(simList,mctrackh);
	        
	        
	        //Vertex Activity
	       /* cccpion.VertexAct(simList);
	        ccQE.VertexAct(simList);
	        ccres.VertexAct(simList);
	        NCDis.VertexAct(simList);
	        NCRes.VertexAct(simList);
          */
	        //Histogram info
          cccpion.Histogram_Fill(G4_Particles,Energy);
          ccQE.Histogram_Fill(G4_Particles,Energy);
          ccres.Histogram_Fill(G4_Particles,Energy);
	        NCRes.Histogram_Fill(G4_Particles,Energy);
          NCDis.Histogram_Fill(G4_Particles,Energy);
	  
	  
      }
    }
  }
}


void CoherentPionAna::beginJob()
{ 
  
 // Setting Histogram variables
  cccpion.h1D_name={"Q_Squared_ccpi","Abs_t_ccpi","ConeAngle_ccpi","MudotPi_ccpi","DOCA_ccpi","VertexAct_ccpi"};
  cccpion.h1D_title={"Q_Squared for CC Coherent Pion Interaction","|t| for CC Coherent Pion Interaction","Cone Angle(ccpi)","Angle Between Muon's and Pion's Momentum","DOCA of ccpi","Vertex Activity for cccpion"};
  ccQE.h1D_name={"Q_Squared_ccqe","Abs_t_ccqe","ConeAngle_ccqe","MudotPi_ccqe","DOCA_ccqe","VertexAct_ccqe"};
  ccQE.h1D_title={"Q_Squared for CC QE Type 1 Pion Interaction","|t| for CC QE Type 1 Interaction","Cone Angle(ccqe)","Angle Between Muon's and Pion's Momentum","DOCA(ccqe)","Vertex Activity for CCQE"};
  

  ccres.h1D_name={"Q_Squared_ccres","Abs_t_ccres","ConeAngle_ccres","MudotPi_ccres","DOCA_ccres","Vertex_ccres"};
  ccres.h1D_title={"Q_Squared for CC Ress Pion Interaction","|t| for CC Ress  Interaction","ConeAngle(ccres)","Angle Between Muon's and Pion's Momentum","DOCA(ccres)","Vertex Activity for CCRes"}; 
	
  NCRes.h1D_name={"Q_Squared_ncres","Abs_t_ncres","ConeAngle_ncres","MudotPi_ncres","DOCA_ncres","VertexAct_ncres"};
  NCRes.h1D_title={"Q_Squared for NCRes","|t| for NCRes","ConeAngle(ncres)","Angle Between Pi(+) and Pi(-)","DOCA(NCRes)","Vertex Act for NCRes"};
 
  NCDis.h1D_name={"Q_Squared_dis","Abs_tdis","ConeAngle_ncdis","MudotPi_dis","DOCA_dis","VertexAct_NCDis"};
  NCDis.h1D_title={"Q_Squared for NCDis","|t| for NCDis","ConeAngle(ncdis)","Angle Between Pi(+) and Pi(-)","DOCA(dis)","Vertex Act for NCDsi"};
	
  
  
  // Histogram 
 
  
  
  /////////////////////////////////////////////
  /////// #### DEFINING HISTOGRAMS #### ///////
  /////////////////////////////////////////////
  
  cccpion.Histogram_Set();
  ccQE.Histogram_Set();
  ccres.Histogram_Set();
  NCRes.Histogram_Set();
  NCDis.Histogram_Set();
 
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
  // printing in the end of job
  std::cout<<"Total Events<> " << TotalEvents <<std::endl;
  std::cout<<"TotalNuInteractions<> " << TotalNuInteractions <<std::endl;
  std::cout<<"TotalNuInteractionsinVertex<> " << TotalNuInteractionsinVertex <<std::endl;
  std::cout<<"CCNuInteractions<>" << CCNuInteractions <<std::endl;
  std::cout<<"NCNuInteractions<>" << NCNuInteractions <<std::endl;
 

  
  // if CCpionNuInteractions Exists
  if(CCNuInteractions>0)
  {
    cccpion.Print();
    ccQE.Print();
    ccres.Print();
    NCRes.Print();
    NCDis.Print();
    for(auto i=Interactions.begin();i!=Interactions.end();++i)
    { 
      
      if(i->first==1001)
        std::cout<<"CCQECount<>"<< i->second<<std::endl;
      if(i->first==1097)
        std::cout<<"CC-COH<>"<< i->second<<std::endl;
      if(i->first>=1003 && i->first<=1090 && (i->first!=1007 || i->first!=1009 || i->first!=1014 || i->first!=1016))
        CCRes+=i->second;
      if(i->first==1007 || i->first==1009 || i->first==1014 || i->first==1016 )
        ncres+=i->second;
      if(i->first==1092)
        std::cout<<"NCDISCount<>"<< i->second<<std::endl;
      
      
    }
    std::cout<<"CCResCount_<>"<< CCRes<<std::endl;
    std::cout<<"NCResCount_<>"<< ncres<<std::endl;

  }
  
 
  
}

void CoherentPionAna::endRun(art::Run const & r)
{

}
void CoherentPionAna::endSubRun(art::SubRun const & sr)
{
  // Implementation of optional member function here.

/*
  art::Handle< sumdata::POTSummary > potListHandle;
  if(sr.getByLabel(fPOTModuleLabel,potListHandle))
    pot=potListHandle->totpot;
      
  //POT->Fill(pot);
  TotPOT.push_back(pot);
  */
}

void CoherentPionAna::reconfigure(fhicl::ParameterSet const & p)
{
 fGenieGenModuleLabel = p.get < std::string >("GenieGenModuleLabel");
 fG4ModuleLabel = p.get < std::string >("G4ModuleLabel");
 fLArGeantModuleLabel = p.get < std::string >("LArGeantModuleLabel");
 fMCShowerModuleLabel= p.get < std::string >("MCShowerModuleLabel");
 fMCTrackModuleLabel= p.get < std::string >("MCTrackModuleLabel");
 fPOTModuleLabel =p.get < std::string >("POTModuleLabel");
 fSimChanModuleLabel =p.get < std::string >("SimChanModuleLabel");
  
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
