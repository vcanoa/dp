#include <TH1.h>
#include <TH2.h>
#include "THmulf.h"
#include <TFile.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TGraph.h>
#include <vector>
#include <TVector3.h>
#include <TLorentzVector.h>

#include "THmulf.h"

#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "Reconstruction.h"
#include "QA.h"

#include "phool.h"

//Global Includes
#include "getClass.h"
#include "PHGlobal.h"
#include "recoConsts.h"
#include "Fun4AllServer.h"
#include "Fun4AllReturnCodes.h"
#include "ReactionPlaneObject.h"
#include "ReactionPlaneSngl.h"
#include "RpConst.h"

//Trigger/Header Includes
#include "utiCentrality.h"
#include "TrigLvl1.h"
#include "TriggerHelper.h"
#include "RunHeader.h"

//Track Includes
#include "PHSnglCentralTrack.h"
#include "PHCentralTrack.h"
#include "PHCentralTrackv24.h"
#include "PHSnglCentralTrackv24.h"

//EMCal includes
#include "emcClusterContainer.h"
#include "emcClusterContent.h"

typedef PHIODataNode <PHObject>   PHObjectNode_t;
typedef PHIODataNode <PHCentralTrack>   PHCentralTrackNode_t;

#include "Math/Polynomial.h"
#include "Math/Interpolator.h"

static const double Me = 0.000510998918;  // Particle Data Group 7-25-2004
static const double Me2 = Me*Me;

static const unsigned int BBCLL1_narr = 0x00000004;




QA::QA(const char *outfile) : 
  SubsysReco("QA"), centReco(0)
            
    {  
        f = NULL;
        rc = NULL;
        se = NULL;
        
	phg = NULL;
        _Trig_ptr = NULL;
        trk = NULL;
        emccont = NULL;
        emc     = NULL;
        rp      = NULL;
        rpsngl  = NULL;
        header  = NULL;
        
        eoverp = NULL;
        dep = NULL;
        emcdphi = NULL;
        emcdz = NULL;

        bbczdc = NULL;
        bbcrp = NULL;
        h_ntrk = NULL;
        h_ncharged = NULL;
        h_nclust = NULL;
        h_nemc = NULL;
        
        map_dc = NULL;
        map_emc = NULL;
        // map_rich = NULL;
        runnumber = 0;
        runselection = 0;
        systemselection = 0;

        OutFileName = outfile;
        
        return;
    }

    QA::~QA()
    { 
        delete f;  
    }

    int QA::Init(PHCompositeNode *topNode)
    {        
        se = Fun4AllServer::instance();	
        
        std::cout << "QA::Init: " << "Book bbczdc multi-histogram" << std::endl;
        bbczdc = new THmulf("bbczdc","Bbc Zdc Distr");
        bbczdc->AddAxis("cent","centrality",100, -0.5, 99.5);
        bbczdc->AddAxis("bbcs","bbc S charge",15,0.,1500.);
        bbczdc->AddAxis("bbcn","bbc N charge",15,0.,1500.);
        bbczdc->AddAxis("vtxz","vertex z position",70,-35.,35.);
        se->registerHisto( bbczdc->GetName(), bbczdc );

        std::cout << "QA::Init: " << "Book rp multi-histogram" << std::endl;
        bbcrp = new THmulf("bbcrp","BBC Reaction Plane");
        bbcrp->AddAxis("psi_bbcn","rp from BBC N",100,-TMath::Pi()/2,TMath::Pi()/2);
        bbcrp->AddAxis("psi_bbcs","rp from BBC S",100,-TMath::Pi()/2,TMath::Pi()/2);
        bbcrp->AddAxis("psi_bbc","rp from BBC N+S",100,-TMath::Pi()/2,TMath::Pi()/2);
        se->registerHisto( bbcrp->GetName(), bbcrp );
        
        std::cout << "QA::Init: " << "Book eoverp multi-histogram" << std::endl;
        eoverp = new THmulf("eoverp","E over p");
        eoverp->AddAxis("ep","e/p",100, 0, 1.5);
        eoverp->AddAxis("pt","pt",100,0.,10.);
        eoverp->AddAxis("sect","sect",8,-0.5,7.5);
        eoverp->AddAxis("charge","charge",3,-1.5,1.5);
        se->registerHisto( eoverp->GetName(), eoverp );

        std::cout << "QA::Init: " << "Book dep multi-histogram" << std::endl;
        dep = new THmulf("dep","E over p");
        dep->AddAxis("sep","sigmalized e/p",100, -5, 5);
        dep->AddAxis("pt","pt",100,0.,10.);
        dep->AddAxis("sect","sect",8,-0.5,7.5);
        dep->AddAxis("charge","charge",3,-1.5,1.5);
        se->registerHisto( dep->GetName(), dep );

        std::cout << "QA::Init: " << "Book emcdphi multi-histogram" << std::endl;
        emcdphi = new THmulf("emcdphi","emcdphi");
        emcdphi->AddAxis("emcdphi","emcdphi",100, -0.1, 0.1);
        emcdphi->AddAxis("pt","pt",100,0.,10.);
        emcdphi->AddAxis("sect","sect",8,-0.5,7.5);
        emcdphi->AddAxis("charge","charge",3,-1.5,1.5);
        se->registerHisto( emcdphi->GetName(), emcdphi );

        std::cout << "QA::Init: " << "Book emcdz multi-histogram" << std::endl;
        emcdz = new THmulf("emcdz","emcdz");
        emcdz->AddAxis("emcdz","emcdz",100, -40, 40);
        emcdz->AddAxis("pt","pt",100,0.,10.);
        emcdz->AddAxis("sect","sect",8,-0.5,7.5);
        emcdz->AddAxis("charge","charge",3,-1.5,1.5);
        se->registerHisto( emcdz->GetName(), emcdz );
        
        std::cout << "QA::Init: " << "Book hits histograms" << std::endl;     
        h_ntrk = new TH1F("h_ntrk","Track Multiplicity",200,0,200);
        se->registerHisto( h_ntrk->GetName(), h_ntrk );
        h_ncharged = new TH1F("h_ncharged","Electron/Positron Multiplicity",3,-1.5,1.5);
        se->registerHisto( h_ncharged->GetName(), h_ncharged );
        h_nclust = new TH1F("h_nclust","Cluster Multiplicity",2000,0,2000);
        se->registerHisto( h_nclust->GetName(), h_nclust );
        h_nemc = new TH1F("h_nemc","Emcal Photon Multiplicity",8,-0.5,7.5);
        se->registerHisto( h_nemc->GetName(), h_nemc );

        // DC dead map (including PC1 effect)
        std::cout << "QA::Init: " << "Book multi-histogram for DC map" << std::endl;
        map_dc = new THmulf("map_dc","map_dc");
        map_dc->AddAxis("board","board", 82*4, -1.5, 80.5);
        map_dc->AddAxis("alpha","alpha", 100, -0.9, 0.9);
        map_dc->AddAxis("arm","arm", 2, -0.5, 1.5);
        map_dc->AddAxis("side","side", 2, -0.5, 1.5);
        se->registerHisto( map_dc->GetName(), map_dc );

        // Emcal dead map
        std::cout << "QA::Init: " << "Book multi-histogram for Emcal map" << std::endl;
        map_emc = new THmulf("map_emc","map_emc");
        map_emc->AddAxis("ypos","ypos", 50, 0, 50);
        map_emc->AddAxis("zpos","zpos", 100, 0, 100);
        map_emc->AddAxis("sect","sect", 8, -0.5, 7.5);
        se->registerHisto( map_emc->GetName(), map_emc );

        // // RICH dead map
        // std::cout << "QA::Init: " << "Book multi-histogram for RICH map" << std::endl;
        // map_rich = new THmulf("map_rich","map_rich");
        // map_rich->AddAxis("phi_pmt","phi_pmt", 80, -0.5, 79.5);
        // map_rich->AddAxis("z_pmt","z_pmt", 16, -0.5, 15.5);
        // map_rich->AddAxis("arm","arm", 2, -0.5, 1.5);
        // map_rich->AddAxis("side","side", 2, -0.5, 1.5);
        // se->registerHisto( map_rich->GetName(), map_rich );
        
        rc = recoConsts::instance();
        runselection = rc->get_IntFlag("RD_RUN_SELECTION", 0);
        systemselection = rc->get_IntFlag("RD_SYSTEM_SELECTION", 0);
        
        cout<<"run "<<runselection<<" @ system "<<systemselection<<endl;

	// if ( runselection==14 && systemselection==0 ) cut = new Run14AuAuCut();

        return 0;
    }

    int QA::InitRun(PHCompositeNode *topNode)
    {
        return 0;
    }

    int QA::process_event(PHCompositeNode *topNode)
    {
        //real data
        static int ncalls = 0;
        ncalls++;
        if ((ncalls)%1000==0)
        {
            std::cout << "Event process = " << ncalls << std::endl;
        }
        
        //===================================================================
        //                          Check Nodes
        //=================================================================== 

        if ( GetNodes(topNode) == DISCARDEVENT ) return DISCARDEVENT;
        
        //===================================================================
        //                        Event selection
        //===================================================================

        if ( !trigger_selection() ) return DISCARDEVENT; 
        // if ( !cut->applyEventCut(topNode) ) return DISCARDEVENT;       
    	  
        //===================================================================
        //             calibrations for e/p, emcdphi, emcdz
        //===================================================================

        calibrations(); 

        //===================================================================
        //                         run by run QA
        //===================================================================

        run_by_run();      
        
        //===================================================================
        //                            DC map
        //===================================================================

        DC_map();

        //===================================================================
        //                           Emcal map
        //===================================================================

        Emcal_map();        

        //===================================================================
        //                           RICH map
        //===================================================================

        // RICH_map(); 

        return 0;
    }

    int QA::GetNodes(PHCompositeNode *topNode)
    {
        // std::cout<<"QA::Reading Nodes "<<std::endl;   
        phg = findNode::getClass<PHGlobal>(topNode,"PHGlobal");
        if( !phg )
        {
            std::cout<<"can't find PHGlobal "<<std::endl;
            return DISCARDEVENT;//exit(1);
        }    
        
        _Trig_ptr = findNode::getClass<TrigLvl1>(topNode, "TrigLvl1");
        if ( !_Trig_ptr )
        {
            std::cout << PHWHERE << "can't find TrigLvl1 Node" << std::endl;
            return DISCARDEVENT;;//exit(1);
        }
        
        trk = findNode::getClass<PHCentralTrack>(topNode,"PHCentralTrack");
        if( !trk )
        {
            std::cout<<"can't find PHCentraltrack "<<std::endl;
            return DISCARDEVENT;//exit(1);
        }
        
        emccont = findNode::getClass<emcClusterContainer>(topNode, "emcClusterContainer");
        if( !emccont )
        {
            std::cout<<"can't find emcClusterContainer "<<std::endl;
            return DISCARDEVENT;//exit(1);
        }

        rp = findNode::getClass<ReactionPlaneObject>(topNode, "ReactionPlaneObject");
        if( !rp )
        {
            std::cout<<"can't find ReactionPlaneObject "<<std::endl;
            return DISCARDEVENT;//exit(1);
        }

        header = findNode::getClass<RunHeader>(topNode,"RunHeader");
        if(!header) cout<<"can't find RunHeader"<<endl;
        
        return EVENT_OK;        
    }

    bool QA::trigger_selection()
    {
        // Run14@AuAu
        // 4 - BBCLL1(>1 tubes) narrowvtx
        // 5 - BBCLL1(>1 tubes)
        // 6 - BBCLL1(>1 tubes) novertex
        int TRIGGERBIT = 4;

        int trigscaled_on = _Trig_ptr->get_lvl1_trigscaled_bit(TRIGGERBIT);
 
        if ( !trigscaled_on ) return 0;

        return 1;
    }

    void QA::run_by_run()
    {
        // fill bbczdc
        float bbcqn = phg->getBbcChargeN();
        float bbcqs = phg->getBbcChargeS();
        // float bbcq = bbcqn + bbcqs;
        float bbcz = phg->getBbcZVertex();
        float percent = phg->getCentrality();
        if(percent < 0 || percent > 100) percent = -9999.; //just a safeguard
        bbczdc->Fill(1.0,percent,bbcqs,bbcqn,bbcz);

        // fill bbc rp
        rpsngl = rp->getReactionPlane(RP::calcIdCode(RP::ID_BBC, 0, 1));
        float psi_BBCS = (rpsngl) ? rpsngl->GetPsi() : -9999;
        rpsngl = rp->getReactionPlane(RP::calcIdCode(RP::ID_BBC, 1, 1));
        float psi_BBCN = (rpsngl) ? rpsngl->GetPsi() : -9999;
        rpsngl = rp->getReactionPlane(RP::calcIdCode(RP::ID_BBC, 2, 1));
        float psi_BBC = (rpsngl) ? rpsngl->GetPsi() : -9999;
        bbcrp->Fill(1.0,psi_BBCN,psi_BBCS,psi_BBC);

        // fill h_ntrk & h_ncharged
        int ntrk = 0; // # of trks pass quality=63 cut
        for (int itrk_reco = 0; itrk_reco < (int)(trk->get_npart()); ++itrk_reco)         
        {
            if ( trk->get_quality(itrk_reco)==63 ) ++ntrk;
            //cuts for single electron cut
	    float Z_GLOBAL=75,EMCDZ=20,EMCDPHI=0.05,N0=1,DISP=5,CHI2_NEP0=10,PROB=0.01,DEP0=-2,DEP1=2;
	    if ( ((trk->get_quality(itrk_reco)==31) || (trk->get_quality(itrk_reco)==51) || (trk->get_quality(itrk_reco)==63 ))
		 && (fabs(trk->get_zed(itrk_reco))< Z_GLOBAL)
		 && (fabs(trk->get_emcdz(itrk_reco))< EMCDZ)
		 && (fabs(trk->get_emcdphi(itrk_reco))< EMCDPHI)
		 && (trk->get_n0(itrk_reco)> N0)
		 && (trk->get_disp(itrk_reco)<DISP)
		 && ((trk->get_chi2(itrk_reco)/trk->get_npe0(itrk_reco))<CHI2_NEP0)
		 && (trk->get_prob(itrk_reco)>PROB)
		 && (trk->get_dep(itrk_reco)>DEP0) 
		 && (trk->get_dep(itrk_reco)< DEP1) )
	    //if ( cut->applySingleCut(trk, itrk_reco) )
            {//quality cut and loose single cuts
                h_ncharged->Fill( trk->get_charge(itrk_reco) );
            }
        }
        h_ntrk->Fill( ntrk ); // total counts in one file will be the total number of events

        // fill h_nclust & h_nemc
        h_nclust->Fill( emccont->size() ); // total counts in one file will be the total number of events
        for (int iclus = 0; iclus < (int)(emccont->size()); iclus++)
        {
            emc = emccont->getCluster(iclus);
	    //cluster cuts
	    if ( emc->chi2()< 3 )
	    // if ( cut->applyClusterCut(emc) )
            {// photon id cut
                int sect = -9999;
                //emc->arm() == 0 ? sect = emc->sector() : sect = 7 - emc->sector(); // not sure...
		sect = emc->sector();
                h_nemc->Fill(sect); 
            }
        }    
    }

    bool QA::cut_for_calibration(int itrk_reco)
    {
        // no dep cut, no E/p cut, not sure about emcdphi & emcdz
        float Z_GLOBAL=75, N0 = 1, DISP = 5, CHI2_NPE0 = 10, EMCDZ = 40, EMCDPHI = 0.1, PROB = 0.01;
        if ( (trk->get_quality(itrk_reco)!=31) && (trk->get_quality(itrk_reco)!=51) && (trk->get_quality(itrk_reco)!=63) ) return 0;
        if ( fabs(trk->get_zed(itrk_reco))>= Z_GLOBAL ) return 0;
        if ( fabs(trk->get_emcdphi(itrk_reco))>=EMCDPHI ) return 0;
        if ( fabs(trk->get_emcdz(itrk_reco))>=EMCDZ ) return 0;
        if ( trk->get_n0(itrk_reco)<=N0 ) return 0;
        if ( trk->get_disp(itrk_reco)>=DISP ) return 0;
        if ( (trk->get_chi2(itrk_reco)/trk->get_npe0(itrk_reco))>=CHI2_NPE0 ) return 0;
        if ( trk->get_prob(itrk_reco)<=PROB ) return 0;

        return 1;
    }

    void QA::calibrations()
    {
        for (int itrk_reco = 0; itrk_reco < int(trk->get_npart()); ++itrk_reco)         
        {
            if ( cut_for_calibration(itrk_reco) )
            {
                double px = trk->get_px(itrk_reco);
                double py = trk->get_py(itrk_reco);
                double pz = trk->get_pz(itrk_reco);
                double ecore = trk->get_ecore(itrk_reco);
                double p = sqrt( px*px + py*py + pz*pz );
                double pT = sqrt(px*px + py*py);
                int sect = -9999;
                trk->get_dcarm(itrk_reco) == 0 ? sect = trk->get_sect(itrk_reco) : sect = 7 - trk->get_sect(itrk_reco); 
                eoverp->Fill(1.0, ecore/p, pT, sect, trk->get_charge(itrk_reco));
                dep->Fill(1.0, trk->get_dep(itrk_reco), pT, sect, trk->get_charge(itrk_reco));
                emcdphi->Fill(1.0, trk->get_emcdphi(itrk_reco), pT, sect, trk->get_charge(itrk_reco));
                emcdz->Fill(1.0, trk->get_emcdz(itrk_reco), pT, sect, trk->get_charge(itrk_reco));
            }
        }
    }

    static float get_board(float phi_dc)
    {
        if ( phi_dc>1.57 ) return ((3.72402-phi_dc+0.008047*cos(phi_dc+0.87851))/0.01963496);
        return ((0.573231+phi_dc-0.0046*cos(phi_dc+0.05721))/0.01963496);
    }

    void QA::DC_map()
    {
        int Z_GLOBAL=75;
        for (int itrk_reco = 0; itrk_reco < (int)(trk->get_npart()); ++itrk_reco)         
        {
            if ( trk->get_quality(itrk_reco)==63 && fabs(trk->get_zed(itrk_reco))<Z_GLOBAL )
            {
                float board = get_board( trk->get_phi(itrk_reco) );
                float alpha = trk->get_alpha(itrk_reco);
                int arm  = trk->get_dcarm(itrk_reco); // arm 0 <--> e/w not sure yet
                int side = (int)( trk->get_zed(itrk_reco)>0 ); // side 0 <--> z>0
                map_dc->Fill(1.0, board, alpha, arm, side);             
            }
        }
    }

    void QA::Emcal_map()
    {
        for (int iclus = 0; iclus < (int)(emccont->size()); iclus++)
        {
            emc = emccont->getCluster(iclus);
            //if ( cut->applyClusterCut(emc) )
	    if ( emc->chi2()< 3 )
            {// photon id cut
                int sect = -9999;
                emc->arm() == 0 ? sect = emc->sector() : sect = 7 - emc->sector(); // not sure...
                int ypos = emc->iypos();
                int zpos = emc->izpos();
                map_emc->Fill(1.0, ypos, zpos, sect);
            }
        }  
    }

    // void QA::RICH_map() {   }

    int QA::End(PHCompositeNode *topNode){
        // write the output trees
        std::cout<<"creating histos for QA"<<std::endl;
        
        f = new TFile(OutFileName.c_str(),"RECREATE");
        bbczdc->Write();
        bbcrp->Write();
        eoverp->Write();
        dep->Write();
        emcdphi->Write();
        emcdz->Write();

        h_ntrk->Write();
        h_ncharged->Write();
        h_nclust->Write();
        h_nemc->Write();

        map_dc->Write();
        map_emc->Write();
        // map_rich->Write();

        f->Write();
        f->Close();

        std::cout<<"end of End"<<std::endl;
        return 0;
    }

