
#include <MINDplotter.h>


class sortTrajByHits{
public:
  bool  operator()(const Trajectory* t1,const Trajectory* t2 ){
    if (t1->nmeas() > t2->nmeas()) return true;
    return false;
  }

};
class sortTrajByLength{
public:
  bool  operator()(const Trajectory* t1,const Trajectory* t2 ){
    if (t1->length() > t2->length()) return true;
    return false;
  }

};

//*************************************************************************************
MINDplotter::MINDplotter() {
  //*************************************************************************************

}

//*************************************************************************************
MINDplotter::~MINDplotter() {
  //*************************************************************************************

}

//*************************************************************************************
void MINDplotter::initialize(string outFileName, bool patRec, bool clust,
			     bhep::prlevel vlevel) {
  //*************************************************************************************

  //bool ok = true;
  _patR = patRec;
  _clu = clust;

  level = vlevel;

  m = bhep::messenger(level);

  m.message("++Creating output root file++",bhep::VERBOSE);

  outFile = new TFile(outFileName.c_str(), "recreate");

  statTree = new TTree("tree", "Tree with pattern rec and fit data for MIND");

  define_tree_branches();

  m.message("++plotter initialized",bhep::VERBOSE);

  

  //return ok;
}

//*************************************************************************************
void MINDplotter::execute(fitter& Fit, const bhep::event& evt) {
  //*************************************************************************************
  std::cout<<"++MINDplotter execute"<<std::endl;  
  bool ok1, ok2;


  ///event informations  
  _evNo = evt.event_number();
  bhep::Point3D evVert = evt.vertex();
  _vert[0] = evVert.x(); _vert[1] = evVert.y(); _vert[2] = evVert.z();
 

  if (_clu){
    get_tru_int_type( evt );
    if ( evt.find_iproperty("Charm") )
      _charm[0] = evt.fetch_iproperty("Charm");
    else _charm[0] = 0;
    if ( evt.find_iproperty("CharmHad") )
      _charm[1] = evt.fetch_iproperty("CharmHad");
    else _charm[1] = 0;
    if ( evt.find_dproperty("Q2") )
      _Q2 = evt.fetch_dproperty("Q2");
    else _Q2 = 0.;
    if ( evt.find_dproperty("EngTrans") )
      _engTrans = evt.fetch_dproperty("EngTrans");
    if ( evt.find_iproperty("nuType") )
      _pdg[0] = evt.fetch_iproperty("nuType");
    else _pdg[0] = 0;
    if ( evt.find_iproperty("intpart") )
      _pdg[1] = evt.fetch_iproperty("intpart");
    else _pdg[1] = 0;
    if ( evt.find_iproperty("nucType") )
      _pdg[2] = evt.fetch_iproperty("nucType");
    else _pdg[2] = 0;
  }



  
  ///Event Info
  _failEvent = 0;
  _failEvent = Fit.GetFailEvent();

  // Initializing the  variable

  for (int i = 0;i<10;i++){
    _engTraj[i] = 0;
    _Fit[i] = 0;
    _reFit[i] = 0;
    _fail[i] = 0;
    _vertZ[i] = 0;
    _leng[i] = 0;
   
    ///_hadE[1] = 0;
    _reseed_count[i][0]=0;///
    _Xtent[i]=0;///

    _plns[0][i] = 0; _plns[1][i] = 0;

    _nhad[0][i] = 0; _nhad[1][i] = 0;

    _Q[0][i] = 0; _Q[1][i] = 0;  _Q[2][i] = 0;
    
    _rangqP[0][i] = 0;  _rangqP[1][i] = 0;  _rangqP[2][i] = 0;

    _hitType[0][i] = 0;  _hitType[1][i] = 0;  _hitType[2][i] = 0;  _hitType[3][i] = 0;  _hitType[4][i] = 0;

    _qP[0][i] = 0; _qP[1][i] = 0; _qP[2][i] = 0; _qP[3][i] = 0; _qP[4][i] = 0; _qP[5][i] = 0;

    /// For hadrons
    _hadE[0][i] = 0.; _hadE[1][i] = 0.;

    _hadP[0][i] = 0.;  _hadP[1][i] = 0.;  _hadP[2][i] = 0.;

    _chadP[0][i] = 0.;  _chadP[1][i] = 0.;  _chadP[2][i] = 0.;  _chadP[3][i] = 0.;

    _nhadP[0][i] = 0.;  _nhadP[1][i] = 0.;  _nhadP[2][i] = 0.;  _nhadP[3][i] = 0.;

    // _XPos[i].clear();
  }
  

  ///
  

  //_plns[0].clear();
  //_plns[1].clear();


  ///
  _XPos.clear();
  _YPos.clear();
  _ZPos.clear();
  _Edep.clear();



  for (int i = 0;i<6;i++){
    for (int j=0; j<2; j++){
      for (int k=0; k<10; k++){
	_X[i][j][k]=0;
	_Th[i][j][k]=0;
      }
    }
  }
  

  
  // Also initialize the hadron momentum and energy 
  /* for (int i = 0;i<5;i++){
    for (int j=0; j<10; j++){
    // if ( i < 4 ) { _chadP[i] = 0.; _nhadP[i] = 0.; }
    //     //else _hadE[0] = 0.;
      
      if ( i < 3 ) _hadP[i][j] = 0.;
      else _hadE[0][j] = 0.;
    }
    }*/
  
 



  
  ///vector of trajectories
  std::vector<Trajectory*> &trajs = Fit.GetTrajs();
  
  _trajs_no = trajs.size();


  
  //loop over trajectories
  for(int i=0; i<trajs.size(); i++  ){
    
    

    // _trajs_no = i; 
    //cout<<"inside MINDplotter:: Trajectory "<<i<<" and meas"<<*trajs[i]<<endl;



    ///fit and refit info  
    bool success = trajs[i]->quality("fitted");
    
    _Fit[i] = (int)success;
    _fail[i] = (int)(trajs[i]->quality("failType"));
    _reFit[i] = (int)trajs[i]->quality("reseed");
    _vertZ[i] =  trajs[i]->quality("vertZ");
    //cout<<"_Fit ="<<_Fit[i]<<" _reFit ="<<_reFit[i]<<"  fail="<<_fail[i]<<"   "<<_vertZ[i]<<endl;    
    


    ///no of planes and freeplanes for each traj  
    
    //if ( _fail[i] != 7 ){
    if ( _failEvent != 7 ){
      _visEng = Fit.get_classifier().get_vis_eng();    
      //      _plns[i][0] = (int)(trajs[i]->quality("nplanes")) ;
      //      _plns[i][1] = (int)(trajs[i]->quality("freeplanes")) ;

      _plns[0][i] = (int)(trajs[i]->quality("nplanes")) ;
      _plns[1][i] = (int)(trajs[i]->quality("freeplanes")) ;
      
      // _plns[0].push_back((int)(trajs[i]->quality("nplanes"))) ;
      // _plns[1].push_back((int)(trajs[i]->quality("freeplanes"))) ;

      //cout<<"nplanes= "<<_plns[0][i]<<" freeplanes = "<<_plns[1][i]<<endl;
      // } else { _plns[0].clear(); _plns[1].clear(); _visEng = 0; }
    } else { _plns[0][i] = 0; _plns[1][i] = -1; _visEng = 0; }
    


    ///interaction type
    //if (_fail[i] != 7)
     if ( _failEvent != 7 )
      _intType[i] = (int)(trajs[i]->quality("intType")) ;
    else _intType[i] = 7;
    
    

    
    ///true particle informations for each traj
    if ( _clu )
      ok1 = extract_true_particle2(evt, Fit, *trajs[i], i);
    else
      ok1 = extract_true_particle1(evt, Fit, *trajs[i], i);
    

    
    ///for fitting success
    if (success) {
      State ste; 
      fill_kinematics(*trajs[i], ste, i);///
      

      
      ///extrapolate upto vertex
      ok2 = extrap_to_vertex(*trajs[i], evt.vertex(), Fit, ste, i);


      
      
      ///length calculation
      _leng[i] = trajs[i]->length();
      
     
      

      ///range of q/p calculation ?
      // _rangqP[0] = Fit.GetInitialqP();///
      _rangqP[0][i] = trajs[i]->quality("initialqP");
      _rangqP[1][i] = _rangqP[0][i] - _qP[0][i];
      _rangqP[2][i] = _qP[1][i] - _rangqP[0][i];
      


      ///chi2 value of traj
      max_local_chi2( *trajs[i], i );



      ///if extrapolate to vertex is successful then fill the parameters
      if (ok2) {
	position_pulls(i);
	direction_pulls(i);
	momentum_pulls(i);
      }
    
      //hadron_direction(Fit);
    
    }
  
    //If fit not successful set rec values to zero.
    if (!success){
      _X[1][0][i] = 0; _X[1][1][i] = 0; 
      _X[2][0][i] = 0; _X[2][1][i] = 0; 
      _X[3][0][i] = 0; _X[3][1][i] = 0; 
      _X[4][0][i] = 0; _X[4][1][i] = 0; 
      _X[5][0][i] = 0; _X[5][1][i] = 0;

      _Th[1][0][i] = 0; _Th[1][1][i] = 0;
      _Th[2][0][i] = 0; _Th[2][1][i] = 0;
      _Th[3][0][i] = 0; _Th[3][1][i] = 0;
      _Th[4][0][i] = 0; _Th[4][1][i] = 0;
      _Th[5][0][i] = 0; _Th[5][1][i] = 0;

      _qP[1][i] = 0; _qP[2][i] = 0;_qP[3][i] = 0; _qP[4][i] = 0; _qP[5][i] = 0;

      _Chi[0][i] = 0; _Chi[1][i] = 0;

      _Q[1][i] = 0; 
      _Q[2][i] = 0;

    }


    /// pattern recognition informations
    if (_patR){
       
      if ( _clu )
	patternStats2( Fit, *trajs[i],i);
      //  else
      // patternStats1( Fit, *trajs[i],i );
    }


  }//end of traj loop
  


  //Fill tree event with the values.
  int fillcatch = statTree->Fill();

  //return ok1;
}


/*****************************************************************/
void MINDplotter::get_tru_int_type(const bhep::event& evt){
  
  string intName = evt.fetch_sproperty("IntType");
  
  if ( intName=="CCQE" || intName=="<QES - Weak[CC]>" )
    _truInt = 1;
  else if ( intName=="NCQE" || intName=="<QES - Weak[NC]>" )
    _truInt = 2;
  else if ( intName=="CCDIS" || intName=="<DIS - Weak[CC]>" )
    _truInt = 3;
  else if ( intName=="NCDIS" || intName=="<DIS - Weak[NC]>" )
    _truInt = 4;
  else if ( intName=="1piRes" )
    _truInt = 5;
  else if ( intName=="miscRes" || intName=="<RES - Weak[CC]>" || intName=="<RES - Weak[NC]>" )
    _truInt = 6;
  else if ( intName=="eEl" || intName=="<NuEEL - Weak[NC]>" || intName=="<NuEEL - Weak[CC]>" )
    _truInt = 7;
  else if ( intName=="muINVe" || intName=="<IMD - Weak[CC]>" )
    _truInt = 7;
  else
    _truInt = 8;
  
}

//*************************************************************************************
void MINDplotter::finalize() {
  //*************************************************************************************

  //bool ok = true;

  m.message("++Finalizing Output++",bhep::VERBOSE);
  
  outFile->Write();
  outFile->Close();

  // if(_intType==5){
  // _CAFile->Write();
  // _CAFile->Close();//}

  //return ok;
}

//*************************************************************************************
void MINDplotter::define_tree_branches() {
  //*************************************************************************************

  ///for each event
  statTree->Branch("Evt", &_evNo, "EventNo/I");
  statTree->Branch("evVertex", &_vert, "vert[3]/D");
  statTree->Branch("EvtFail", &_failEvent, "EvtFail/I");
  if (_clu) statTree->Branch("TrueInteraction",&_truInt,"truInt/I");
  statTree->Branch("NeuEng", &_nuEng, "NuEng/D");
  if (_clu) {
    statTree->Branch("Charm", &_charm, "isCharm/I:pdg/I");
    statTree->Branch("InteractionQ2", &_Q2, "Q2/D");
    statTree->Branch("EventPdg",&_pdg, "initnu/I:partPDG/I:nucType/I");
  }
  statTree->Branch("visibleEng", &_visEng, "visEng/D");

 
  ///for each traj 
  statTree->Branch("TrajectoryNo", &_trajs_no, "trajNo/I");
  statTree->Branch("TrajVertex", &_vertZ, "trajVert[trajNo]/D");
  statTree->Branch("Fitted", &_Fit, "success[trajNo]/I");
  statTree->Branch("backFit",&_reFit,"backFit[trajNo]/I");
  statTree->Branch("Fail", &_fail, "FailType[trajNo]/I");
  statTree->Branch("interaction",&_intType,"Inter[trajNo]/I");

  ///not working yet
  statTree->Branch("visEngTraj",&_engTraj, "engTraj[trajNo]/D");
  statTree->Branch("trajEngDep",&_engvar[0],"meanDep[trajNo]/D");
  statTree->Branch("trajEngVar",&_engvar[1],"engVar[trajNo]/D");
  
  //statTree->Branch("NoPlanes", &_plns, "nplanes[trajNo]/I:fplanes[trajNo]/I");
  //statTree->GetLeaf("planes")->SetAddress(_plns[0]);
  statTree->Branch("Planes", &_plns[0], "nplanes[trajNo]/I");
  statTree->Branch("FPlanes",&_plns[1], "freeplanes[trajNo]/I");
 

  ///position  
  //statTree->Branch("Position", &_X, "truPos[trajNo][2]/D:recPos[trajNo][2]/D:ErrPos[trajNo][2]/D:recPos_WVE[trajNo][2]/D:ErrPos_WVE[trajNo][2]/D:pull_X[trajNo][2]/D");///
  
  statTree->Branch("TrueXPosition", &_X[0][0], "truXPos[trajNo]/D");
  statTree->Branch("TrueYPosition", &_X[0][1], "truYPos[trajNo]/D");
  statTree->Branch("RecXPosition", &_X[1][0],"recXPos[trajNo]/D");
  statTree->Branch("RecYPosition", &_X[1][1], "recYPos[trajNo]/D");
  statTree->Branch("ErrXPosition", &_X[2][0], "ErrXPos[trajNo]/D");
  statTree->Branch("ErrYPosition", &_X[2][1], "ErrYPos[trajNo]/D");
  statTree->Branch("RecXPosition_WVE", &_X[3][0], "recXPos_WVE[trajNo]/D");
  statTree->Branch("RecYPosition_WVE", &_X[3][1], "recYPos_WVE[trajNo]/D");
  statTree->Branch("ErrXPosition_WVE", &_X[2][0], "ErrXPos_WVE[trajNo]/D");
  statTree->Branch("ErrYPosition_WVE", &_X[2][1], "ErrYPos_WVE[trajNo]/D");
  statTree->Branch("XPull", &_X[5][0], "pull_x[trajNo]/D");
  statTree->Branch("YPull", &_X[5][1], "pull_y[trajNo]/D");

  ///Directions
  //statTree->Branch("Direction", &_Th, "truTh[trajNo][2]/D:recTh[trajNo][2]/D:ErrTh[trajNo][2]/D:recTh_WVE[trajNo][2]/D:ErrTh_WVE[trajNo][2]/D:pull_Th[trajNo][2]/D");
  statTree->Branch("TrueXDirection", &_Th[0][0], "truXTh[trajNo]/D");
  statTree->Branch("TrueYDirection", &_Th[0][1], "truYTh[trajNo]/D");
  statTree->Branch("RecXDirection", &_Th[1][0], "recXTh[trajNo]/D");
  statTree->Branch("RecYDirection", &_Th[1][1], "recYTh[trajNo]/D");
  statTree->Branch("ErrXDirection", &_Th[2][0], "ErrXTh[trajNo]/D");
  statTree->Branch("ErrYDirection", &_Th[2][1], "ErrYTh[trajNo]/D");
  statTree->Branch("RecXDirection_WVE", &_Th[3][0], "recXTh_WVE[trajNo]/D");
  statTree->Branch("RecYDirection_WVE", &_Th[3][1], "recYTh_WVE[trajNo]/D");
  statTree->Branch("ErrXDirection_WVE", &_Th[4][0], "ErrXTh_WVE[trajNo]/D");
  statTree->Branch("ErrYDirection_WVE", &_Th[4][1], "ErrYTh_WVE[trajNo]/D");
  statTree->Branch("XPullTh", &_Th[5][0], "pull_Thx[trajNo]/D");
  statTree->Branch("YPullTh", &_Th[5][1], "pull_Thy[trajNo]/D");


  ///Momentum
  //statTree->Branch("Momentum", &_qP, "truqP[trajNo]/D:recqP[trajNo]/D:ErrqP[trajNo]/D:recqP_WVE[trajNo]/D:ErrqP_WVE[trajNo]/D:pull_mom[trajNo]/D");
  statTree->Branch("TruMom", &_qP[0], "truqP[trajNo]/D");
  statTree->Branch("RecMom", &_qP[1], "recqP[trajNo]/D");
  statTree->Branch("ErrMom", &_qP[2], "ErrqP[trajNo]/D");
  statTree->Branch("RecMom_WVE", &_qP[3], "recqP_WVE[trajNo]/D");
  statTree->Branch("ErrMom_WVE", &_qP[4], "ErrqP_WVE[trajNo]/D");
  statTree->Branch("MomPull", &_qP[5], "pull_mom[trajNo]/D");



  ///Charge 
  //statTree->Branch("Charge", &_Q, "truQ[trajNo]/I:recQ[trajNo]/I:ID[trajNo]/B");
  statTree->Branch("TruCharge", &_Q[0], "truQ[trajNo]/I");
  statTree->Branch("RecCharge", &_Q[1], "recQ[trajNo]/I");
  statTree->Branch("ID", &_Q[2], "ID[trajNo]/I");



  statTree->Branch("length", &_leng,"lenTraj[trajNo]/D");

  ///range of qP calculations only for success ?
  statTree->Branch("rangqP", &_rangqP[0],"rangqP[trajNo]/D");
  statTree->Branch("rangErr", &_rangqP[1],"rangErr[trajNo]/D");
  statTree->Branch("rangdiff", &_rangqP[2],"recrangediff[trajNo]/D");

  statTree->Branch("FitChiInfo", &_Chi[0], "trajChi[trajNo]/D");
  statTree->Branch("MaxFitChiInfo", &_Chi[1], "MaxLoc[trajNo]/D");

  statTree->Branch("Reseed", &_reseed_count, "reseed_count[trajNo]/I");
  statTree->Branch("extent", &_Xtent, "xtent[trajNo]/D");

  /*statTree->Branch("YPositions", &_YPos, "Y[trajNo][nhits]/D");
  statTree->Branch("ZPositions", &_ZPos, "Z[trajNo][nhits]/D");
  statTree->Branch("EngDeposit", &_Edep, "E[trajNo][nhits]/D");*/

  statTree->Branch("NoHits", &_nhits, "nhits[trajNo]/I");
  //worked statTree->Branch("XPositions", &_XPos[_trajs_no]);
  statTree->Branch("XPositions", &_XPos);
  statTree->Branch("YPositions", &_YPos);
  statTree->Branch("ZPositions", &_ZPos);
  statTree->Branch("EngDeposit", &_Edep);


  ///for hadrons
  statTree->Branch("NChadnNhad", &_nhad[0], "nChad[trajNo]/I");
  statTree->Branch("Nhad", &_nhad[1], "nNhad[trajNo]/I");
  statTree->Branch("hadronPx", &_hadP[0], "hadPx[trajNo]/D");
  statTree->Branch("hadronPy", &_hadP[1], "hadPy[trajNo]/D");
  statTree->Branch("hadronPx", &_hadP[2], "hadPz[trajNo]/D");
  statTree->Branch("EnergyTransfer", &_engTrans, "nu/D");
  statTree->Branch("hadTruEng", &_hadE[0], "truE[trajNo]/D");
  statTree->Branch("hadRecEng", &_hadE[1], "recE[trajNo]/D");

  statTree->Branch("ChadMom", &_chadP[0], "chadPx[trajNo]/D");
  statTree->Branch("ChadMom", &_chadP[1], "chadPy[trajNo]/D");
  statTree->Branch("ChadMom", &_chadP[2], "chadPz[trajNo]/D");
  statTree->Branch("ChadMom", &_chadP[3], "chadE[trajNo]/D");

  statTree->Branch("NhadMom", &_nhadP[0], "nhadPx[trajNo]/D");
  statTree->Branch("NhadMom", &_nhadP[1], "nhadPy[trajNo]/D");
  statTree->Branch("NhadMom", &_nhadP[2], "nhadPz[trajNo]/D");
  statTree->Branch("NhadMom", &_nhadP[3], "nhadE[trajNo]/D");
    
    //statTree->Branch("hadDir", &_haddot, "dotProd/D");
    //statTree->Branch("truMuHitIndex", &_truMuHitIndex, "truMuHitInd[2]/I");


  ///HitBreakUp
  //statTree->Branch("HitBreakDown", &_hitType, "nTruMu/I:nInMu/I:nMuInMu/I:nFitN/I:nhad/I");
  statTree->Branch("TruMu", &_hitType[0], "nTruMu[trajNo]/I");
  statTree->Branch("InMu", &_hitType[1], "nInMu[trajNo]/I");
  statTree->Branch("MuInMu", &_hitType[2], "nMuInMu[trajNo]/I");
  statTree->Branch("FitN", &_hitType[3], "nFitN[trajNo]/I");
  statTree->Branch("hadN", &_hitType[4], "nhad[trajNo]/I");
  
  

  statTree->Branch("MuHits", &_mus, "truMu[trajNo][nhits]/B");
  statTree->Branch("CandHits", &_cand, "inMu[trajNo][nhits]/B");
  statTree->Branch("FittedNodes",&_node,"fitNode[trajNo][nhits]/B");
  statTree->Branch("HadronNodes",&_had,"hadN[trajNo][nhits]/B");
  statTree->Branch("PatRecChi", &_pChi, "maxChiMu/D:MinChiHad/D:MaxConsecHol/D");
  
  //statTree->Branch("hitsInTraj", &_traj_hits, "traj_hits/I");

}

//*************************************************************************************
bool MINDplotter::fill_kinematics(const Trajectory& traj, State& ste, const int trajno){ 
  //*************************************************************************************

  //State ste;

  m.message("++Fill vertex kinematics++",bhep::VERBOSE);


  ste = traj.node(traj.first_fitted_node()).state();

  //Grab fitted vertex information.
  EVector v = ste.hv().vector();
  EMatrix C = ste.hv().matrix();

  //Reconstructed x,y position.
  _X[3][0][trajno] = v[0]; _X[3][1][trajno] = v[1];
 

  //Corresponding Error.
  if (C[0][0]>0)
    _X[4][0][trajno] = sqrt(C[0][0]);
  if (C[1][1]>0)
    _X[4][1][trajno] = sqrt(C[1][1]);
  
 
  //Reconstructed q/P.
  if (v[5] !=0) _Q[1][trajno] = (int)( v[5]/fabs(v[5]) );
  _qP[3][trajno] = v[5];
  
  //cout<<"recqP_WVE="<< _qP[3][trajno]<<endl;
  //Corresponding Error.
  if (C[5][5]>0){
    _qP[4][trajno] = sqrt(C[5][5]);
    
    //momentum pull
    _qP[5][trajno] = (_qP[3][trajno] - _qP[0][trajno])/_qP[4][trajno];
  }
  
  //Correctly ID'd charge?.
  if (_Q[0][trajno] == _Q[1][trajno]) _Q[2][trajno] = true;
  else _Q[2][trajno] = false;
  
  //cout<<" fill_kinematics: truQ="<<_Q[0][trajno]<<"  recQ="<<_Q[1][trajno]<<"  ID="<<_Q[2][trajno]<<endl;
  //Reconstructed direction.
  _Th[3][0][trajno] = v[3]; _Th[3][1][trajno] = v[4];
  
  //Corresponding error.
  if (C[3][3]>0){
    _Th[4][0][trajno] = sqrt(C[3][3]);
    //direction pull
    _Th[5][0][trajno] = (_Th[3][0][trajno] - _Th[0][0][trajno])/_Th[4][0][trajno] ;}
  
  
  if (C[4][4]>0){
    _Th[4][1][trajno] = sqrt(C[4][4]);
    //direction pull
    _Th[5][1][trajno] = (_Th[3][1][trajno] - _Th[0][1][trajno])/_Th[4][1][trajno] ;
  }
  
  
}

//*************************************************************************************
bool MINDplotter::extrap_to_vertex(const Trajectory& traj, 
				   const bhep::Point3D& vertexLoc,
				   fitter& fitObj, State& ste, const int trajno) {
  //*************************************************************************************

  m.message("++Extrapolation function, Finding best fit to vertex++",bhep::VERBOSE);
 

  ///set the state
  ste = traj.node(traj.first_fitted_node()).state();

  EVector pos(3,0); pos[2] = vertexLoc.z();
  EVector axis(3,0); axis[2] = 1;

  double R1 = 1000000000;
  double R2 = 0;
  double l;

  Ring surf(pos, axis, R1, R2);

  fitObj.man().geometry_svc().setup().add_surface("Detector","vertex",&surf);
  bool ok = fitObj.man().navigation_svc().propagate(surf,ste,l);
  fitObj.man().geometry_svc().setup().remove_surface("vertex");

  //Convert to slopes representation.
  fitObj.man().model_svc().model(RP::particle_helix)
    .representation().convert(ste, RP::slopes_curv_z);

  //Grab fitted vertex information.
  vert = ste.hv().vector();
  vertMat = ste.hv().matrix();
  // cout<<"true  vertex ="<<vert[2]<<"\n"<<endl;

  return ok;
}


//**************************************************************************************
void MINDplotter::position_pulls(const int trajno) {
  //**************************************************************************************

  //Function to calculate the pull for a measurement
  m.message("++Calculating position pulls++",bhep::VERBOSE);
  
  //Reconstructed x,y position.
  _X[1][0][trajno] = vert[0]; _X[1][1][trajno] = vert[1];
  
  // cout<<"recPos ="<<_X[1][0][trajno]<<endl;
  
  
  //Corresponding Error.
  if (vertMat[0][0]>0){
    _X[2][0][trajno] = sqrt(vertMat[0][0]);
    
    
    //pull of X
    _X[5][0][trajno] = (_X[1][0][trajno] - _X[0][0][trajno])/_X[2][0][trajno] ;
  }
  if (vertMat[1][1]>0){
    _X[2][1][trajno] = sqrt(vertMat[1][1]);
    
    
    //pull of Y
    _X[5][1][trajno] = (_X[1][1][trajno] - _X[0][1][trajno])/_X[2][1][trajno] ;
  }
  
  //cout<<"rec position ="<< _X[1][0][trajno]<<"   "<< _X[1][1][trajno]<<endl;
  
}

//**************************************************************************************
void MINDplotter::momentum_pulls(const int trajno) {
  //**************************************************************************************

  ///Function to calculate momentum pulls.
  m.message("++Calculating momentum pulls++",bhep::VERBOSE);

  //Reconstructed q/P.
  //tapasi if (vert[5] !=0) _Q[1] = (int)( vert[5]/fabs(vert[5]) );
  _qP[1][trajno] = vert[5];

  //Corresponding Error.
  if (vertMat[5][5]>0)
    _qP[2][trajno] = sqrt(vertMat[5][5]);

  /*tapasi //Correctly ID'd charge?.
    if (_Q[0] == _Q[1]) _Q[2] = true;
    else _Q[2] = false;*/
  
}

//**************************************************************************************
void MINDplotter::direction_pulls(const int trajno) {
  //**************************************************************************************

  //Reconstructed direction.
  _Th[1][0][trajno] = vert[3]; _Th[1][1][trajno] = vert[4];

  //Corresponding error.
  if (vertMat[3][3]>0)
    _Th[2][0][trajno] = sqrt(vertMat[3][3]);
 
  if (vertMat[4][4]>0)
    _Th[2][1][trajno] = sqrt(vertMat[4][4]);

}

//*************************************************************************************
bool MINDplotter::extract_true_particle1(const bhep::event& evt, fitter &Fit, const Trajectory& traj, const int trajno) {
  //*************************************************************************************

  /* sets true particle momentum for calculation and returns a reference
     to the particle */

  //_nuEng = atof( evt.fetch_property("Enu").c_str() ) * bhep::GeV;
  _nuEng = evt.fetch_dproperty("nuEnergy") * bhep::MeV;
  // for (int jj=0;jj<4;jj++)
  //     _hadP[jj] = 0.;

  const vector<bhep::particle*> Pospart = evt.true_particles();

  int count = 0;
  for (int iParts=0;iParts < (int)Pospart.size();iParts++){
    if (Pospart[iParts]->name().compare("mu-")==0){
      _Q[0][trajno] = -1;
      _truPart = Pospart[iParts];
      count++;
    } 
    else if (Pospart[iParts]->name().compare("mu+")==0){
      _Q[0][trajno] = 1;
      _truPart = Pospart[iParts];
      count++;
    } 
    if (Pospart[iParts]->name().compare("Hadronic_vector")==0){
      _hadP[0][trajno] = Pospart[iParts]->px();
      _hadP[1][trajno] = Pospart[iParts]->py();
      _hadP[2][trajno] = Pospart[iParts]->pz();
      _hadE[0][trajno] = atof( Pospart[iParts]->fetch_property("HadE").c_str() ) * bhep::GeV;
      // _hadP[3] = Pospart[iParts]->fetch_dproperty("HadE") * bhep::GeV;
    }
  }
  //STUFF ABOVE FOR REDESIGN.!!!!!
  _nhits[trajno] = traj.nmeas();
  
  for (int iHits = 0;iHits < _nhits[trajno];iHits++){
    

    //_XPos[trajno][iHits] = Fit.GetMeas(iHits)->vector()[0];
    _YPos[trajno][iHits] = Fit.GetMeas(iHits)->vector()[1];
    _ZPos[trajno][iHits] = Fit.GetMeas(iHits)->position()[2];

    if (!_patR)
      if ( traj.node(iHits).status("fitted") )
	_hitType[3][trajno]++;
  }

  if (count == 0) {
    //cout << "No particles of muon or antimuon type in file" << endl;
    _Q[0][trajno] = 0;
    _qP[0][trajno] = 0;
    _Th[0][0][trajno] = _Th[0][1][trajno]= 0;
    _X[0][0][trajno] = _X[0][1][trajno] =0;///
    return false;
  }
  
  //Set true values of muon mom. etc.
  //True q/P.
  _qP[0][trajno] = _Q[0][trajno]/_truPart->p();
  //True direction.
  _Th[0][0][trajno]= _truPart->px()/_truPart->pz();
  _Th[0][1][trajno]= _truPart->py()/_truPart->pz();
  //True x,y position.
  _X[0][0][trajno] = evt.vertex().x(); _X[0][1][trajno] = evt.vertex().y();
  
 
  return true;
}




//**********************************************************************************
bool MINDplotter::extract_true_particle2(const bhep::event& evt, fitter &Fit, const Trajectory& traj, const int trajno) {
  //*************************************************************************************

  bool primNu = false;
  _nuEng = evt.fetch_dproperty("nuEnergy") * bhep::MeV;
  
  const vector<bhep::particle*> Pospart = evt.true_particles();

  int count = 0;
  for (int iParts=(int)Pospart.size()-1;iParts >= 0;iParts--){
    
    if ( Pospart[iParts]->fetch_sproperty("CreatorProcess") == "none" &&
	 Pospart[iParts]->find_sproperty("PrimaryLepton") ){
      
      if ( Pospart[iParts]->name() == "mu+" ){
	_truPart = Pospart[iParts];
	_Q[0][trajno] = 1;
	count++;
      } else if ( Pospart[iParts]->name() == "mu-" ){
	_truPart = Pospart[iParts];
	_Q[0][trajno] = -1;
	count++;
      }
    } else if ( Pospart[iParts]->fetch_sproperty("CreatorProcess") == "none" &&
		!Pospart[iParts]->find_sproperty("PrimaryLepton") )
      add_to_hads( *Pospart[iParts], trajno );///
  }



  std::vector<double> hadInf = evt.fetch_dvproperty("had4vec");
 


  for (int itn = 0;itn < 4;itn++)
    _hadP[itn][trajno] = hadInf[itn];

  /*_nhits = Fit.get_nMeas();
    for (int iHits = 0;iHits < _nhits;iHits++){
    
    _XPos[iHits] = Fit.GetMeas(iHits)->vector()[0];
    _YPos[iHits] = Fit.GetMeas(iHits)->vector()[1];
    _ZPos[iHits] = Fit.GetMeas(iHits)->position()[2];
    _Edep[iHits] = Fit.GetMeas(iHits)->get_eng();*/
    
    
  ///
  _nhits[trajno] = traj.nmeas();



  ///temporary vectors
  std::vector<double> vxpos;
  std::vector<double> vypos;
  std::vector<double> vzpos;
  std::vector<double> vedep;

  for (int iHits = 0;iHits < _nhits[trajno];iHits++){
    

    
    //cluster *mnt = static_cast<cluster*>(&(traj.measurement(iHits)));
    //_Edep[trajno][iHits] = mnt->get_eng();

    /*
    _XPos[trajno][iHits] = traj.measurement(iHits).vector()[0];
    _YPos[trajno][iHits] = traj.measurement(iHits).vector()[1];
    _ZPos[trajno][iHits] = traj.measurement(iHits).position()[2];
    _Edep[trajno][iHits] = meas.hv("energy").vector()[0];

    //_XPos[trajno].push_back(meas.vector()[0]);
    _XPos[trajno][iHits] = meas.vector()[0];
    _YPos[trajno][iHits] = meas.vector()[1];
    _ZPos[trajno][iHits] = meas.position()[2];
    */


    const Measurement& meas = traj.measurement(iHits);

    
    vxpos.push_back(meas.vector()[0]);
    vypos.push_back(meas.vector()[1]);
    vzpos.push_back(meas.position()[2]);
    vedep.push_back(meas.hv("energy").vector()[0]);

    
    //cout<<"_Xpos ="<<  _XPos[trajno][iHits]<<"  _Ypos ="<<  _YPos[trajno][iHits]<<"   Edep="<<_Edep[trajno][iHits]<<endl;    
    
    if (!_patR)
      if ( traj.node(iHits).status("fitted") )
	_hitType[3][trajno]++;
  }
  

  ///Fill the Tree variables for hits
  _XPos.push_back(vxpos);
  _YPos.push_back(vypos);
  _ZPos.push_back(vzpos);
  _Edep.push_back(vedep);
  


  if (count == 0) {
    //cout << "No particles of muon or antimuon type in file" << endl;
    _Q[0][trajno] = 0;
    _qP[0][trajno] = 0;
    _Th[0][0][trajno] = _Th[0][1][trajno]= 0;
    _X[0][0][trajno] = _X[0][1][trajno] = 0;///

     
    return false;
  }
  
 

  //Set true values of muon mom. etc.

  //True q/P.
  _qP[0][trajno] = _Q[0][trajno]/_truPart->p();

  //True direction.
  _Th[0][0][trajno]= _truPart->px()/_truPart->pz();
  _Th[0][1][trajno]= _truPart->py()/_truPart->pz();

  //True x,y position.
  _X[0][0][trajno] = evt.vertex().x(); 
  _X[0][1][trajno] = evt.vertex().y();
 
  //_X cout<<"tru position ="<< _X[0][0][trajno]<<"   "<< _X[0][1][trajno]<<endl;
  // std::cout<<"end of extract_true_particle"<<std::endl;
  return true;
}



/**************************************************************/
void MINDplotter::add_to_hads(const bhep::particle& part, const int trajno){
  /**************************************************************/

  if ( part.charge() != 0 ){
    _chadP[0][trajno] += part.px();
    _chadP[1][trajno] += part.py();
    _chadP[2][trajno] += part.pz();
    _chadP[3][trajno] += part.e();
    
    _nhad[0][trajno]++;
  } else {
    _nhadP[0][trajno] += part.px();
    _nhadP[1][trajno] += part.py();
    _nhadP[2][trajno] += part.pz();
    _nhadP[3][trajno] += part.e();
    
    _nhad[1][trajno]++;
  }
  // _hadP[0] += part.px();
  //   _hadP[1] += part.py();
  //   _hadP[2] += part.pz();

  //   _hadE[0] += part.e();

}

//*************************************************************************************
void MINDplotter::hadron_direction(fitter& fit) {
  //*************************************************************************************
  
  // double normal;
  //   EVector fitunit;
  //   normal = sqrt(pow(_hadP[0],2)+pow(_hadP[1],2)+pow(_hadP[2],2));
  
  //   if ( _nhits >= 2 ){
  //     fitunit = fit.get_had_unit();
    
  //     _haddot = fitunit[0]*(_hadP[0]/normal)
  //       +fitunit[1]*(_hadP[1]/normal)
  //       +fitunit[2]*(_hadP[2]/normal);
    
  //   } else _haddot = 99;
  
  //_hadE[1] = fit.get_had_eng();
  
}

//*************************************************************************************
void MINDplotter::max_local_chi2(const Trajectory& traj, const int trajno) {
  //*************************************************************************************

  m.message("++Finding trajectory local chi max++",bhep::VERBOSE);
  
  size_t nNodes = traj.size();
  double trajMax = 0;
  


  for (size_t iNode = 0;iNode < nNodes;iNode++){
    
    if ( traj.node(iNode).qualitymap().has_key("predicted") )
      trajMax = TMath::Max(trajMax, traj.node(iNode).quality("predicted") );

  }

  _Chi[0][trajno] = traj.quality();
  _Chi[1][trajno] = trajMax;

}

//****************************************************************************************
void MINDplotter::patternStats1(fitter& Fit, const Trajectory& traj, const int trajno) {
  //****************************************************************************************
  //Event classifier version.
  //_nhits = Fit.get_nMeas();
  const dict::Key candHit = "inMu";
  const dict::Key engDep = "E_dep";
  int nNode = 0;
  /// if ( Fit.CheckReseed() ) nNode = (int)traj.size()-1;
  bool isMu;
  
  for (int iHits = 0;iHits < _nhits[trajno];iHits++){
  
    HitIndex = -1;//tapasi 
    if (Fit.GetMeas(iHits)->name("MotherParticle").compare("mu+")==0
	|| Fit.GetMeas(iHits)->name("MotherParticle").compare("mu-")==0){
      isMu = true;
      _mus[trajno][iHits] = true;
      _hitType[0][trajno]++;

      //tapasi
      HitIndex++;
      if(HitIndex == 0) _truMuHitIndex[0] = iHits;
    
      
       
    }
    else {_mus[trajno][iHits] = false; isMu = false;}
    
    // if ( _fail[trajno] != 7 && Fit.GetMeas(iHits)->names().has_key(candHit) ){
    if ( _failEvent != 7 && Fit.GetMeas(iHits)->names().has_key(candHit) ){
      if ( Fit.GetMeas(iHits)->name(candHit).compare("True")==0 ){//has_key(candHit) ){
	_cand[trajno][iHits] = true;
	_hitType[1][trajno]++;
	_engTraj[trajno] += bhep::double_from_string( Fit.GetMeas(iHits)->name(engDep) ) * bhep::GeV;
	if ( isMu ) _hitType[2][trajno]++;
	//tapasi
	HitIndex++;     
	if(HitIndex == 0) _truMuHitIndex[1] = iHits;

	
	
	if ( traj.node(nNode).status("fitted") && _fail[trajno]!=1 && _fail[trajno]<4 ){	
	  _node[trajno][iHits] = true; _hitType[3][trajno]++; }
	else _node[trajno][iHits] = false;

	nNode++;


      }
      else { _node[trajno][iHits] = false; _cand[trajno][iHits] = false; }
    } else if ( _failEvent != 7) { _node[trajno][iHits] = false; _cand[trajno][iHits] = false; }
    //else if ( _fail[trajno] != 7) { _node[trajno][iHits] = false; _cand[trajno][iHits] = false; }

  }
  
  // if (_fail[trajno] != 7){
   if (_failEvent != 7){
    _pChi[0] = Fit.get_classifier().get_PatRec_Chis()[0];
    _pChi[1] = Fit.get_classifier().get_PatRec_Chis()[1];
    _pChi[2] = Fit.get_classifier().get_PatRec_Chis()[2];

    _engvar[1][trajno] = 0;
    for(int ii=0;ii<_hitType[1][trajno];ii++)
      _engvar[1][trajno] += pow( bhep::double_from_string( traj.node(ii).measurement().name(engDep) ) * bhep::GeV - _engTraj[trajno]/_hitType[1][trajno], 2);
    _engvar[1][trajno] /= (_hitType[1][trajno]-1);
    _engvar[0][trajno] = _engTraj[trajno]/_hitType[1][trajno];
  } else {_engvar[1][trajno] = -1; _engvar[0][trajno] = -1;}
  
  // for (int iclear = _nhits;iclear<300;iclear++){
  //     _mus[iclear] = false; _cand[iclear] = false;
  //     _node[iclear] = false;
  //   }

}


/**********************************************************************************/
void MINDplotter::patternStats2(fitter& Fit, const Trajectory& traj, int trajno) {
  /**********************************************************************************/ 


  const dict::Key candHit = "inMu";
  const dict::Key hadHit = "inhad";
  int nNode = 0;
  HitIndex=-1; ///

  ///
  if ( _reFit[trajno] )_reseed_count[trajno] [0]++;
  else _reseed_count[trajno] [0]=0;



  bool isMu;
  //std::vector<cluster*> inMuC;
  std::vector<Measurement*> inMuC;///

  //std::vector<cluster*> hits = Fit.GetMeasVec();
  

  ///creat the vector<cluster> from Measurements
  std::vector<Measurement*> hits = traj.measurements();


  ///loop over hits
  for (int iHits = 0;iHits < _nhits[trajno];iHits++){
    
    
    ///create a measurement
    //const Measurement& meas = traj.measurement(iHits);

    
   
    //if ( hits[iHits]->get_mu_prop() > 0.8 ){//good number??
    if ( hits[iHits]->hv("MuonProp").vector()[0] > 0.8 ){//good number??
      
      isMu = true;
      _mus[trajno][iHits] = true;
      _hitType[0][trajno]++;



      ///
      HitIndex++; 
      if(HitIndex == 0) _truMuHitIndex[0] = iHits;
      
       
      
    }
    else {_mus[trajno][iHits] = false; isMu = false;}

    // if ( _fail[trajno] != 7 && hits[iHits]->names().has_key(candHit) ){
     if ( _failEvent != 7 && hits[iHits]->names().has_key(candHit) ){
    
      
      if ( hits[iHits]->name(candHit).compare("True")==0 ){
     
	_cand[trajno][iHits] = true;
	_had[trajno][iHits] = false;
	_hitType[1][trajno]++;
	///_engTraj[trajno] += hits[iHits]->get_eng() * bhep::MeV;
	_engTraj[trajno] +=(hits[iHits]->hv("energy").vector()[0])* bhep::MeV;
	cout<<"	_engTraj[trajno]="<<	_engTraj[trajno]<<endl;
	inMuC.push_back( hits[iHits] );
	if ( isMu ) _hitType[2][trajno]++;

	///
	HitIndex++;	
	if(HitIndex == 0) _truMuHitIndex[1] = iHits;
	
    	
	if ( traj.node(nNode).status("fitted") && _fail[trajno]!=1 && _fail[trajno]<4 ){	
	  _node[trajno][iHits] = true; _hitType[3][trajno]++; }
	else _node[trajno][iHits] = false;
	nNode++;
      } else {
	_node[trajno][iHits] = false; _cand[trajno][iHits] = false;
	if ( hits[iHits]->names().has_key(hadHit) ){
	  _hitType[4][trajno]++;
	  _had[trajno][iHits] = true;
	} else _had[trajno][iHits] = false;

      }

    } else if ( _failEvent != 7) {
      _node[trajno][iHits] = false; _cand[trajno][iHits] = false;
      if ( hits[iHits]->names().has_key(hadHit) ){
	_hitType[4][trajno]++;
	_had[trajno][iHits] = true;
      } else _had[trajno][iHits] = false;
    }
  }


  ///pattern Recognition Chi2
  if (_failEvent != 7){
    _pChi[0] = Fit.get_classifier().get_PatRec_Chis()[0];
    _pChi[1] = Fit.get_classifier().get_PatRec_Chis()[1];
    _pChi[2] = Fit.get_classifier().get_PatRec_Chis()[2];

    _engvar[1][trajno] = 0;



    for(int ii=0;ii<_hitType[1][trajno];ii++)
      //_engvar[1][trajno] += pow( inMuC[ii]->get_eng()*bhep::MeV - _engTraj[trajno]/_hitType[1][trajno], 2);
      _engvar[1][trajno] += pow(( inMuC[ii]->hv("energy").vector()[0])*bhep::MeV - _engTraj[trajno]/_hitType[1][trajno], 2);
    _engvar[1][trajno] /= (_hitType[1][trajno]-1);
    _engvar[0][trajno] = _engTraj[trajno]/_hitType[1][trajno];
  } else {_engvar[1][trajno] = -1; _engvar[0][trajno] = -1;}


  ///clear the vectors
  hits.clear();
  inMuC.clear();

  ///
  _Xtent[trajno] =  (double)(traj.quality("xtent"));  ;
  // _traj_hits[trajno] = traj.nmeas();
  
  }

//***********************************************************************/
/*void MINDplotter::patternStats2(fitter& Fit, const Trajectory& traj, const int trajno) {
  
const dict::Key candHit = "inMu";
const dict::Key hadHit = "inhad";
int nNode = 0;
HitIndex=-1; ///
if ( _reFit[trajno] ) nNode = (int)traj.size()-1;

///
if ( _reFit[trajno] )_reseed_count[trajno][0]++;
else _reseed_count[trajno][0]=0;

bool isMu;
int true_muon = 0;///

std::vector<cluster*> hits = Fit.GetMeasVec();
std::vector<cluster*> inMuC;
  

//total true muon
for (int iHits = 0;iHits < _nhits[trajno];iHits++){
if ( traj.measurement(iHits).name(candHit).compare("True")==1 )
true_muon ++;
}
  


///loop over hits
for (int iHits = 0;iHits < _nhits[trajno];iHits++){
//  cout<<"   hit no="<<iHits<<endl;  
if ( hits[iHits]->get_mu_prop() > 0.8 ){//good number??
   
isMu = true;
_mus[iHits] = true;
_hitType[0][trajno]++;
      
HitIndex++; 
if(HitIndex == 0){ _truMuHitIndex[0] = iHits;
//    std::cout<<"HitIndex for true mu, 2 ="<<_truMuHitIndex[0]<<std::endl;
}    
 
}
else {_mus[iHits] = false; isMu = false;}

///  if ( _fail != 7 && hits[iHits]->names().has_key(candHit) ){

if ( _fail[trajno] != 7 && traj.measurement(iHits).names().has_key(candHit) ){
if ( hits[iHits]->name(candHit).compare("True")==0 ){
_cand[iHits] = true;
_had[iHits] = false;
_hitType[1][trajno]++;
_engTraj[trajno] += hits[iHits]->get_eng() * bhep::MeV;
inMuC.push_back( hits[iHits] );
if ( isMu ) _hitType[2][trajno]++;

//tapasi
HitIndex++;	
if(HitIndex == 0) _truMuHitIndex[1] = iHits;
//std::cout<<"HitIndex for cand mu, 2 ="<<_truMuHitIndex[1]<<std::endl;
    	
if ( traj.node(nNode).status("fitted") && _fail[trajno]!=1 && _fail[trajno]<4 ){	
_node[iHits] = true; _hitType[3][trajno]++; }
else _node[iHits] = false;
if ( _reFit[trajno] ) nNode--;
else nNode++;
} else {
_node[iHits] = false; _cand[iHits] = false;
if ( hits[iHits]->names().has_key(hadHit) ){
_hitType[4][trajno]++;
_had[iHits] = true;
} else _had[iHits] = false;

}

} else if ( _fail[trajno] != 7) {
_node[iHits] = false; _cand[iHits] = false;
/// if ( hits[iHits]->names().has_key(hadHit) ){
if ( traj.measurement(iHits).names().has_key(hadHit) ){
_hitType[4][trajno]++;
_had[iHits] = true;
} else _had[iHits] = false;
}
}

if (_fail[trajno] != 7){
_pChi[0] = Fit.get_classifier().get_PatRec_Chis()[0];
_pChi[1] = Fit.get_classifier().get_PatRec_Chis()[1];
_pChi[2] = Fit.get_classifier().get_PatRec_Chis()[2];

_engvar[1][trajno] = 0;
for(int ii=0;ii<_hitType[1][trajno];ii++)
_engvar[1][trajno] += pow( inMuC[ii]->get_eng()*bhep::MeV - _engTraj[trajno]/_hitType[1][trajno], 2);
_engvar[1][trajno] /= (_hitType[1][trajno]-1);
_engvar[0][trajno] = _engTraj[trajno]/_hitType[1][trajno];
} else {_engvar[1][trajno] = -1; _engvar[0][trajno] = -1;}

hits.clear();
inMuC.clear();

///
_Xtent[trajno] =  (double)(traj.quality("xtent"));  ;
// _traj_hits = traj.nmeas();
  
}*/

