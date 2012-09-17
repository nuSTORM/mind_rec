
#include <fitter.h>
#include <CLHEP/Units/PhysicalConstants.h>
#include <TMath.h>


//*************************************************************
fitter::fitter(const bhep::gstore& pstore,bhep::prlevel vlevel){
  //*************************************************************
  
  level = vlevel;
  
  store = pstore;
  
  m = bhep::messenger(level);

  m.message("++fitter Messenger generated++",bhep::VERBOSE);


}

//*************************************************************
fitter::~fitter() {
  //*************************************************************
 
  //Reset();
  //if ( _doClust ) delete _clusters;
  
}

//*************************************************************
//bool fitter::initialize(const bhep::sstore& run_store) {
void fitter::Initialize() {
  //*************************************************************
    
  m.message("+++ fitter init  function ++++",bhep::VERBOSE);
  
  _hadunit = EVector(3,0);
  _hadunit[2] = 1.;
  _detect = store.fetch_sstore("detect");


  // read parameters
  ReadParam();

  
  // initialize geometry
  geom.init(store, level);

  //Instantiate recpack manager.
  MINDfitman::instance().set_man_parameters( store, geom.setup() );
  
  if (_X0 == 0){
    man().model_svc().enable_noiser(model, RP::ms, false);
  }



  //If required make the clustering object.
  if ( _doClust )
    _clusters = new hit_clusterer( store );



  ///initialize classifier
  get_classifier().initialize( store, level, geom.get_Fe_prop() );

  m.message("+++ End of init function ++++",bhep::VERBOSE);
  
  //return true;
}

//*************************************************************
bool fitter::Execute(bhep::particle& part,int evNo){
  //*************************************************************
  std::cout<<"+++I am fitter execute"<<std::endl;    
  m.message("+++ fitter execute function ++++",bhep::VERBOSE);
  
  bool ok = true;/// 
  _reseed_ok = false;
  _reseed_called = false;
  _fitted = false;
  _pr_count = 0;
  
 
  ///create clusters or fill measurement vector
  ok = CreateMeasurements(part);
  

  ///if pattern recognition 
  if (patternRec) {
    if ((int)_meas.size() < _min_seed_hits) {
      _failEvent = 7; 
      ok = false; 
    }
  } 
  else{ 
    ///for single traj
    Trajectory* straj = new Trajectory();
    ok = CreateSingleTrajectory(*straj);
    _trajs.push_back(straj);
  }
  
  if (!ok){
    cout<<"CreateMeasurements not ok"<<endl;
    return true;
  }

  //Sort in increasing z here when classifier up and running.!!!
  sort( _meas.begin(), _meas.end(), forwardSorter() );

  ///if pattern recognition and recTrajectory is ok
  if (patternRec){
        
    /// execute event classification
    get_classifier().execute( _meas, _trajs, _hadmeas);
     
    ///sort the hadrons
    sort( _hadmeas.begin(), _hadmeas.end(), reverseSorter() );
    
    ///PR seeds for all the trajectories from classifier
    _vPR_seed = get_classifier().get_patRec_seed_vector();
    
  }
  /// for non PR track need to set tracks infos separately
  /*else if(ok) _trajs.push_back(_traj);
    
  */
  
  /// loop over trajectories 
  for (unsigned int i=0; i< _trajs.size(); i++){ 
    
    m.message("inside traj loop::if classifier ok, traj no =",i,"*********",bhep::DETAILED); 
    ///
    _fitted = false;
    _reseed_ok = false;
    _reseed_called = false;
    _failType = 0;
    _intType = 0; //set to 'success' before run to avoid faults in value.
    ok = true;///track finded by PR or CA ??
    
    
    
    /// Get the trajectory
    Trajectory& traj = *(_trajs[i]);
    
    //cout<<"from classifier traj.nmeas ="<<traj<<endl;
    cout<<"fitter::PR, size = "<< _vPR_seed.size()<<" & trajno="<<i<<"  nmeas ="<<traj.nmeas()<<endl;
    
    //get traj informations
    int nplanes = 0, freeplanes = 0;
    double xtent = 0, vertZ =0;
    
    //Get traj info before fitting
    nplanes = (int)(traj.quality("nplanes"));
    freeplanes = (int)(traj.quality("freeplanes"));
    _intType = (int)(traj.quality("intType"));
    _failType = (int)(traj.quality("failType"));
    xtent = (double)(traj.quality("xtent"));
    vertZ = (double)(traj.quality("vertZ"));
    
    
    // cout<<"inside fitter:: from classifier traj.nmeas ="<<*_traj<<endl;
    
    cout<<" fitter: For traj no ="<<i<<"  intType ="<<_intType<<"  failType="<<_failType<<endl;
                
    ///sort the nodes in increasing Z (event when PR is not running) 
    traj.sort_nodes(1);
    //	sort( _hadmeas.begin(), _hadmeas.end(), reverseSorter() );
        
    
    
    ///if the traj finding fails during event_classification CA/PR anyone 
    if(_failType==4 || _failType==5 || _failType==6) ok = false;
    
    
    //SetFit mode to manager.   
    if(ok) 
      MINDfitman::instance().fit_mode();
  
    
    ///track found by finder (CA or PR not failed)
    State seedState;
    if (ok) {
      ok = CheckValidTraj(traj);
      
      /// seed for Fit 
      if ( ok ) ComputeSeed(traj,seedState);
    }
    
    
    ///fit the trajectory 
    if (ok) {
      /// if (ok)cout<<"if classifier3="<<endl; 
      _fitted = FitTrajectory(traj,seedState,i);
            

      cout<<"- traj node0="<<*(traj.nodes()[0])<<endl;
      
      cout<<"- copied trajectory =" << traj<<endl;


      ///length of the traj
      double length;
      ok = man().matching_svc().compute_length(traj, length);///
                  
      if (_fitted) 
	if (_failType!=3) _failType = 0;      
    }
    
    
    ///assign quality for each trajectory
    
    traj.set_quality("failType",_failType);
    traj.set_quality("intType",_intType);
    traj.set_quality("nplanes",nplanes);
    traj.set_quality("freeplanes",freeplanes);
    traj.set_quality("reseed",_reseed_ok);
    traj.set_quality("xtent",xtent);
    traj.set_quality("initialqP",_initialqP);
    traj.set_quality("fitted",_fitted);
    traj.set_quality("vertZ", vertZ);
    
    cout<<"traj ="<<i<<" failType ="<<_failType<<" fitted ="<<(int)traj.quality("fitted")<<endl;;  
    
    
    // std::cout<<"*********************************************"<<std::endl;
    //std::cout<<" ++++++ Trajectory "<<i<<"  fitted+++++++++++"<<std::endl;
    //std::cout<<"*********************************************\n"<<std::endl;
  }
  
  cout<<" I am ******************fitter end"<<endl; 
  
  //return _fitted;
  return true;///signifies fitter executed
 
}

//*************************************************************
void fitter::Reset() {
  //*************************************************************
  
  m.message("+++ Reset function +++",bhep::VERBOSE);
  
  //Reset trajectory 
  
   stc_tools::destroy(_meas);
   stc_tools::destroy(_trajs);
   _hadmeas.clear();
   _hadEng = 0;
   _failEvent = 0;///
   _pr_count = 0;///

   ///   
  _trajs.clear();
  _vPR_seed.clear();
     
}

//*************************************************************
bool fitter::FitTrajectory(Trajectory& traj, const State& seedState0, const int trajno) {
  //*************************************************************
  
  m.message("+++ FitTrajectory function ++++",bhep::VERBOSE);

  Trajectory traj1 = traj;
  Trajectory traj2 = traj;

  bool ok; 
  bool ok0, ok_quality;
  bool ok1 = false;

  /// fit the trajectory                
  ok0 = man().fitting_svc().fit(seedState0,traj);
  
  // Check the quality if the traj is fitted
  if(ok0) ok_quality = CheckQuality(traj); 
  
  ///refit the trajectory only when the quality is not good
  if (_refit && ok0 && !ok_quality){    
    State seedState1;
    ComputeSeedRefit(traj, seedState1);
    ok1 = man().fitting_svc().fit(seedState1,traj1);
  }

  ///check number of fitted nodes in traj
  int fitCheck =0;
  vector<Node*>::iterator nDIt;
  for (nDIt = traj.nodes().begin();nDIt!=traj.nodes().end();nDIt++)
    if ( (*nDIt)->status("fitted") )
      fitCheck++;
  // cout<<"fitCheck="<<fitCheck<<endl;
  
  
  ///check for reseeding 
  double low_fit_cut;
  if (_intType == 2)
    low_fit_cut = _lowFit2;//store.fetch_dstore("low_fit_cut2");
  else
    low_fit_cut = _lowFit1;//store.fetch_dstore("low_fit_cut0");



  //**********disallow backfit on cell auto tracks for now.
  
  ///reseed the trajectory
  if (_intType!= 5){   // not CA
    if (_intType != 2){ // not all planes are single occ  
      if (!ok0 || !ok1 || (double)fitCheck/(double)traj.nmeas() < low_fit_cut)
	_reseed_ok = ReseedTrajectory(traj2,trajno);
    } 
    else if ((double)fitCheck/(double)traj.nmeas() < low_fit_cut){      
      ///if traj contains all sinle occ planes
      _reseed_ok = ReseedTrajectory(traj2,trajno);      
    }
    else _pr_count++;    
  }
  // cout<<"pr_count="<<_pr_count<<endl; 
  // cout<<"***before :inside FitTrajectory, ok="<<ok<<" /reseed_called="<<_reseed_called<<" /reseed_ok="<< _reseed_ok<<endl;
  
  //if reseed successful
  ok=true; 
  
  if (_reseed_ok){    
    traj.clone(traj2);
    cout<<"traj2 =" << traj2<<endl;
  }
  else if (ok1){
    traj.clone(traj1);
    cout<<"traj1 =" << traj1<<endl;
  }
  else if (!ok0) 
    ok=false;

  cout<<"copied trajectory =" << traj<<endl;

  cout<<"***inside FitTrajectory, ok="<<ok<<endl;
  
  return ok;
  
}
  
//*************************************************************
bool fitter::ReseedTrajectory(Trajectory& traj,const int trajno){
  //*************************************************************
  m.message("****inside ReseedTrajectory *****************\n", bhep::VERBOSE);
  cout<<"****inside ReseedTrajectory *****************traj.nmeas = "<<traj.nmeas()<<"\n"<<endl;
  
  bool ok, ok1;
  _reseed_called = true;
  
  State backSeed = _vPR_seed[_pr_count];
  
  //Want to re-seed with diagonal matrix.
  HyperVector HV1 = backSeed.hv();//.keepDiagonalMatrix();
  HV1.keepDiagonalMatrix();
  backSeed.set_hv( HV1 );
    
  ///sort nodes in reverse order
  traj.sort_nodes(-1);

  // a copy of the track
  Trajectory traj1 = traj;
  
  ///fit the traj
  ok = man().fitting_svc().fit(backSeed,traj);
  
  /// compute seed for refitting and refit the traj
  if (ok && _refit){
    // Check the quality if the traj is fitted
    if (!CheckQuality(traj)){ 
      State seedState1;
      ComputeSeedRefit(traj,seedState1);
      ok1 = man().fitting_svc().fit(seedState1,traj1);
    }
  }

  if (ok1){
    traj.clone(traj1);
    ok=true;
  }

  // sort nodes back
  traj.sort_nodes(1);
    
  ///increase the count to set backseed from PR vector
  _pr_count++;
  
   return ok;
}

//*************************************************************
void fitter::ComputeSeedRefit(const Trajectory& traj, State& seedState) {
  //*************************************************************
  
  
  m.message("Going to calculate seed for refit...",bhep::VERBOSE);
  
  //--------- refit using a new seed --------//
  /// what is the differance in this new state than the earlier seed??	
  seedState = traj.state(traj.first_fitted_node());
  
  EVector v = seedState.vector();
  EMatrix C0 = seedState.matrix();
    
  ApplyCovarianceFactor(_facRef,C0);
  
  HyperVector HV(v,C0);
  HV.keepDiagonalMatrix();
  
  /// seedstate.set_hv( HV );
  seedState.set_hv( HV );
  
}

//*************************************************************
bool fitter::CheckQuality(const Trajectory& traj){
  //*************************************************************
    
  bool ok = true;
    
  if (traj.quality()>_chi2fit_max) ok=false;
       
  return ok;

}

//*************************************************************
bool fitter::CreateMeasurements(const bhep::particle& p) {
  //*************************************************************
 
  m.message("+++ CreateMeasurements function ++++",bhep::VERBOSE);
  
  Reset();
  
  bool ok = true;

  //string detect = store.fetch_sstore("detect");
  
  const vector<bhep::hit*> hits = p.hits( _detect ); 
  
  //Cluster or directly make measurements.
  if ( _doClust && hits.size() != 0 ){

    // Make clusters
    _clusters->execute( hits, _meas );}
  
  else {
    
    // Create a cluster of each hit (without clustering)
    for(size_t j=0; j< hits.size(); j++){
      
      //---------- create measurement ---------------//
      
      cluster* mnt = GetMeasurement(*hits[j]);
      
      _meas.push_back(mnt); 
      
      m.message("Measurement added:",*mnt,bhep::VVERBOSE);
    }//end of loop over hits
    
  }
  return ok;
}

//*************************************************************
bool fitter::CreateSingleTrajectory(Trajectory& traj) {
  //*************************************************************
 
  m.message("+++ CreateSingleTrajectory function ++++",bhep::VERBOSE);
  
  
  //--------- add measurements to trajectory --------//
     
  ///create the trajectory  
  std::vector<cluster*>::iterator it1;
  for (it1 = _meas.begin();it1 != _meas.end();it1++)
    traj.add_measurement( *(*it1) );
  
  m.message("Trajectory created:",traj,bhep::VVERBOSE);

  return true;
}


//*************************************************************
bool fitter::CheckValidTraj(const Trajectory& traj) {
  //*************************************************************
  //cout<<" +++++inside fitter:: CheckValidTraj func "<<endl;  
 

  //--------- Reject too many hits --------//
  
  if ((int)traj.nmeas() < _lowPass) { 
    _failType = 1;
    return false;
  }
  return true;
}

//*****************************************************************************
double fitf(Double_t *x,Double_t *par) { 
  //*****************************************************************************

  double z = x[0]; 
  double fitval = par[0]+par[1]*z+par[2]*z*z;

  return fitval ;

}

//*************************************************************
int fitter::GetQ(const Trajectory& traj){
  //*************************************************************
  
  if (model.compare("particle/helix")!=0) return 0;
  double q;
  
  /// 
  q = traj.state(traj.last_fitted_node()).vector()[dim-1];
  
  if (q<0) q=-1; else q=1;
  
  return (int) q;
  
}

//*************************************************************
cluster*  fitter::GetMeasurement(bhep::hit& hit){
  //*************************************************************
    
  m.message("+++ getMeasurement function ++++",bhep::VERBOSE);
    
  //---- generate a virtual plane to hold the hit ----//
    
  bhep::Point3D bhit_pos = hit.x(); 

  string meastype = geom.getMeasType();
  EMatrix cov = geom.getCov();
  //pnumber++;

  //----- generate repack hit from bhep one ----//
    
  EVector hit_pos(2,0);
  hit_pos[0]=bhit_pos[0];
  hit_pos[1]=bhit_pos[1];

  EVector meas_pos(3,0);
  meas_pos[0] = hit_pos[0];
  meas_pos[1] = hit_pos[1];
  meas_pos[2] = bhit_pos[2];
    

  cluster* me = new cluster();
  me->set_name(meastype);
  me->set_hv(HyperVector(hit_pos,cov));
  me->set_name("volume", "Detector");
  me->set_position( meas_pos );
  //Add the hit energy deposit as a key to the Measurement.
  const dict::Key Edep = "E_dep";
  const dict::Key EdepVal = bhep::to_string( hit.ddata("EnergyDep") );
  me->set_name(Edep, EdepVal);
  if (patternRec){
    const dict::Key motherP = "MotherParticle";
    //const dict::Key mothName = hit.mother_particle().name();
    const dict::Key mothName = hit.sdata( "true_moth" );
    me->set_name(motherP, mothName);
  }
  
  return me;
  
}

//*************************************************************
void fitter::Finalize() {
  //*************************************************************
   
  get_classifier().finalize();
  Reset();

  if ( _doClust ) delete _clusters;
  
}

//*************************************************************
void fitter::ComputeSeed(const Trajectory& traj, State& seedState, int firsthit) {
  //*************************************************************

  m.message("+++ computeSeed function ++++",bhep::VERBOSE);

  //use position slightly offset from first meas as seed 
  ///_lastIso is the total no of candidate muon hits inside muonTraj in free section before incremental filtering
  
  if ( (double)(traj.quality("lastIso"))/(double)traj.nmeas() > _min_iso_prop )
    firsthit = (int)traj.nmeas() - (int)(traj.quality("lastIso"));

  EVector v(6,0), v2(1,0);
  EMatrix C(6,6,0), C2(1,1,0);
    
  // take the position from the first hit
  v[0] = traj.nodes()[firsthit]->measurement().vector()[0];
  v[1] = traj.nodes()[firsthit]->measurement().vector()[1];
  v[2] = traj.nodes()[firsthit]->measurement().position()[2];   

  // Estime the momentum from range
  ComputeMomFromRange( traj, (int)traj.nmeas(), firsthit, v);

  double pSeed;
  double wFe = geom.get_Fe_prop();
  //Approximate p from plot of p vs. no. hits, then approx. de_dx from this.
  if (v[5] == 0) { //pSeed = (double)(0.060*traj.nmeas())*bhep::GeV;
    pSeed = (13300-11200*wFe) + (-128+190*wFe)*(double)traj.nmeas();
    v[5] = 1.0/pSeed;
  }
  
  // But use a larger covariance matrix
  // diagonal covariance matrix
  C[0][0] = C[1][1] = 9.*bhep::cm*bhep::cm;
  C[2][2] = EGeo::zero_cov()/2;
  C[3][3] = C[4][4] = 1.;
  C[5][5] = pow(v[5],2)*3;
  
  seedState.set_name(RP::particle_helix);
  seedState.set_name(RP::representation,RP::slopes_curv_z);
  
  v2[0] = 1;
  seedState.set_hv(RP::sense,HyperVector(v2,C2));
  seedState.set_hv(HyperVector(v,C));

  m.message("++ Seed estate after setSeed() in fitter:",seedState,bhep::VERBOSE);
}

//*************************************************************
void fitter::ApplyCovarianceFactor(double factor, EMatrix& C0){
  //*************************************************************
    
  //--- a large diagonal covariance matrix ---//
    
  C0 *= factor;
}

//*****************************************************************************
double fitf2(Double_t *x,Double_t *par) { 
  //*****************************************************************************

  double z = x[0];

  double fitval = par[0] + par[1]*z+par[2]*z*z+par[3]*z*z*z+par[4]*z*z*z*z;

  return fitval;
}

//*****************************************************************************
void fitter::ComputeMomFromParabola(const Trajectory& traj, int nplanes, int firsthit, EVector& V){
  //*****************************************************************************

  //Some catchers for pointless returns.
  int fitcatch;
  //
  int nfit, sign;
  int fitRange[3];
  const int fitpoints = nplanes - firsthit;
  
  double xpos[fitpoints], ypos[fitpoints], zpos[fitpoints];
  double upos[fitpoints], vpos[fitpoints];

  int pos = 0;

  EVector currentpos = EVector(3,0);
  EVector currentB   = EVector(3,0);
  EVector z = EVector(3,0);
  z[2] = 1;
  double Bmean=0;
  for( int ipoint=firsthit; ipoint < nplanes; ipoint ++ ){
    xpos[pos] = traj.measurement(ipoint).vector()[0];
    ypos[pos] = traj.measurement(ipoint).vector()[1];
    zpos[pos] = traj.measurement(ipoint).position()[2]
      - traj.measurement(firsthit).position()[2];
    currentpos[0] = traj.measurement(ipoint).vector()[0];
    currentpos[1] = traj.measurement(ipoint).vector()[1];
    currentpos[2] = 0.;
    currentB = geom.getBField(currentpos);
    upos[pos] = xpos[pos] > 0 ? asin(ypos[pos]/currentpos.norm())
      : -asin(ypos[pos]/currentpos.norm());
    vpos[pos] = dot(currentpos,crossprod(z,currentB))/currentB.norm();
    Bmean += currentB.norm();
    ++pos;
  }
  Bmean /= pos;
  Bmean /= bhep::tesla;
  
  if (fitpoints <= 15) { nfit = 1; fitRange[0] = fitpoints;}
  else if (fitpoints <= 40) { 
    nfit = 2;
    fitRange[0] = 15; fitRange[1] = (int)(0.7*fitpoints);
  }
  else if (fitpoints > 40) { 
    nfit = 3;
    fitRange[0] = 15; fitRange[1] = (int)(fitpoints/2); fitRange[2] = (int)(0.7*fitpoints);
  }
  for (int ifit = 0;ifit < nfit;ifit++) {
    TGraph *trajFitXZ = new TGraph(fitRange[ifit],zpos, xpos);
    TGraph *trajFitYZ = new TGraph(fitRange[ifit],zpos, ypos);
    TGraph *trajFitUZ = new TGraph(fitRange[ifit],zpos, upos);
    TGraph *trajFitVZ = new TGraph(fitRange[ifit],zpos, vpos);
    
    TF1 *func = new TF1("fit",fitf2,-3,3,3);
    func->SetParameters(0.,0.,0.001,0.0001,0.0001);
    func->SetParNames("a", "b", "c", "d", "e");
    
    TF1 *func2 = new TF1("fit2",fitf2,-3,3,3);
    func2->SetParameters(0.,0.,0.001,0.0001,0.0001);
    func2->SetParNames("f", "g", "h", "i", "j");

    TF1 *func3 = new TF1("fit3",fitf2,-3,3,3);
    func->SetParameters(0.,0.,0.001,0.0001,0.0001);
    func->SetParNames("a1", "b1", "c1", "d1", "e1");
    
    TF1 *func4 = new TF1("fit4",fitf2,-3,3,3);
    func2->SetParameters(0.,0.,0.001,0.0001,0.0001);
    func2->SetParNames("f1", "g1", "h1", "i1", "j1");

    fitcatch = trajFitXZ->Fit("fit", "QN");
    fitcatch = trajFitYZ->Fit("fit2", "QN");
    fitcatch = trajFitUZ->Fit("fit3", "QN");
    fitcatch = trajFitVZ->Fit("fit4", "QN");
    
    double b = func->GetParameter(1);
    double c = func->GetParameter(2);
    double g = func2->GetParameter(1);
    /*double f = func2->GetParameter(0);
      double a = func->GetParameter(0);
      double h = func2->GetParameter(2);  
      double a1 = func3->GetParameter(0);
      double b1 = func3->GetParameter(1);
      double c1 = func3->GetParameter(2);  
      double f1 = func4->GetParameter(0);*////
    double g1 = func4->GetParameter(1);
    double h1 = func4->GetParameter(2);  
    
    if (ifit == 0) {

      V[4] = g;   //func2->GetParameter(1);
      V[3] = b;

      if (h1!=0) {
	V[5] = 1./(-0.3*Bmean*pow((1+g1*g1),3./2.)/
		   (2*h1)*0.01);
	V[5] /= bhep::GeV;
	sign = (int)( V[5]/fabs( V[5] ));
      } else V[5] = 0;
    } else {
      if ((int)(-c/fabs(c)) == sign) {
	V[4] = g;
	V[3] = b;
	V[5] = 1/(-0.3*Bmean*pow((1+g1*g1),3./2.)/(2*h1)*0.01);
	V[5] /= bhep::GeV;
      } else break;
    }
    
    delete trajFitXZ;
    delete trajFitYZ;
    delete trajFitUZ;
    delete trajFitVZ;
  
    delete func;
    delete func2;
    delete func3;
    delete func4;
  }
  
  //std::cout<<"Momentum guess from polynomial fit: p/q = "<<1./V[5]<<std::endl;
}

//*****************************************************************************
void fitter::ComputeMomFromRange(const Trajectory& traj, int nplanes, int firsthit, EVector& V){
  //*****************************************************************************

  //Some catchers for pointless returns.
  int fitcatch;
  //
  /// int nfit;
  /// int fitRange[3];
  const int fitpoints = nplanes - firsthit;
  ///double meanchange = 0;
  double xpos[fitpoints], ypos[fitpoints], zpos[fitpoints];
  double upos[fitpoints];/// wpos[fitpoints];
  std::vector<EVector> dr;
  std::vector<EVector> B;
  bool isContained = true, cuspfound = false;



  double Xmax = geom.getPlaneX() - 1*bhep::cm;
  double Ymax = geom.getPlaneY() - 1*bhep::cm;
  /// double Zmax = geom.getPlaneZ() - 1*bhep::cm;


  //double dx[fitpoints-1], dy[fitpoints-1], dz[fitpoints-1];
  // double ax[fitpoints-2], ay[fitpoints-2], az[fitpoints-2];
  // double bx[fitpoints-2], by[fitpoints-2], bz[fitpoints-2];



  ///double ds0=0, ds1=0;
  double Bmean=0;
  double pathlength=0;
  int Npts=0;
  ///double initR = 0;
  double sumDR = 0;
  int minindex = nplanes - firsthit;
  ///double minR = 999999.9999;
  double pdR = 0.0;



  EVector Z = EVector(3,0); Z[2] = 1;
  for (int ipoint=firsthit;ipoint < nplanes;ipoint++){
    
    xpos[ipoint-firsthit] = traj.measurement(ipoint).vector()[0];
    ypos[ipoint-firsthit] = traj.measurement(ipoint).vector()[1];
    zpos[ipoint-firsthit] = traj.measurement(ipoint).position()[2]
      - traj.measurement(firsthit).position()[2];
    if(fabs(xpos[ipoint-firsthit]) > Xmax || fabs(ypos[ipoint-firsthit]) > Ymax)
      isContained = false;
    else if(fabs(ypos[ipoint-firsthit]) > 
	    (1 + tan(atan(1.)/2.)) * Xmax - fabs(xpos[ipoint-firsthit])) 
      isContained = false;
    EVector pos0 = EVector(3,0);
    pos0[0] = xpos[ipoint-firsthit];
    pos0[1] = ypos[ipoint-firsthit];
    pos0[2] = zpos[ipoint-firsthit];
    EVector B0 = geom.getBField(pos0);
    B.push_back(B0);
    Bmean += B0.norm();
    upos[ipoint-firsthit] = // sqrt(pos0[0]*pos0[0] + pos0[1]*pos0[1]);
      dot(pos0,crossprod(Z, B0))/crossprod(Z, B0).norm();
    //if(!cuspfound)
    //  if(ipoint == firsthit) initR = upos[ipoint-firsthit];
    //  else {
    //sumDR += initR - upos[ipoint-firsthit];
    //initR = upos[ipoint - firsthit];
    //  }
    Npts++;
    if ( ipoint > firsthit){
      EVector drtemp = EVector(3,0);
      drtemp[0] = xpos[ipoint-firsthit] - xpos[ipoint-firsthit-1];
      drtemp[1] = ypos[ipoint-firsthit] - ypos[ipoint-firsthit-1];
      drtemp[2] = zpos[ipoint-firsthit] - zpos[ipoint-firsthit-1];
      dr.push_back(drtemp);      
      pathlength +=  drtemp.norm();
      if ( ipoint > firsthit + 1 ) {
	int k = ipoint-firsthit-1;
	EVector dr0 = dr[k-1];
	EVector dr1 = dr[k];
	EVector ddr = dr1 + dr0;
	EVector Ddr = dr1 - dr0;
	EVector pos = EVector(3,0);
	pos[0] = xpos[k-1]; pos[1] = ypos[k-1]; pos[2] = zpos[k-1]; 
	EVector B = geom.getBField(pos);
	double dR = dot(ddr, crossprod(Z, B0))/ (crossprod(Z,B0).norm());
	double DR = dot(Ddr, crossprod(Z, B0))/ (crossprod(Z,B0).norm());
	if(pdR != 0.0){
	  if(!cuspfound && DR/fabs(DR) == pdR/fabs(pdR)){
	    // sumDR += fabs(dR) > 0.0 ? dR/fabs(dR):0.0;
	    sumDR += dR;
	    // pdR = dR;
	    pdR = dR;
	  }
	  else if(dR/fabs(dR) != pdR/fabs(pdR)){
	    // cuspfound = true;
	    minindex = ipoint - firsthit - 1;
	    pdR = dR;
	    // std::cout<<"At cusp, sumDR = "<<sumDR<<std::endl;
	  }
	}
	else if(!cuspfound && fabs(dR) > 0){
	  // sumDR += fabs(DR) > 0.0 ? DR/fabs(DR) : 0.0;
	  sumDR += dR;
	  pdR = dR;
	}
	/*
	  if(pdR != 0){
	  if(minR < upos[ipoint - firsthit - 1] &&
	  (xpos[ipoint-firsthit]/fabs(xpos[ipoint-firsthit]) != 
	  xpos[ipoint-firsthit-1]/fabs(xpos[ipoint-firsthit-1])) ||
	  (ypos[ipoint-firsthit]/fabs(ypos[ipoint-firsthit]) != 
	  ypos[ipoint-firsthit-1]/fabs(ypos[ipoint-firsthit-1]))
	  && !cuspfound){
	  minR = upos[ipoint - firsthit - 1];
	  minindex = ipoint - firsthit - 1;
	  cuspfound = true;
	  }
	  }*/
      }
    }
  }
  Bmean /=Npts;
  
  double wFe = geom.get_Fe_prop();
  double p = (wFe*(0.017143*bhep::GeV/bhep::cm * pathlength - 1.73144*bhep::GeV)
	      + (1- wFe)*(0.00277013*bhep::GeV/bhep::cm * pathlength + 1.095511*bhep::GeV));
  double meansign = 1;
  if(sumDR != 0) {
    //std::cout<<"sumDR = "<<sumDR<<std::endl;
    meansign = sumDR/fabs(sumDR);
  }
  
  ///double planesep  = fabs(zpos[1] - zpos[0]);



  // Assume that the magnetic field does not change very much over 1 metre
  // (terrible assumption by the way)
  const int sample = minindex < 20 ? (const int)minindex: 20;
  
  V[3] = dr.at(0)[0]/dr.at(0)[2];
  V[4] = dr.at(0)[1]/dr.at(0)[2];
  if(isContained && p != 0)
    V[5] = meansign/fabs(p);
  else{


    // meansign = 0;
    // Consider a fit to a subset of points at the begining of the track
    //for(int j=0; j<fitpoints-2; j++){
    TGraph* localcurveUW = new TGraph(sample, zpos, upos);
    TF1 *func = new TF1("fit",fitf2,-30,30,4);
    func->SetParameters(0.,0.,0.001,0.0001,0.0001);
    func->SetParNames("a", "b", "c", "d", "e");
    fitcatch = localcurveUW->Fit("fit", "QN");
    double b = func->GetParameter(1);
    double c = func->GetParameter(2);
    
    p = 300.0 * B[1].norm() * pow(1. + b*b,3./2.) /2./c;
    //double wt = TMath::Gaus(fabs(pt), p, 0.25*p, true);
    // std::cout<<pt<<std::endl;
    // meansign += wt * pt/fabs(pt);
    
    delete localcurveUW;      
    delete func;
    if(p!=0){
      meansign = p/fabs(p);
      V[5] = 1./p;
    }
  
  }
  int sign = 1;
  if(meansign==meansign){
    if(meansign!=0)
      sign = int(meansign/fabs(meansign));
    else
      sign = 0;
  }
  else
    sign = 1;




  // std::cout<<"Pathlength is "<<pathlength // <<" or "<<pathlength0
  //	   <<" with charge "<<meansign<<std::endl;
  
  _initialqP = V[5];
  //cout<<"_initialqP ="<<_initialqP <<endl;
  
}

//*****************************************************************************
void fitter::ReadParam(){
  //*****************************************************************************
    
  m.message("+++ ReadParam function of fitter ++++",bhep::VERBOSE);
        
  model = store.find_sstore("model");//"particle/helix"; 
  dim=6; // ??????
    

  if ( store.find_istore("refit") )
    _refit=store.fetch_istore("refit");
  else _refit=false;

  if ( store.find_istore("patRec") )
    patternRec=store.fetch_istore("patRec");
  else patternRec=false;

  if ( store.find_istore("do_clust") )
    _doClust = store.fetch_istore("do_clust");
  else _doClust = false;

  _facRef = store.fetch_dstore("facRef");

  _min_seed_hits = store.fetch_istore("min_seed_hits");
  _min_iso_prop = store.fetch_dstore("min_iso_prop");

  _chi2fit_max = store.fetch_dstore("chi2fit_max");
    
  _X0 = store.fetch_dstore("x0Fe") * bhep::mm;
  //_tolerance = store.fetch_dstore("pos_res") * bhep::cm;
  _highPass = store.fetch_istore("high_Pass_hits");
  _lowPass = store.fetch_istore("low_Pass_hits");
  _lowFit1 = store.fetch_dstore("low_fit_cut0");
  _lowFit2 = store.fetch_dstore("low_fit_cut2");

      
}


