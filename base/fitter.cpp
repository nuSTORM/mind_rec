
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
 
  //reset();
  //if ( _doClust ) delete _clusters;
  
}

//*************************************************************
//bool fitter::initialize(const bhep::sstore& run_store) {
void fitter::initialize() {
//*************************************************************
    
  m.message("+++ fitter init  function ++++",bhep::VERBOSE);
  
  _hadunit = EVector(3,0);
  _hadunit[2] = 1.;
  _detect = store.fetch_sstore("detect");
  // read parameters
  readParam();
  
  // initialize geometry
  
  geom.init(store, level);

  //Instantiate recpack manager.
  MINDfitman::instance().set_man_parameters( store, geom.setup() );
  
  if (X0 == 0){
    man().model_svc().enable_noiser(model, RP::ms, false);
  }

  //If required make the clustering object.
  if ( _doClust )
    _clusters = new hit_clusterer( store );
  //
  //get_classifier().initialize( store, level, geom.setup(), geom.get_Fe_prop() );
  get_classifier().initialize( store, level, geom.get_Fe_prop() );

  m.message("+++ End of init function ++++",bhep::VERBOSE);
  
  //return true;
}

//*************************************************************
bool fitter::execute(bhep::particle& part,int evNo){
//*************************************************************
    
  m.message("+++ fitter execute function ++++",bhep::VERBOSE);
  
  bool ok;//, checker; 
  bool fitted=false;
  
  _failType = 0; //set to 'success' before run to avoid faults in value.
  ok = readTrajectory(part);
  
  reseed_ok = false;

  if ( ok ) computeSeed();
  
  if (ok) {
    
    fitted = fitTrajectory(seedstate);
    
    // if ( fitted )
//       checker = fitHadrons();
//     if (!checker) std::cout << "Had fit fail" << std::endl;
    
    double length;
    if ( reseed_ok )
      ok = man().matching_svc().compute_length(_traj2, length);
    else ok = man().matching_svc().compute_length(_traj, length);
    
    // bool noprot = true;
    if (reseed_ok && fitted)
      noprot = test_de_dx(_traj2.state(_traj2.first_fitted_node()));
    else if(fitted) noprot = test_de_dx(_traj.state(_traj.first_fitted_node()));
    

    // if (fitted) m.message("++ Particle fitted",bhep::VERBOSE);
    //     else m.message("++ Particle not fitted !!!",bhep::VERBOSE);
    
    // if (!fitted)
    //       m.message("++Failed fit trajectory++",bhep::NORMAL);
  }
  // else m.message("++ Particle lies outside no. hit cuts!",bhep::VERBOSE);
  
  if (fitted) 
    if (_failType!=3) _failType = 0;
  
  return fitted;
  
}

//*************************************************************
void fitter::reset() {
//*************************************************************
  
  m.message("+++ reset function +++",bhep::VERBOSE);
  
  //reset trajectory 
  
  _traj.reset();
  _traj1.reset();
  _traj2.reset();
  stc_tools::destroy(_meas);
  _hadmeas.clear();
  _hadEng = 0;
  
  
}

// //*************************************************************
// void fitter::calculate_len_mom(double len, double *mom){
// //*************************************************************

//   double weight = geom.get_Fe_prop();

//   //parameters form fit to range data for Scint and Fe.

//   mom[0] = exp( ( log( len ) - 8.35 + 1.72 * weight ) / ( 0.981 - 0.023 * weight ) );
//   mom[0] *= bhep::GeV;
//   //All considered?? 5% from range fit and 1cm error on length?
//   mom[1] = sqrt( pow( 0.05, 2) + pow( 10.0/len, 2) ) / ( 0.981 - 0.023 * weight );
//   mom[1] *= mom[0]; 
  
// }

//*************************************************************
bool fitter::fitTrajectory(State seed) {
//*************************************************************
    
  m.message("+++ fitTrajectory function ++++",bhep::VERBOSE);
  
  bool ok = man().fitting_svc().fit(seed,_traj);
  
  if (ok && refit){
    
    ok = checkQuality(); 
    if (ok) {
      
      m.message("Going to refit...",bhep::VERBOSE);
      
      //--------- refit using a new seed --------//	
      State newstate = _traj.state(_traj.first_fitted_node());

      EVector v = newstate.vector();
      EMatrix C0 = newstate.matrix();
      
    
      set_de_dx( fabs(1./v[5])/bhep::GeV );
      
      EMatrix C = setSeedCov(C0,facRef);
      
      HyperVector HV(v,C);
      HV.keepDiagonalMatrix();
      //seedstate.set_hv(HyperVector(v,C)); 
      seedstate.set_hv( HV );
      
      ok = man().fitting_svc().fit(seedstate,_traj);
      
    }
    
  }
  
  int fitCheck = 0;
  vector<Node*>::iterator nDIt;
  for (nDIt = _traj.nodes().begin();nDIt!=_traj.nodes().end();nDIt++)
    if ( (*nDIt)->status("fitted") )
      fitCheck++;
  
  double low_fit_cut;
  if (get_classifier().get_int_type() == 2)
    low_fit_cut = lowFit2;//store.fetch_dstore("low_fit_cut2");
  else
    low_fit_cut = lowFit1;//store.fetch_dstore("low_fit_cut0");
  //disallow backfit on cell auto tracks for now.
  
  if (get_classifier().get_int_type() != 5){
    if (get_classifier().get_int_type() != 2){
      if (!ok || (double)fitCheck/(double)_traj.nmeas() < low_fit_cut)
	reseed_ok = reseed_traj();
    } else if ((double)fitCheck/(double)_traj.nmeas() < low_fit_cut)
      reseed_ok = reseed_traj();
  }
  
  if (ok && !reseed_ok) ok = checkQuality();
  else if ( reseed_ok ) ok = true;
  
  return ok;
  
}

//*************************************************************
bool fitter::reseed_traj(){
//*************************************************************
  bool ok;
    
  State backSeed = get_classifier().get_patRec_seed();
  //Want to re-seed with diagonal matrix.
  HyperVector HV1 = backSeed.hv();//.keepDiagonalMatrix();
  HV1.keepDiagonalMatrix();
  backSeed.set_hv( HV1 );
  
  vector<Node*>::iterator nIt;
  for (nIt = _traj.nodes().end()-1;nIt >= _traj.nodes().begin();nIt--)
    _traj2.add_measurement( (*nIt)->measurement() );
  
  ok = man().fitting_svc().fit(backSeed,_traj2);
  
  if (ok && refit){
     
    if (_traj2.quality() <= chi2fit_max) {

      m.message("Going to refit...",bhep::VERBOSE);
      
      //--------- refit using a new seed --------//	
      State newstate = _traj2.state(_traj2.first_fitted_node());

      EVector v = newstate.vector();
      EMatrix C0 = newstate.matrix();
      
      set_de_dx( fabs(1./v[5])/bhep::GeV );
      
      EMatrix C = setSeedCov(C0,facRef);

      HyperVector HV(v,C);
      HV.keepDiagonalMatrix();
      //seedstate.set_hv(HyperVector(v,C)); 
      seedstate.set_hv( HV );

      ok = man().fitting_svc().fit(seedstate,_traj2);
    }
    
  }
  
  return ok;
}

// //*************************************************************
// bool fitter::fitHadrons(){
// //*************************************************************
  
//   size_t nhadhits = _hadmeas.size();

//   measurement_vector::iterator engIt;
//   const dict::Key Edep = "E_dep";

//   for (engIt = _hadmeas.begin();engIt != _hadmeas.end();engIt++)
//     _hadEng += bhep::double_from_string( (*engIt)->name( Edep ) );

//   _hadEng = eng_scale( _hadEng );

//   if (nhadhits<2) {
//     _hadunit[0] = 0; _hadunit[1] = 0; return false;}
//   const int nplanes = ( (int)_hadmeas[0]->position()[2]
// 			- (int)_hadmeas[nhadhits-1]->position()[2] + 50 )/50;
  
//   if (nplanes<2) {
//     _hadunit[0] = 0; _hadunit[1] = 0; return false;}
//   double x[nplanes], y[nplanes], z[nplanes];
//   double EngPlane, EngHit, testZ, curZ;
  
//   size_t hits_used = 0, imeas = 0;
//   int ientry = nplanes-1;

//   do {  
    
//     int count = 0;
//     EngPlane = 0;
//     EngHit = bhep::double_from_string( _hadmeas[imeas]->name( Edep ) )*bhep::GeV;
//     x[ientry] = _hadmeas[imeas]->vector()[0] * EngHit;
//     y[ientry] = _hadmeas[imeas]->vector()[1] * EngHit;
//     z[ientry] = _hadmeas[imeas]->position()[2];
//     testZ = z[ientry];
//     hits_used++;
//     count++;
//     EngPlane += EngHit;
    
//     for (size_t i=hits_used;i < nhadhits;i++){
//       curZ = _hadmeas[i]->position()[2];
//       if (curZ >= testZ-_tolerance){
// 	EngHit = bhep::double_from_string( _hadmeas[i]->name( Edep ) )*bhep::GeV;
// 	x[ientry] += _hadmeas[i]->vector()[0] * EngHit;
// 	y[ientry] += _hadmeas[i]->vector()[1] * EngHit;
// 	z[ientry] += curZ;
// 	count++;
// 	hits_used++;
// 	EngPlane += EngHit;

//       } else break;
//     }
    
//     x[ientry] /= EngPlane; y[ientry] /= EngPlane; 
//     z[ientry] /= (double)count;
    
//     ientry--;
//     imeas+=count;
    
//   } while (hits_used != nhadhits);
  
//   TGraph *gr1 = new TGraph(nplanes, z, x);
//   TGraph *gr2 = new TGraph(nplanes, z, y);
  
//   TF1 *fitfunc = new TF1("fitfun","[0]+[1]*x",-3,3);
  
//   gr1->Fit("fitfun","QN");
//   _hadunit[0] = fitfunc->GetParameter(1);

//   gr2->Fit("fitfun","QN");
//   _hadunit[1] = fitfunc->GetParameter(1);
  
//   _hadunit /= _hadunit.norm();

//   delete gr1;
//   delete gr2;
//   delete fitfunc;

//   return true;
// }

// //*************************************************************
// double fitter::eng_scale(double visEng){
// //*************************************************************

//   double factor1, factor2, scaledEng;
//   //Equations/error calc. in log. root from parafit at mo.
//   if ( visEng==0 ) return -99;
//   else if ( visEng > 0.03263 ) visEng = 0.03263;

//   factor1 = (-0.0024 + sqrt( pow(-0.0024, 2) - 4*(0.00288-visEng)*-0.0000484))
//     /(2*-0.0000484);
//   factor2 = (-0.0024 - sqrt( pow(-0.0024, 2) - 4*(0.00288-visEng)*-0.0000484))
//     /(2*-0.0000484);
  
//   if ( factor1 < 0 )
//     scaledEng = factor2;
//   else if ( factor2 < 0 )
//     scaledEng = factor1;
//   else scaledEng = TMath::Min(factor1, factor2);
  
//   return scaledEng * bhep::GeV;
// }

//*************************************************************
bool fitter::checkQuality(){
//*************************************************************

    
    bool ok = true;
    
    if (getChi2()>chi2fit_max) ok=false;
       
    return ok;

}

//*************************************************************
bool fitter::readTrajectory(const bhep::particle& part){
//*************************************************************
    
  m.message("+++ readTrajectory function ++++",bhep::VERBOSE);
  
  //---------- read trajectory, i.e, particle hits ---------//
  
  //generate rec trajectory
  bool ok = recTrajectory(part);
  
  if (patternRec && ok){
    
    ok = get_classifier().execute( _meas, _traj, _hadmeas);
    
    _traj.sort_nodes(1);
    sort( _hadmeas.begin(), _hadmeas.end(), reverseSorter() );
    
    if (!ok){
      _failType = get_classifier().get_fail_type();
      //patFail++;
    } else {
      //SetFit mode to manager.
      MINDfitman::instance().fit_mode();
    }
  }
  // // ******HARDWIRE FAIL******** only interested in likelihood info at the mo
//   ok = false;
//   _failType = 5;
  
  // Check that the 'muon' can be fitted.
  // Is the lowest z hit in fiducial volume?
  // Are there enough hits? Are there too many?
  
  if (ok)
    ok = check_valid_traj();
  
  return ok;
}


//*************************************************************
bool fitter::recTrajectory(const bhep::particle& p) {
//*************************************************************

  m.message("+++ recTrajectory function ++++",bhep::VERBOSE);
 
    reset();

    //string detect = store.fetch_sstore("detect");

    const vector<bhep::hit*> hits = p.hits( _detect ); 
    
    //Cluster or directly make measurements.
    if ( _doClust && hits.size() != 0 )
      _clusters->execute( hits, _meas );
    else {
 
      for(size_t j=0; j< hits.size(); j++){
        	        
	//---------- create measurement ---------------//
	
	cluster* mnt = getMeasurement(*hits[j]);
	
	_meas.push_back(mnt); 
	
	m.message("Measurement added:",*mnt,bhep::VVERBOSE);
	
      }//end of loop over hits

    }
    
    //--------- add measurements to trajectory --------//
    //Sort in increasing z here when classifier up and running.!!!
    if (patternRec) {

      if ((int)_meas.size() < min_seed_hits) {
      _failType = 7; return false;
      }
      
      sort( _meas.begin(), _meas.end(), forwardSorter() );
      
    } else {
      //_traj.add_measurements(_meas);
      std::vector<cluster*>::iterator it1;
      for (it1 = _meas.begin();it1 != _meas.end();it1++)
	_traj.add_measurement( *(*it1) );

      m.message("Trajectory created:",_traj,bhep::VVERBOSE);
    }
    
    return true;
}

//*************************************************************
bool fitter::check_valid_traj() {
//*************************************************************

  // double zCut = store.fetch_dstore("z_cut");
//   double xCut = store.fetch_dstore("x_cut");
//   double yCut = store.fetch_dstore("y_cut");

  //--------- Reject too many hits --------//
  // int highPass = store.fetch_istore("high_Pass_hits");
//   int lowPass = store.fetch_istore("low_Pass_hits");

  // if ((int)_traj.nmeas() > highPass) { 
//     _failType = 2; 
//     return false; }

  if ((int)_traj.nmeas() < lowPass) { 
    _failType = 1;
    return false;
  }
  
  //---- Reject if initial meas outside fid. Vol ----//
  // if (_traj.nodes()[0]->measurement().position()[2] > geom.getPlaneZ()/2-zCut*cm
//       || fabs(_traj.nodes()[0]->measurement().vector()[0]) > geom.getPlaneX()/2-xCut*cm
//       || fabs(_traj.nodes()[0]->measurement().vector()[1]) > geom.getPlaneY()/2-yCut*cm)
//     {
//     _failType = 3;}
    
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
int fitter::getQ(){
//*************************************************************
  
  if (model.compare("particle/helix")!=0) return 0;
  double q;
  
  if ( reseed_ok )
    q = _traj2.state(_traj2.last_fitted_node()).vector()[dim-1];
  else q = _traj.state(_traj.last_fitted_node()).vector()[dim-1];
  
  if (q<0) q=-1; else q=1;

  return (int) q;

}


//*************************************************************
cluster*  fitter::getMeasurement(bhep::hit& hit){
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
void fitter::finalize() {
//*************************************************************
   
  get_classifier().finalize();
  reset();

  if ( _doClust ) delete _clusters;
  
  //return true;
}


//*************************************************************
void fitter::computeSeed(int firsthit) {
//*************************************************************
    
  m.message("+++ computeSeed function ++++",bhep::VERBOSE);
    
    //use position slightly offset from first meas as seed 

    EVector v(3,0); 
    if ( (double)get_classifier().get_last_iso()/(double)_traj.nmeas() > min_iso_prop )
      firsthit = (int)_traj.nmeas() - get_classifier().get_last_iso();

    v[0] = _traj.nodes()[firsthit]->measurement().vector()[0];
    v[1] = _traj.nodes()[firsthit]->measurement().vector()[1];
    v[2] = _traj.nodes()[firsthit]->measurement().position()[2];   

    setSeed(v, firsthit);
    
    m.message("++ Seed estate:",seedstate,bhep::VERBOSE);

}

//*************************************************************
void fitter::setSeed(EVector r, int firsthit){
//*************************************************************
  
  m.message("+++ setSeed function ++++",bhep::VERBOSE);
  
  EVector v(6,0), v2(1,0);
  EMatrix C(6,6,0), C2(1,1,0);

  v[0]=r[0];
  v[1]=r[1];
  v[2]=r[2];
  
  mom_from_range( (int)_traj.nmeas(), firsthit, v);
  // 
  
  std::cout<<"Momentum guess from range: p/q = "<<1./v[5]<<std::endl;

  double pSeed;
  double wFe = geom.get_Fe_prop();
  //Approximate p from plot of p vs. no. hits, then approx. de_dx from this.
  if (v[5] == 0) { //pSeed = (double)(0.060*_traj.nmeas())*bhep::GeV;
    pSeed = (13300-11200*wFe) + (-128+190*wFe)*(double)_traj.nmeas();
    v[5] = 1.0/pSeed;
  }
  
  set_de_dx( fabs(1./v[5])/bhep::GeV );

  m.message("reset energy loss to approx value",bhep::VERBOSE);

  // But use a larger covariance matrix
  // diagonal covariance matrix
  C[0][0] = C[1][1] = 9.*bhep::cm*bhep::cm;
  C[2][2] = EGeo::zero_cov()/2;
  C[3][3] = C[4][4] = 1.;
  C[5][5] = pow(v[5],2)*3;
  
  seedstate.set_name(RP::particle_helix);
  seedstate.set_name(RP::representation,RP::slopes_curv_z);
  
  v2[0] = 1;
  seedstate.set_hv(RP::sense,HyperVector(v2,C2));
  seedstate.set_hv(HyperVector(v,C));
  

}

//*************************************************************
EMatrix fitter::setSeedCov(EMatrix C0, double factor){
//*************************************************************
    
    //--- a large diagonal covariance matrix ---//
    
  EMatrix c = factor*C0;


    return c;
}

//*****************************************************************************
double fitf2(Double_t *x,Double_t *par) { 
//*****************************************************************************

  double z = x[0];

  double fitval = par[0] + par[1]*z+par[2]*z*z+par[3]*z*z*z+par[4]*z*z*z*z;

  return fitval;
}

//*****************************************************************************
void fitter::mom_from_parabola(int nplanes, int firsthit, EVector& V){
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
  for( int ipoint=firsthit; ipoint < nplanes; ipoint ++ )
    {
      xpos[pos] = _traj.measurement(ipoint).vector()[0];
      ypos[pos] = _traj.measurement(ipoint).vector()[1];
      zpos[pos] = _traj.measurement(ipoint).position()[2]
	- _traj.measurement(firsthit).position()[2];
      currentpos[0] = _traj.measurement(ipoint).vector()[0];
      currentpos[1] = _traj.measurement(ipoint).vector()[1];
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
    double a = func->GetParameter(0);
    double b = func->GetParameter(1);
    double c = func->GetParameter(2);  
    double f = func2->GetParameter(0);
    double g = func2->GetParameter(1);
    double h = func2->GetParameter(2);  
    double a1 = func3->GetParameter(0);
    double b1 = func3->GetParameter(1);
    double c1 = func3->GetParameter(2);  
    double f1 = func4->GetParameter(0);
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
  
  std::cout<<"Momentum guess from polynomial fit: p/q = "<<1./V[5]<<std::endl;
}


//*****************************************************************************
void fitter::mom_from_range(int nplanes, int firsthit, EVector& V){
//*****************************************************************************

//Some catchers for pointless returns.
  int fitcatch;
  //
  int nfit;
  // int fitRange[3];
  const int fitpoints = nplanes - firsthit;
  double meanchange = 0;
  double xpos[fitpoints], ypos[fitpoints], zpos[fitpoints];
  double upos[fitpoints], rpos[fitpoints];
  std::vector<EVector> dr;
  std::vector<EVector> B;
  bool isContained = true, cuspfound = false;
  double Xmax = geom.getPlaneX() - 1*bhep::cm;
  double Ymax = geom.getPlaneY() - 1*bhep::cm;
  double Zmax = geom.getPlaneZ() - 1*bhep::cm;
  //double dx[fitpoints-1], dy[fitpoints-1], dz[fitpoints-1];
  // double ax[fitpoints-2], ay[fitpoints-2], az[fitpoints-2];
  // double bx[fitpoints-2], by[fitpoints-2], bz[fitpoints-2];
  double ds0=0, ds1=0;
  double Bmean=0;
  double pathlength=0;
  int Npts=0;
  double initR = 0;
  double sumDR = 0;
  int minindex = nplanes - firsthit;
  double min_r = 999999.9999;
  double pdR = 0.0;
  EVector Z = EVector(3,0); Z[2] = 1;
  bool bounce=false;
  for (int ipoint=firsthit;ipoint < nplanes;ipoint++){
    
    xpos[ipoint-firsthit] = _traj.measurement(ipoint).vector()[0];
    ypos[ipoint-firsthit] = _traj.measurement(ipoint).vector()[1];
    zpos[ipoint-firsthit] = _traj.measurement(ipoint).position()[2]
      - _traj.measurement(firsthit).position()[2];
    rpos[ipoint-firsthit] = sqrt(pow(xpos[ipoint-firsthit],2) + pow(ypos[ipoint-firsthit],2));
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
    // std::cout<<pos0<<B0<<std::endl;
    B.push_back(B0);
    Bmean += B0.norm();
    upos[ipoint-firsthit] = dot(pos0,crossprod(Z, B0))/B0.norm();
    // int usign = upos[ipoint-firsthit]/fabs(upos[ipoint-firsthit]); 
    // Find the minimum radial position.
    if ( min_r > rpos[ipoint-firsthit] && !bounce){
      min_r = rpos[ipoint-firsthit];
      minindex = ipoint-firsthit;
    }
    else if(minindex > 0 && min_r < rpos[ipoint-firsthit] &&
	    !bounce) bounce = true;
    
    if ( ipoint > firsthit){
      EVector drtemp = EVector(3,0);
      drtemp[0] = xpos[ipoint-firsthit] - xpos[ipoint-firsthit-1];
      drtemp[1] = ypos[ipoint-firsthit] - ypos[ipoint-firsthit-1];
      drtemp[2] = zpos[ipoint-firsthit] - zpos[ipoint-firsthit-1];
      dr.push_back(drtemp);      
      pathlength +=  drtemp.norm();
    }
    Npts++;
  }
  // Now to determine the charge, get the first upos and the last
  // upos.  
  // If the last rpos is greater than the first rpos, it is
  // likely that the magnetic field focusses the charge, and it is
  // possible that there is a "bounce" in the track.
  
  double firstu = upos[1];
  double lastu  = upos[fitpoints-1];
  double firstz = zpos[1];
  double lastz  = zpos[fitpoints-1];
  int umax = Npts;
  
  if(rpos[1] > rpos[Npts-1] && bounce ){
    lastu = upos[minindex];
    lastz = zpos[minindex];
    umax = minindex;
  }
  // Now define a line in u-z space 
  double uslope = (lastu - firstu)/(lastz - firstz);
  double uintercept = (-firstz)/(lastz-firstz) + firstu;
  std::cout<<firstz<<"\t"<<lastz<<"\t"<<firstu<<"\t"<<lastu<<"\t"<<uslope<<"\t"<<uintercept<<std::endl;
  for (int iu = 0; iu < umax; iu++){
    double dr = upos[iu] - ( uslope * zpos[iu] + uintercept ); 
    //if(dr != dr) 
    // else 
    sumDR -= dr;
  }
  /*
      if ( ipoint > firsthit + 1 ) {
	int k = ipoint-firsthit-1;
	EVector dr0 = dr[k-1];
	EVector dr1 = dr[k];
	EVector ddr = dr1 + dr0;
	EVector Ddr = dr1 - dr0;
	EVector pos = EVector(3,0);
	pos[0] = xpos[k-1]; pos[1] = ypos[k-1]; pos[2] = zpos[k-1]; 
	EVector B = geom.getBField(pos);
	double dR = dot(Ddr, crossprod(Z, B0));
	// double DR = dot(Ddr, crossprod(Z, B0));
	sumDR += fabs(dR) > 0.0 ? dR/fabs(dR):0.0;
	// pdR = dR;
      }
      
    }
  }
  */
  Bmean /=Npts;
  
  double p = get_classifier().RangeMomentum(pathlength);
  double meansign = 1;
  std::cout<<"sumDR = "<<sumDR<<", Mean B field is "<<Bmean<<std::endl;
  if(sumDR != 0 && sumDR == sumDR) {
    meansign = sumDR/fabs(sumDR);
  }
  
  double planesep  = fabs(zpos[1] - zpos[0]);
  // Assume that the magnetic field does not change very much over 1 metre
  // (terrible assumption by the way)
  const int sample = minindex < 20 ? (const int)minindex: 20;
  
  V[3] = dr.at(1)[0]/dr.at(1)[2];
  V[4] = dr.at(1)[1]/dr.at(1)[2];
  // V[3] = (xpos[1]  - xpos[0])/(zpos[1] - zpos[0]);
  // V[4] = (ypos[1] - ypos[0])/(zpos[1] - zpos[0]);
  if(isContained && p != 0) // && fabs(meansign)==1.0)
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
    
    double p1 = 300.0 * B[1].norm() * pow(1. + b*b,3./2.) /2./c;
    //double wt = TMath::Gaus(fabs(pt), p, 0.25*p, true);
    // std::cout<<pt<<std::endl;
    // meansign += wt * pt/fabs(pt);
    
    delete localcurveUW;      
    delete func;
    if(p1!=0){
      meansign = p1/fabs(p1);
    }
    
  }
  V[5] = meansign/p;
  std::cout<<"Pathlength is "<<pathlength <<" for p="<<p
	   <<" GeV/c, with charge "<<meansign<<std::endl;
    //<<" and acceleration in bending plane "
    //   <<meanacceleration/count<<" from "
    //   <<count<<" points."<<std::endl;

  //std::cout<<"Weight of iron is "<<wFe<<std::endl;
  //std::cout<<"Numbers taken from interpolations of Range tables between 1 and 80 bhep::GeV. \n";
    
  // std::cout<<"Momentum is "<<p<<std::endl;
  _initialqP = V[5];
  // std::cout<<"Momentum pull is "<<sign/p<<std::endl;
  
  // V[5] = (sign)/p;
  // return 1./p;
}

//*****************************************************************************
void fitter::set_de_dx(double mom){
//*****************************************************************************

  double weight = geom.get_Fe_prop();
  
  double de_dx = -( 12.37 * weight * pow( mom, 0.099)
		    + (1 - weight) * 2.16 * pow( mom, 0.075) );
  de_dx *= bhep::MeV/bhep::cm;
  
  geom.setDeDx(de_dx);

}


//*****************************************************************************
bool fitter::test_de_dx(State baseseed){
//*****************************************************************************
// Test whether the track produces a better fit if other particles are assumed

  m.message("+++ test_de_dx function +++",bhep::VERBOSE);
  
  double basechi2 = _traj.quality();
  State proseed = baseseed;
  proseed.set_name(RP::PID,"Proton");
  
  
  vector<Node*>::iterator nIt;
  for (nIt = _traj.nodes().begin();nIt < _traj.nodes().end();nIt++)
    _traj1.add_measurement( (*nIt)->measurement() );


  bool ok = man().fitting_svc().fit(proseed,_traj1);
  /* 
  if (ok && refit){
    
    ok = checkQuality(); 
    if (ok) {
      
      m.message("Going to refit...",bhep::VERBOSE);
      
      //--------- refit using a new seed --------//	
      State newstate = _traj1.state(_traj1.first_fitted_node());

      EVector v = newstate.vector();
      EMatrix C0 = newstate.matrix();
      
    
      set_de_dx( fabs(1./v[5])/bhep::GeV );
      
      EMatrix C = setSeedCov(C0,facRef);
      
      HyperVector HV(v,C);
      HV.keepDiagonalMatrix();
      //seedstate.set_hv(HyperVector(v,C)); 
      seedstate.set_hv( HV );
      
      ok = man().fitting_svc().fit(seedstate,_traj1);
      
      }
      
      }
      
  int fitCheck = 0;
  vector<Node*>::iterator nDIt;
  for (nDIt = _traj.nodes().begin();nDIt!=_traj.nodes().end();nDIt++)
    if ( (*nDIt)->status("fitted") )
      fitCheck++;
  
  double low_fit_cut;
  if (get_classifier().get_int_type() == 2)
    low_fit_cut = lowFit2;//store.fetch_dstore("low_fit_cut2");
  else
    low_fit_cut = lowFit1;//store.fetch_dstore("low_fit_cut0");
  //disallow backfit on cell auto tracks for now.
  
  if (get_classifier().get_int_type() != 5){
    if (get_classifier().get_int_type() != 2){
      if (!ok || (double)fitCheck/(double)_traj.nmeas() < low_fit_cut)
	reseed_ok = reseed_traj();
    } else if ((double)fitCheck/(double)_traj.nmeas() < low_fit_cut)
      reseed_ok = reseed_traj();
  }
  */
  if (ok) ok = basechi2 < _traj1.quality();
  _chi2diff = _traj1.quality() - basechi2;

  return ok;
  
}  

//*****************************************************************************
void fitter::readParam(){
//*****************************************************************************
    
  m.message("+++ readParam function of fitter ++++",bhep::VERBOSE);
        
    model = store.find_sstore("model");//"particle/helix"; 
    dim=6; // ??????
    

    if ( store.find_istore("refit") )
      refit=store.fetch_istore("refit");
    else refit=false;

    if ( store.find_istore("patRec") )
      patternRec=store.fetch_istore("patRec");
    else patternRec=false;

    if ( store.find_istore("do_clust") )
      _doClust = store.fetch_istore("do_clust");
    else _doClust = false;

    facRef = store.fetch_dstore("facRef");

    min_seed_hits = store.fetch_istore("min_seed_hits");
    min_iso_prop = store.fetch_dstore("min_iso_prop");

    chi2fit_max = store.fetch_dstore("chi2fit_max");
    
    X0 = store.fetch_dstore("x0Fe") * bhep::mm;
    //_tolerance = store.fetch_dstore("pos_res") * bhep::cm;
    highPass = store.fetch_istore("high_Pass_hits");
    lowPass = store.fetch_istore("low_Pass_hits");
    lowFit1 = store.fetch_dstore("low_fit_cut0");
    lowFit2 = store.fetch_dstore("low_fit_cut2");
    
}




