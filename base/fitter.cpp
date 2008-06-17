
#include <fitter.h>

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
 
  reset();
  
}

//*************************************************************
bool fitter::initialize(const bhep::sstore& run_store) {
//*************************************************************
    
  m.message("+++ fitter init  function ++++",bhep::VERBOSE);
  
  // Initialise fail counters for fit stats.
  totFitAttempts = 0;
  fitSucceed = 0;
  toomany = 0;
  toofew = 0;
  kink = 0;
  nonFid = 0;
  patFail = 0;
  
  // read parameters
  readParam();
  
  // initialize geometry
  
  geom.init(store,run_store,level);
  
  // initialize geometrical limits
  man().geometry_svc().set_zero_length(1e-5 * mm);
  man().geometry_svc().set_infinite_length(1e12 * mm);
  
  // select the fitter
  
  man().fitting_svc().select_fitter(kfitter);
  
  // turn off MS if X0 = 0
  
  if (X0 == 0){
    man().model_svc().enable_noiser(model, RP::ms, false);
  }
  
  // set maximum local chi2
  
  man().fitting_svc().retrieve_fitter<KalmanFitter>(kfitter,model).
    set_max_local_chi2ndf(chi2node_max);
  
  man().fitting_svc().retrieve_fitter<KalmanFitter>(kfitter,model).
    set_number_allowed_outliers(max_outliers);
  
  // create the experimental setup
        
  create_setup();
  
  // don't propagate to surface with no measurement
  
  man().navigation_svc().navigator(model).set_unique_surface(true);
  
  // set verbosity of recpack services 
  
  setVerbosity(vfit,vnav,vmod);
  
  // 
  define_pattern_rec_param();

  m.message("+++ End of init function ++++",bhep::VERBOSE);
  
  return true;
}



//*************************************************************
void fitter::create_setup() {
//*************************************************************
    
   
  // add the setup to the geometry service
  man().geometry_svc().add_setup("main",geom.setup());
  
  // select the setup to be used by the geometry service
  man().geometry_svc().select_setup("main");
  
  
}

//*************************************************************
bool fitter::execute(bhep::particle& part,State seed,bool tklen){
//*************************************************************
  
  m.message("+++ fitter execute function ++++",bhep::VERBOSE);
  
  // create seed state
  
  userseed=true;
  
  setSeed(seed.vector());
  
  return execute(part,tklen);
  
}

//*************************************************************
bool fitter::execute(bhep::particle& part,int evNo, bool tklen){
//*************************************************************
    
  m.message("+++ fitter execute function ++++",bhep::VERBOSE);
  
  bool ok; 
  bool fitted=false;
  totFitAttempts++;
  _failType = 0; //set to 'success' before run to avoid faults in value.
  ok = readTrajectory(part);
  // if (evNo==1292 || evNo==2440 || evNo==4636 || evNo==5344){
//     std::string fuck = "VERBOSE";
//     Messenger::Level maw = Messenger::str(fuck);
//     man().fitting_svc().fitter(model).set_verbosity(maw);
  if (!userseed && ok) computeSeed();
  
  if (ok) {
    
    fitted = fitTrajectory(seedstate);
    
    addFitInfo(part,fitted);
    
    if (tklen) addTrackLength(part,_traj);
    
    if (fitted) m.message("++ Particle fitted",bhep::VERBOSE);
    
    else m.message("++ Particle not fitted !!!",bhep::VERBOSE);
    if (!fitted) {
      
      if (_traj.last_fitted_node() < (int)_meas.size()-1 && _failType!=3)
	{  kink++;
	_failType = 4; }
      
      m.message("++Failed fit trajectory++",bhep::NORMAL);
    }
  }
  else m.message("++ Particle lies outside no. hit cuts!",bhep::VERBOSE);
  // std::string shite = "MUTE";
//   Messenger::Level paw = Messenger::str(shite);
//   man().fitting_svc().fitter(model).set_verbosity(paw);
//   }
  userseed=false;
  
  if (fitted) {
    fitSucceed++; 
    if (_failType!=3) _failType = 0;
  } 
  
  return fitted;
  
}

//*************************************************************
void fitter::reset() {
//*************************************************************
  
  m.message("+++ reset function +++",bhep::VERBOSE);
  
  //reset trajectory 
  
  _traj.reset(); stc_tools::destroy(_meas);
  
  //reset virtual planes
  
  resetVirtualPlanes();
  
  
}


//*************************************************************
void fitter::resetVirtualPlanes() {
//*************************************************************
  
  m.message("+++ resetVirtualPlanes function +++",bhep::VERBOSE);
  
  for (size_t i=0; i<pnumber; i++){
    
    const dict::Key sname = "VPLANE_"+bhep::to_string(i);
    m.message("Deleting virtual surface: ",sname,bhep::VERBOSE);
    geom.setup().remove_surface(sname);
    delete virtual_planes[i];
  }
    
  virtual_planes.clear(); virtual_planes.resize(0); 
  
  m.message("++ Virtual surfaces deleted",bhep::VERBOSE);
  
  pnumber=0;
  
}

//***********************************************************************
void fitter::addFitInfo(bhep::particle& part,bool fitted) {
//***********************************************************************
 
  m.message("+++ addFitInfo function +++",bhep::VERBOSE);
  
  bool hasproperty = part.find_property("fitted");
  
  if (fitted) {
    
    m.message("++ Particle fitted",bhep::VERBOSE);
    m.message("++ Chi2 of fit: ",getChi2(),bhep::VERBOSE);
    
    if (hasproperty) {
      
      part.change_property("fitted","1");  
      part.change_property("fitChi2",bhep::to_string(getChi2())); 
      part.change_property("ndof",bhep::to_string(_traj.ndof())); 
      part.change_property("charge",bhep::to_string(getQ()));
      
    }
    
    else{

      part.add_property("fitted","1");  
      part.add_property("fitChi2",bhep::to_string(getChi2()));
      part.add_property("ndof",bhep::to_string(_traj.ndof())); 
      part.add_property("charge",bhep::to_string(getQ()));
      part.add_property("length","0");  
      
    }

  }
    
  else{ //particle not fitted

    if (hasproperty){
      part.change_property("fitted","0"); 
      part.change_property("fitChi2","-999"); 
      part.change_property("ndof","0"); 
      part.change_property("charge","0"); 
      part.change_property("length","0");  
    }
    else{
      part.add_property("fitted","0");
      part.add_property("fitChi2","-999"); 
      part.add_property("ndof","0"); 
      part.add_property("charge","0"); 
      part.add_property("length","0");  
    }    
  }
    

}

//*************************************************************
bool fitter::fitTrajectory(State seed) {
//*************************************************************
    
    m.message("+++ fitTrajectory function ++++",bhep::VERBOSE);
    
    bool ok = man().fitting_svc().fit(seed,_traj);

    if (ok && refit){
        
        ok = checkQuality(); if (!ok) return ok;

        m.message("Going to refit...",bhep::VERBOSE);
	
        //--------- refit using a new seed --------//	
	State newstate = _traj.state(_traj.first_fitted_node());
	man().model_svc().model(RP::particle_helix).representation()
	  .convert(newstate, RP::slopes_z);

	EVector v = newstate.vector();
	EMatrix C0 = newstate.matrix();

        EMatrix C = setSeedCov(C0,facRef);

	man().model_svc().model(RP::particle_helix).representation()
	  .convert(seedstate, RP::slopes_z);
	seedstate.set_hv(HyperVector(v,C)); 

	man().model_svc().model(RP::particle_helix).representation(RP::slopes_z)
	  .convert(seedstate,RP::default_rep);
	//seedstate.keepDiagonalMatrix();

        ok = man().fitting_svc().fit(seedstate,_traj);
	
    }
    
    if (ok) ok = checkQuality();

    return ok;

}



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
    ok = find_muon_pattern();
    if (!ok)
      patFail++;
  }

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

    int lowPass = store.fetch_istore("low_Pass_hits");
    //--------- take hits from particle -------//
    
    const vector<bhep::hit*>& hits = p.hits("MIND");     
    if ((int)hits.size() < lowPass ) { 
      toofew++; 
      _failType = 1;
      return false;
    }
    
    //------------- loop over hits ------------//
 
    for(size_t j=0; j< hits.size(); j++){

      bhep::hit& hit = *hits[j];
        	        
      //---------- create measurament ---------------//
      
      Measurement* mnt = getMeasurement(hit);
      
      //---------end of create measurement-----------//
      
      _meas.push_back(mnt); 
      
      m.message("Measurement added:",*mnt,bhep::VVERBOSE);
      
    }//end of loop over hits
    
    
    //--------- add measurements to trajectory --------//
   
    if (patternRec) {
      if ((int)hits.size() < min_seed_hits) {toofew++; _failType = 1; return false;}
      sort( _meas.begin(), _meas.end(), reverseSorter() );
      //sortingReverseByZ() defined in recpack/Trajectory.h>
    } else {
      _traj.add_measurements(_meas);

      m.message("Trajectory created:",_traj,bhep::VVERBOSE);
    }

    return true;
}

//*************************************************************
bool fitter::check_valid_traj() {
//*************************************************************

  double zCut = store.fetch_dstore("z_cut");
  double xCut = store.fetch_dstore("x_cut");
  double yCut = store.fetch_dstore("y_cut");

  //--------- Reject too many hits --------//
  int highPass = store.fetch_istore("high_Pass_hits");
  
  if ((int)_traj.nmeas() > highPass) { 
    toomany++;
    _failType = 2; 
    return false; }
  
  //---- Reject if initial meas outside fid. Vol ----//
  if (_traj.nodes()[0]->measurement().surface().position()[2] > geom.getPlaneZ()/2-zCut*cm
      || fabs(_traj.nodes()[0]->measurement().vector()[0]) > geom.getPlaneX()/2-xCut*cm
      || fabs(_traj.nodes()[0]->measurement().vector()[1]) > geom.getPlaneY()/2-yCut*cm)
    { nonFid++; 
    _failType = 3;}
    
  return true;
}

//*************************************************************
int fitter::getQ(){
//*************************************************************
  
  if (model.compare("particle/helix")!=0) return 0;
  
  double q = _traj.state(_traj.last_fitted_node()).vector()[dim-1];
  
  if (q<0) q=-1; else q=1;

  return (int) q;

}


//*************************************************************
Measurement*  fitter::getMeasurement(bhep::hit& hit){
//*************************************************************
    
    m.message("+++ getMeasurement function ++++",bhep::VERBOSE);
    
    //---- generate a virtual plane to hold the hit ----//
    
    bhep::Point3D bhit_pos = hit.x(); 
    EVector pos(3,0); pos[2] = bhit_pos[2];
    
    EVector zaxis = geom.getZaxis();
    EVector xaxis = geom.getXaxis();
    double height = geom.getPlaneX();
    double width  = geom.getPlaneY();
    string meastype = geom.getMeasType();
    EMatrix cov = geom.getCov();
    
    const dict::Key surf_name = "VPLANE_"+bhep::to_string(pnumber);
    Surface* surf = new Rectangle(pos,zaxis,xaxis,height/2,width/2);
    geom.setup().add_surface("Detector",surf_name,surf);
    geom.setup().set_surface_property(surf_name,"measurement_type",meastype);
    geom.setup().set_surface_property(surf_name,"resolution",cov);
    m.message("++ Adding virtual plane: ",surf_name,bhep::VERBOSE);
    pnumber++;
    
    virtual_planes.push_back(surf);

    //----- generate repack hit from bhep one ----//
    
    EVector hit_pos(2,0);
    hit_pos[0]=bhit_pos[0];
    hit_pos[1]=bhit_pos[1];
    
    Measurement* me = new Measurement();
    me->set_name(meastype);
    me->set_hv(HyperVector(hit_pos,cov));
    me->set_surface(geom.setup().surface(surf_name));
    me->set_position(geom.setup().surface(surf_name).position());
    //Add a key to the measurement with the true mother particle for PatRec.
    if (patternRec){
      const dict::Key motherP = "MotherParticle";
      const dict::Key mothName = hit.mother_particle().name();
      me->set_name(motherP, mothName);
    }
    
      
    return me;

}

//*************************************************************
void fitter::addTrackLength(bhep::particle& p,const Trajectory& t){
//*************************************************************
   
  m.message("+++ addTrackLength function +++",bhep::VERBOSE);
   
  if (p.find_property("length")){ 
    
    p.change_property("length",bhep::to_string(trackLength(t)));  
  }
  else p.add_property("length",bhep::to_string(trackLength(t)));  
  
  
}

//*************************************************************
double fitter::trackLength() {
//************************************************************
  
  return trackLength(_traj);
  
}

//*************************************************************
double fitter::trackLength(const Trajectory& t) {
//************************************************************
  
  m.message("+++ trackLength function +++",bhep::VERBOSE);
  
  double tot_length=0;
  double length=0;
  
  vector<Node*> nodes = t.nodes();
  
  if (!nodes.size()) return length;
  
  size_t current = t.first_fitted_node();
  size_t last = t.last_fitted_node();
  
  if (current==last) return length;

  size_t next = current + 1;
  State cState = State(nodes[current]->state());
  
  for (size_t in = next; in<=last; in++){ 
    if (!nodes[in]->status("fitted")) continue; 
    
    Measurement me = nodes[in]->measurement(); 
    
    State tempState = State(cState); 
    man().navigation_svc().propagate(me.surface(),tempState,length);
   
    cState = State(nodes[in]->state());
    
    tot_length += length;  
  }
    
  return abs(tot_length);

}


//*************************************************************
bool fitter::finalize() {
//*************************************************************
  
  ofstream fitstats;
  fitstats.open("MindFitStats.txt");

  fitstats << "Type of fail \t No. Events" << endl
	   << "Total Fits: \t" << totFitAttempts << endl
	   << "Successes: \t"  << fitSucceed     << endl
	   << "Too few hit: \t"<< toofew         << endl
	   << "Too many: \t"   << toomany        << endl
	   << "With Kink: \t"  << kink           << endl
	   << "Patrec. fail: \t"<< patFail       << endl
	   << "Outside Fid: \t"<< nonFid         << endl;
  fitstats.close();

  return true;
}


//*************************************************************
void fitter::computeSeed() {
//*************************************************************
    
    m.message("+++ computeSeed function ++++",bhep::VERBOSE);
    
    //use position slightly offset from first meas as seed 

    EVector v(3,0); 
    
    v[0] = _traj.nodes()[0]->measurement().vector()[0];
    v[1] = _traj.nodes()[0]->measurement().vector()[1];
    v[2] = _traj.nodes()[0]->measurement().surface().position()[2];   

    setSeed(v);
    
    m.message("++ Seed estate:",seedstate,bhep::VERBOSE);

}

//*************************************************************
void fitter::setSeed(EVector r, double factor){
//*************************************************************
  
  m.message("+++ setSeed function ++++",bhep::VERBOSE);
  
  /*
  EMatrix C(dim,dim,0);
  C = setSeedCov(r,factor);
    
  double p = 5 * GeV;
  double q = 1;
  
  int sens = 1; // ????? 
  
  EVector u(3,0); u[0]=0.0; u[1]=0.0; u[2]=sens; u/=u.norm();
  EVector eu(3,0); eu[0]=C[3][3]; eu[1]=C[4][4];
  EVector er(3,0); er[0]=C[0][0]; er[1]=C[1][1]; er[2]=C[2][2]; 
  double ep = 100 * GeV; // ???? total momentum error
  
  seedstate = ParticleState(r,u,p,q,er,eu,ep,model);

  seedstate.keepDiagonalMatrix();
  */

  // take as seed the state vector
  EVector v(6,0), v2(1,0);
  EMatrix C(6,6,0), C2(1,1,0);
  EVector dr(3,0); EVector dr2(3,0);

  v[0]=r[0];
  v[1]=r[1];
  v[2]=r[2];

  find_directSeed(dr, 1);
  find_directSeed(dr2,2);
  //Compare the two possible direction seeds.
  double Dx = (v[0]+(dr[0]/dr[2])*(_traj.nodes()[4]->measurement().surface().position()[2]-v[2]))
    - (v[0]+(dr2[0]/dr2[2])*(_traj.nodes()[4]->measurement().surface().position()[2]-v[2]));
  double Dy = (v[1]+(dr[1]/dr[2])*(_traj.nodes()[4]->measurement().surface().position()[2]-v[2]))
    - (v[1]+(dr2[1]/dr2[2])*(_traj.nodes()[4]->measurement().surface().position()[2]-v[2]));

  if (fabs(Dx)< 3*cm && fabs(Dy)< 3*cm){
    dr += dr2;
    dr /= dr.norm();
    v[3] = dr[0]/dr[2];// + dr2[0]/dr2[2]);
    v[4] = dr[1]/dr[2];// + dr2[1]/dr2[2]);
  } else {
    v[3] = dr[0]/dr[2];
    v[4] = dr[1]/dr[2];
  }

  //Approximate p from plot of p vs. no. hits, then approx. de_dx from this.
  double pSeed = (double)(0.060*_traj.nmeas())*GeV;
  double de_dx = -7.87*(0.013*(pSeed/GeV)+1.5)*MeV/cm;
  geom.setDeDx(de_dx);
  
  m.message("reset energy loss to approx value",bhep::VERBOSE);

  //v[5]=dr[2];//1;
  v[5]=-1/pSeed;

  // But use a larger covariance matrix
  // diagonal covariance matrix
  C[0][0] = C[1][1] = 9.*cm*cm;
  C[2][2] = EGeo::zero_cov()/2;
  C[3][3] = C[4][4] = 1.;
  //C[5][5] = 10.;
  C[5][5] = pow(v[5],2)*3;

  seedstate.set_name(RP::particle_helix);
  seedstate.set_name(RP::representation,RP::slopes_z);
  //seedstate.set_name(RP::representation,RP::default_rep);
  v2[0] = 1;
  seedstate.set_hv(RP::sense,HyperVector(v2,C2));
  seedstate.set_hv(HyperVector(v,C));
  
  man().model_svc().model(RP::particle_helix).representation(RP::slopes_z)
    .convert(seedstate,RP::default_rep);

}

//*************************************************************
EMatrix fitter::setSeedCov(EVector v, double factor){
//*************************************************************
    
    //--- a large diagonal covariance matrix ---//
    
    EMatrix c(dim,dim,0);

    double p=5*GeV; // ????
  
    if (model.compare("particle/helix")==0){ 

      c[0][0] = c[1][1] = pow(100*100*mm*mm,2)/factor; // ???? 
      
      c[2][2] = EGeo::zero_cov()/2;//no error in z    
   
      c[3][3] = c[4][4] = 3/factor;     
      //      c[6][6] = pow(1/p,2)/factor;   
      c[6][6] = pow(1/p,2)*100/factor;   
      
    }
    

    return c;
}


//*************************************************************
EMatrix fitter::setSeedCov(EMatrix C0, double factor){
//*************************************************************
    
    //--- a large diagonal covariance matrix ---//
    
  EMatrix c = factor*C0;


    return c;
}

//*************************************************************
void fitter::find_directSeed(EVector& R, int sense){
//*************************************************************

  int farPos;

  if (sense==1){
    if (_traj.nmeas()>10) farPos = 9;
    else farPos = 4;
    R[0] = (_traj.nodes()[farPos]->measurement().vector()[0]
	    - _traj.nodes()[0]->measurement().vector()[0]);
    R[1] = (_traj.nodes()[farPos]->measurement().vector()[1]
	    - _traj.nodes()[0]->measurement().vector()[1]);
    R[2] = (_traj.nodes()[farPos]->measurement().surface().position()[2]
	    - _traj.nodes()[0]->measurement().surface().position()[2]);
  }
  if (sense==-1){
    R[0] = _meas[0]->vector()[0] - _meas[1]->vector()[0];
    R[1] = _meas[0]->vector()[1] - _meas[1]->vector()[1];
    R[2] = _meas[0]->surface().position()[2] 
      - _meas[1]->surface().position()[2];
  }
  if (sense==2){
    R[0] = (_traj.nodes()[2]->measurement().vector()[0]
	    - _traj.nodes()[0]->measurement().vector()[0]);
    R[1] = (_traj.nodes()[2]->measurement().vector()[1]
	    - _traj.nodes()[0]->measurement().vector()[1]);
    R[2] = (_traj.nodes()[2]->measurement().surface().position()[2]
	    - _traj.nodes()[0]->measurement().surface().position()[2]);
  }

  R /= R.norm();
  
}

//*****************************************************************************
void fitter::readParam(){
//*****************************************************************************
    
    m.message("+++ readParam function of fitter ++++",bhep::VERBOSE);
        
    model="particle/helix"; 
    dim=7; // ??????
    

    if (store.find_sstore("fitter"))
      kfitter = store.fetch_sstore("kfitter");
    else kfitter="kalman";

    if (store.find_istore("refit"))
      refit=store.fetch_istore("refit");
    else refit=false;

    if (store.find_istore("patRec"))
      patternRec=store.fetch_istore("patRec");
    else patternRec=false;

    facRef = store.fetch_dstore("facRef");

    min_seed_hits = store.fetch_istore("min_seed_hits");
    max_seed_hits = store.fetch_istore("max_seed_hits");

    chi2node_max = store.fetch_dstore("chi2node_max");
    
    max_outliers = store.fetch_istore("max_outliers");

    chi2fit_max = store.fetch_dstore("chi2fit_max");
    
    patRec_maxChi = store.fetch_dstore("pat_rec_max_chi");
    patRec_max_outliers = store.fetch_istore("pat_rec_max_outliers");
    max_consec_missed_planes = store.fetch_istore("max_consec_missed_planes");

    X0 = store.fetch_dstore("x0") * mm;
    
    vfit = store.fetch_istore("vfit");
    vnav = store.fetch_istore("vnav");
    vmod = store.fetch_istore("vmod");
    
}

//*****************************************************************************
void fitter::setVerbosity(int v0,int v1,int v2){
//*****************************************************************************
  
  m.message("+++ setVerbosity function ++++",bhep::VERBOSE);

  std::string info[4]={"MUTE","NORMAL","VERBOSE","VVERBOSE"};

  Messenger::Level l0 = Messenger::str(info[v0]);
  Messenger::Level l1 = Messenger::str(info[v1]);
  Messenger::Level l2 = Messenger::str(info[v2]);
  

  // verbosity levels related with fitting
  man().fitting_svc().fitter(model).set_verbosity(l0);

  // verbosity levels related with navigation
  man().navigation_svc().set_verbosity(l1);
  man().navigation_svc().navigator(model).set_verbosity(l1);
  man().navigation_svc().inspector("ms").set_verbosity(l1);
  man().navigation_svc().navigator(model).master_inspector().set_verbosity(l1);
  man().navigation_svc().inspector("BField").set_verbosity(l1);  
  man().navigation_svc().inspector("eloss").set_verbosity(l1);     
  man().model_svc().model(model).intersector("plane").set_verbosity(l1);

  // verbosity levels related with model operation (soft intersection)
  man().model_svc().model(model).equation().set_verbosity(l2);
  man().model_svc().model(model).propagator().set_verbosity(l2);
  man().model_svc().model(model).tool("noiser/ms").set_verbosity(l2);


}

//*****************************************************************************
void fitter::define_pattern_rec_param() {
//*****************************************************************************

  patman().geometry_svc().set_zero_length(1e-5 * mm);
  patman().geometry_svc().set_infinite_length(1e12 * mm);

  // add the setup to the geometry service
  patman().geometry_svc().add_setup("main",geom.setup());
  
  // select the setup to be used by the geometry service
  patman().geometry_svc().select_setup("main");
  
  patman().navigation_svc().navigator(model).set_unique_surface(true);
  
  patman().fitting_svc().retrieve_fitter<KalmanFitter>(kfitter,model).
    set_max_local_chi2ndf(patRec_maxChi);

  patman().fitting_svc().retrieve_fitter<KalmanFitter>(kfitter,model).
    set_number_allowed_outliers(patRec_max_outliers);

}

//*****************************************************************************
bool fitter::find_muon_pattern() {
//*****************************************************************************

//Algorithms for isolating possible muon from event hits.
  m.message("++ Starting Pattern Recognition ++",bhep::VERBOSE);

  _patRecStat.clear();
  int measCount = (int)_meas.size();
  _patRecStat = vector<bool> (measCount, 0);
  _recChi = EVector(3,0);
  _recChi[1] = 100000; //SetLarge dummy value for minChiHadron bit.

  State patternSeed;
  
  bool ok = get_patternRec_seedtraj();
  if (ok)
    ok = get_patternRec_seed(patternSeed);

  if (ok)
    ok = perform_pattern_rec(patternSeed);
  
  // Sort the measurements identified as muon hits into ascending
  // z for muon fit.
  _traj.sort_nodes(1);

  return ok;
}

//*****************************************************************************
bool fitter::get_patternRec_seedtraj() {
//*****************************************************************************

  double tolerance = 1 * cm;
  int nloops;
  double currentZ, prevZ;

  int nMeas = (int)_meas.size();
  if (nMeas < max_seed_hits) nloops = nMeas;
  else nloops = max_seed_hits;

  prevZ = _meas[0]->surface().position()[2];
  
  for (int Iso = 1;Iso < nloops;Iso++) {

    currentZ = _meas[Iso]->surface().position()[2];
    
    if (currentZ >= (prevZ-tolerance)) break;
    else {
      
      _traj.add_measurement(*_meas[Iso-1]);
      if (_meas[Iso-1]->name("MotherParticle").compare("Hadronic_vector")!=0)
	_patRecStat[Iso-1] = true;
      prevZ = currentZ;
      
    }

  }
  
  if ((int)_traj.nmeas() < min_seed_hits) { 
    cout<<"No seedTraj found"<<endl;
    _failType = 5;
    return false;}
  
  iGroup = (int)_traj.nmeas();

  return true;

}

//*****************************************************************************
bool fitter::get_patternRec_seed(State& seed) {
//*****************************************************************************
  
  int lastMeas = (int)_meas.size() - 1;
  
  EVector V(6,0); EVector V2(1,0);
  EMatrix M(6,6,0); EMatrix M2(1,1,0);
  EVector dr(3,0);

  V[0] = _meas[0]->vector()[0];
  V[1] = _meas[0]->vector()[1];
  V[2] = _meas[0]->surface().position()[2];
  
  find_directSeed(dr, -1);
  V[3] = dr[0]/dr[2];
  V[4] = dr[1]/dr[2];
  
  //Approximate momentum by extent (hard wired for 5cm iron/scint
  //separation at mo) using same empirical functions as in isolated
  //muon fitting.
  double Xtent = (_meas[0]->surface().position()[2]
    - _meas[lastMeas]->surface().position()[2])/5;
  double pSeed = (0.060*Xtent)*GeV;
  double de_dx = -7.87*(0.013*(pSeed/GeV)+1.5)*MeV/cm;
  geom.setDeDx(de_dx);
  
  V[5] = 1./pSeed;
  
  //Errors
  M[0][0] = M[1][1] = 15.*cm*cm;
  M[2][2] = EGeo::zero_cov()/2;//1.*cm*cm;
  M[3][3] = M[4][4] = 1.5;
  M[5][5] = pow(V[5],2)*4;

  //Seedstate fit properties
  seed.set_name(RP::particle_helix);
  seed.set_name(RP::representation,RP::slopes_z);
  V2[0] = 1; 
  seed.set_hv(RP::sense,HyperVector(V2,M2));
  seed.set_hv(HyperVector(V,M));
  
  patman().model_svc().model(RP::particle_helix).representation(RP::slopes_z).convert(seed,RP::default_rep);
  
  //bool ok = perform_least_squares(seed);
  bool ok = perform_kalman_fit(seed);
  if (!ok)
    _failType = 5;

  return ok;
}

//*****************************************************************************
bool fitter::perform_least_squares(State& seed) {
//*****************************************************************************

  m.message("++ About to perform least squares fit ++",bhep::VERBOSE);
  
  std::string YerMaw = "NORMAL";
  Messenger::Level blah = Messenger::str(YerMaw);
  patman().fitting_svc().fitter(model).set_verbosity(blah);
  patman().fitting_svc().select_fitter("lsq");
  
  bool ok = patman().fitting_svc().fit(seed, _traj);
  
  //seed.clear();
  if (ok)
    seed = _traj.state(_traj.first_fitted_node());
  
  return ok;
}

//*****************************************************************************
bool fitter::perform_kalman_fit(State& seed) {
//*****************************************************************************

  m.message("++ About to perform kalman fit ++",bhep::VERBOSE);
  
  std::string YerMaw = "NORMAL";
  Messenger::Level blah = Messenger::str(YerMaw);
  patman().fitting_svc().fitter(model).set_verbosity(blah);
  patman().fitting_svc().select_fitter("kalman");
  
  bool ok = patman().fitting_svc().fit(seed, _traj);
  
  //seed.clear();
  if (ok)
    seed = _traj.state(_traj.first_fitted_node());
  // cout << "Seed Check: "<<endl
//        <<t<<endl
//        <<seed<<endl;
  return ok;
}

//*****************************************************************************
bool fitter::perform_pattern_rec(const State& seed) {
//*****************************************************************************

  m.message("++ About to perform filter on hits ++",bhep::VERBOSE);
  
  measurement_vector NeedFiltered;

  bool ok;
  double curZ, preZ;
  double tolerance = 1 * cm;
  _nConsecHoles = 0;

  //fitter parameters.
  patman().fitting_svc().select_fitter("kalman");

  while ( iGroup < (int)_meas.size() ) {

    if ((int)NeedFiltered.size() == 0){
      NeedFiltered.push_back( _meas[iGroup] );
      iGroup++;
      continue;
    }

    curZ = _meas[iGroup]->surface().position()[2];
    preZ = _meas[iGroup-1]->surface().position()[2];

    if (curZ < (preZ-tolerance)) {
      if (NeedFiltered.size() > 1) {
	
	ok = filter_close_measurements(NeedFiltered, seed);
	if (!ok) {
	  cout<<"Failed in filtering"<<endl; 
	   _failType = 6;
	  return ok;}

      } else {
	
	ok = patman().fitting_svc().filter(*NeedFiltered[0], seed, _traj);

	if (ok){
	  _patRecStat[iGroup-1] = true;
	  
	  if (_meas[iGroup-1]->name("MotherParticle").compare("Hadronic_vector")!=0){
	    
	    _recChi[0] = TMath::Max(_traj.node(_traj.last_fitted_node()).quality(), _recChi[0]);
	  }
	  else _recChi[1] = TMath::Min(_traj.node(_traj.last_fitted_node()).quality(), _recChi[1]);

	  _recChi[2] = TMath::Max(_nConsecHoles, (int)_recChi[2]);
	  _nConsecHoles = 0;

	} else {
	  _nConsecHoles++;

	  if (_nConsecHoles > max_consec_missed_planes) {
	    _failType = 6;
	    _recChi[2] = _nConsecHoles;
	    return ok;
	  }
	}

      }
	
      NeedFiltered.clear();
      NeedFiltered.push_back( _meas[iGroup] );
      iGroup++;

    } else {
      NeedFiltered.push_back( _meas[iGroup] );
      iGroup++;
    }

  }

  return ok;

}

//*****************************************************************************
bool fitter::filter_close_measurements(measurement_vector& Fmeas,
				       const State& seed) {
//*****************************************************************************
  
  bool ok;
  const int nMeas = (int)Fmeas.size();
  double Chi2[nMeas];

  // Filter measurements in Fmeas. Delete all but the lowest chi2
  // from traj and cast the corresponding measurement to the
  // hadron measurement vector.

  
  for (int iMat = 0;iMat < nMeas;iMat++) {

    ok = patman().matching_svc().match_trajectory_measurement(_traj, *Fmeas[iMat],
							      Chi2[iMat]);
  }

  long ChiMin = TMath::LocMin(nMeas, Chi2);

  for (int iFilt = 0;iFilt < nMeas;iFilt++) {

    if (iFilt == (int)ChiMin){

      ok = patman().fitting_svc().filter(*Fmeas[(int)iFilt], seed , _traj);

      if (ok && _meas[iGroup-(nMeas-(int)ChiMin)]->name("MotherParticle")
	  .compare("Hadronic_vector")!=0) {

	_patRecStat[iGroup-(nMeas-(int)ChiMin)] = true;

	_recChi[0] = TMath::Max(Chi2[iFilt], _recChi[0]);

      }
      else if (ok && _meas[iGroup-(nMeas-(int)ChiMin)]->name("MotherParticle")
	       .compare("Hadronic_vector")==0)
	_recChi[1] = TMath::Min(Chi2[iFilt], _recChi[1]);

      else cout<< "Filter failed" <<endl;
    }
    else
      _hadmeas.push_back( Fmeas[(int)iFilt] );

  }

  if (ok) {
    _recChi[2] = TMath::Max(_nConsecHoles, (int)_recChi[2]);
    _nConsecHoles = 0;
  } else {
    _nConsecHoles++;

    if ( _nConsecHoles > max_consec_missed_planes) {
      _recChi[2] = _nConsecHoles;
      _failType = 6;
      return false;
    }
  }
  
  return true;
}
