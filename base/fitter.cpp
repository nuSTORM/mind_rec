
#include <fitter.h>


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

    // read parameters
    totFitAttempts = 0;
    fitSucceed = 0;
    toomany = 0;
    toofew = 0;
    kink = 0;
    unkFail = 0;

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

    //    man().fitting_svc().retrieve_fitter<KalmanFitter>(kfitter,model).
    //      set_number_allowed_outliers(500);

     // create the experimental setup
        
    create_setup();
    
    // don't propagate to surface with no measurement

    man().navigation_svc().navigator(model).set_unique_surface(true);

    // set verbosity of recpack services 

    setVerbosity(vfit,vnav,vmod);
    
    // create seed state
    
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
bool fitter::execute(bhep::particle& part,bool tklen){
//*************************************************************
    
    m.message("+++ fitter execute function ++++",bhep::VERBOSE);

    bool ok; 
    bool fitted=false;
    totFitAttempts++;
    ok = readTrajectory(part);

    if (!userseed && ok) computeSeed();

    if (ok) {
      
      fitted = fitTrajectory(seedstate);

      addFitInfo(part,fitted);
      
      if (tklen) addTrackLength(part,traj);
      
      if (fitted) m.message("++ Particle fitted",bhep::VERBOSE);
      
      else m.message("++ Particle not fitted !!!",bhep::VERBOSE);
      if (!fitted) {
	int nFittedNodes = 0;
	for (int inode=0;inode<(int)meas.size();inode++){
	  if (traj.node(inode).status("fitted"))
	    nFittedNodes++;
	}
	if (nFittedNodes < (int)meas.size())
	  {  kink++; cout << "Traj. failed after kink found" << endl;
	  cout << traj << endl; }
	else { unkFail++; cout << "Unknown Fail" << endl;
	cout << traj << endl; }
	m.message("++Failed fit trajectory++",bhep::NORMAL);
      }
    }
    else m.message("++ Particle lies outside no. hit cuts!",bhep::VERBOSE);

    userseed=false;

    if (fitted) fitSucceed++;  

    return fitted;

}

//*************************************************************
void fitter::reset() {
//*************************************************************
  
  m.message("+++ reset function +++",bhep::VERBOSE);

  //reset trajectory 
   
  traj.reset(); stc_tools::destroy(meas);
  
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
      part.change_property("ndof",bhep::to_string(traj.ndof())); 
      part.change_property("charge",bhep::to_string(getQ()));
      
    }
    
    else{

      part.add_property("fitted","1");  
      part.add_property("fitChi2",bhep::to_string(getChi2()));
      part.add_property("ndof",bhep::to_string(traj.ndof())); 
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
    
    bool ok = man().fitting_svc().fit(seed,traj);

    if (ok && refit){
        
        ok = checkQuality(); if (!ok) return ok;

        m.message("Going to refit...",bhep::VERBOSE);
	
        //--------- refit using a new seed --------//	
	State newstate = traj.state(traj.first_fitted_node());
	man().model_svc().model(RP::particle_helix).representation()
	  .convert(newstate, RP::slopes_z);

	EVector v = newstate.vector();
	EMatrix C0 = newstate.matrix();

        EMatrix C = setSeedCov(C0,50.);

	man().model_svc().model(RP::particle_helix).representation()
	  .convert(seedstate, RP::slopes_z);
	seedstate.set_hv(HyperVector(v,C)); 

	man().model_svc().model(RP::particle_helix).representation(RP::slopes_z)
	  .convert(seedstate,RP::default_rep);
	//seedstate.keepDiagonalMatrix();

        ok = man().fitting_svc().fit(seedstate,traj);
	
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

    bool ok = recTrajectory(part,traj);
  
    return ok;
}


//*************************************************************
bool fitter::recTrajectory(const bhep::particle& p,Trajectory& t) {
//*************************************************************

    m.message("+++ recTrajectory function ++++",bhep::VERBOSE);
 
    reset();

    int lowPass = store.fetch_istore("low_Pass_hits");
    int highPass = store.fetch_istore("high_Pass_hits");

    //--------- take hits from particle -------//
    
    const vector<bhep::hit*>& hits = p.hits("MIND"); 
        
    if (hits.size()==0) { toofew++; return false; }
    if ((int)hits.size()>highPass) { toomany++; return false; }
    //------------- loop over hits ------------//
    
    vector<string> vplanes;
 
    for(size_t j=0; j< hits.size(); j++){

      //if (j>200) break;
      
      bhep::hit& hit = *hits[j];
        	        
      //---------- create measurament ---------------//
      
      Measurement* mnt = getMeasurement(hit);
      
      //---------end of create measurement-----------//
      
      meas.push_back(mnt); 
      
      m.message("Measurement added:",*mnt,bhep::VVERBOSE);
      
    }//end of loop over hits
    
    
    //--------- add measurements to trajectory --------//
   
    t.add_measurements(meas);
    
    m.message("Trajectory created:",t,bhep::VVERBOSE);
    
    if ((int)meas.size()<lowPass) { toofew++; return false; } //not enough dof
    
    ////sort traj measurements
    //double z1=geom.setup().surface(vplanes[0]).position()[2];
    //double z2=geom.setup().surface(vplanes[1]).position()[2];
    //if ( fabs(z1) >fabs(z2) ) t.sort_nodes((int) (z1/fabs(z1)));
    ////

    return true;

}


//*************************************************************
int fitter::getQ(){
//*************************************************************
  
  if (model.compare("particle/helix")!=0) return 0;
  
  double q = traj.state(traj.last_fitted_node()).vector()[dim-1];
  
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
  
  return trackLength(traj);
  
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
	   << "Too few hit: \t"<< toofew           << endl
	   << "Too many: \t"   << toomany          << endl
	   << "With Kink: \t"  << kink           << endl
	   << "Unknown: \t"    << unkFail        << endl;
  fitstats.close();

  return true;
}


//*************************************************************
void fitter::computeSeed() {
//*************************************************************
    
    m.message("+++ computeSeed function ++++",bhep::VERBOSE);
    
    //use position slightly offset from first meas as seed 

    EVector v(3,0); 
        
    v[0] = meas[0]->vector()[0];
    v[1] = meas[0]->vector()[1];
    v[2] = meas[0]->surface().position()[2];   

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
  EVector dr(3,0);

  v[0]=r[0];
  v[1]=r[1];
  v[2]=r[2];

  find_directSeed(dr);
  v[3] = dr[0]/dr[2];
  v[4] = dr[1]/dr[2];

  //Approximate p from plot of p vs. no. hits, then approx. de_dx from this.
  double pSeed = (double)(0.060*meas.size())*GeV;
  double de_dx = -7.87*(0.013*(pSeed/GeV)+1.5)*MeV/cm;
  geom.setDeDx(de_dx);
  //cout<<"pSeed = "<<pSeed<<", de_dx = "<<de_dx<<endl;
  //Empirical formulae: de_dex=density*pprox relation).
  //  geom.setup().set_volume_property("Detector","de_dx",de_dx);
  m.message("reset energy loss to approx value",bhep::NORMAL);

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
  
  man().model_svc().model(RP::particle_helix).representation(RP::slopes_z).convert(seedstate,RP::default_rep);

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
void fitter::find_directSeed(EVector& R){
//*************************************************************

  R[0] = meas[1]->vector()[0] - meas[0]->vector()[0];
  R[1] = meas[1]->vector()[1] - meas[0]->vector()[1];
  R[2] = meas[1]->surface().position()[2] - meas[0]->surface().position()[2];

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

    chi2node_max = store.fetch_dstore("chi2node_max");
    
    chi2fit_max = store.fetch_dstore("chi2fit_max");
    
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

