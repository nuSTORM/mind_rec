/* -*- mode: c++ -*- */
#ifndef _fitter___
#define _fitter___

#include <recpack/RecpackManager.h>
#include <mind/MINDsetup.h>
#include <mind/Utilities.h>
#include <bhep/event.h>
#include <bhep/gstore.h>
#include <recpack/RayTool.h>
#include <recpack/KalmanFitter.h>
#include <recpack/ParticleState.h>

using namespace Recpack;

class fitter{

public:
  
  fitter(const bhep::gstore& pstore,
	    bhep::prlevel vlevel=bhep::NORMAL);
    
  virtual ~fitter();
    
  //------------------ main functions ---------------// 
  //init
  bool initialize(const bhep::sstore&) ;
  //exe
  bool execute(bhep::particle& part,bool tklen=false);
  bool execute(bhep::particle& part,State seed,bool tklen=false);
  //end
  bool finalize() ;
  //-------------------------------------------------// 
    
  const Trajectory& get_traj(){return traj;}
  
  double getChi2(){return traj.quality();}
  
  Measurement* getMeasurement(bhep::hit& hit);
  
  int getQ();
 
  //recpack manager
  
  RecpackManager& man(){return _man;}
  
  //get tracklength
  
  double trackLength();
  double trackLength(const Trajectory&);
  void addTrackLength(bhep::particle&,const Trajectory&);
  
  //fit twice

  void setRefit(bool ok){refit=ok;}
   
protected:
  
  void resetVirtualPlanes(); 

  //generate Recpack Setup
  void create_setup();

  //recpack verbosity levels
  void setVerbosity(int v0,int v1,int v2);

  //read parameters from store
  void readParam();
    
  //seed for fit
  void computeSeed();
  void setSeed(EVector v, double factor=1.);
 
  //seed error
  EMatrix setSeedCov(EVector,double factor=1.);
 
  //fit trajectory
  bool fitTrajectory(State seed);
    
  //-------- get traj from event----------//
  bool readTrajectory(const bhep::particle& part);
  //get unfitted rec traj, i.e, measurements
  bool recTrajectory(const bhep::particle& part,Trajectory& t); 
  string getPlaneName(bhep::hit);
  //--------------------------------------//
  
  bool checkQuality();
   
  void addFitInfo(bhep::particle&,bool);
    
  void reset();

protected:

  bhep::gstore store;
  
  bhep::prlevel level;
    
  bhep::messenger m;
    
  MINDsetup geom;
  
  RecpackManager _man;
    
  vector<Surface*> virtual_planes;
  
  //counter for virtual planes
  size_t pnumber;
  
  bool userseed;
  
  bool refit;

  //------------------ Physics -----------------//
    
  double X0;

  int dim; //dimension of model state
  int meas_dim; //dimension of measurement
  
  State seedstate;   
  EVector qoverp;

  string model;  // fit model
  string kfitter; // kind of fit
    
  double chi2node_max;
  double chi2fit_max;
    
  Trajectory traj;
  measurement_vector meas;
    
  size_t nnodes;

    
  //-------------- verbosity levels ------------//

  int vfit,vnav,vmod,vmat,vsim;
  
};


#endif
