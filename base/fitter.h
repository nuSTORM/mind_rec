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
#include <recpack/LsqFitter.h>
#include <recpack/ParticleState.h>

#include <TH1F.h>
#include <TGraphErrors.h>
#include <TF1.h>

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
  bool execute(bhep::particle& part,int evNo,bool tklen=false);
  bool execute(bhep::particle& part,State seed,bool tklen=false);
  //end
  bool finalize() ;
  //-------------------------------------------------// 
    
  const Trajectory& get_traj(){
    if (reseed_ok == 0) return _traj;
    else return _traj2; }
  EVector& get_had_unit(){return _hadunit;}
  vector<double>& get_fit_tracker(){ return _fitTracker; }
  int get_fail_type(){return _failType;}
  EVector& get_PatRec_Chis(){return _recChi;};

  Measurement* get_meas(int num){return _meas[num];}
  int get_nMeas(){return (int)_meas.size();}

  double getChi2(){
    if (reseed_ok == 0) return _traj.quality();
    else return _traj2.quality(); }
  
  Measurement* getMeasurement(bhep::hit& hit);
  
  int getQ();
 
  //recpack manager
  
  RecpackManager& man(){return _man;}
  RecpackManager& patman(){return _patRecman;}
  
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
  void computeSeed(int firsthit=0);
  void setSeed(EVector v, int firsthit=0);
  void find_directSeed(EVector& R, int sense);
  void mom_from_parabola(int nplanes, int firsthit, EVector& V);

  //seed error
  EMatrix setSeedCov(EVector,double factor=1.);
  EMatrix setSeedCov(EMatrix C0, double factor);
 
  //fit trajectory
  bool fitTrajectory(State seed);
  bool reseed_traj();
  void check_starting_measurements();
  void remove_suspected_hads();
  bool fitHadrons();

  //Pattern recognition functions.
  void define_pattern_rec_param();
  bool find_muon_pattern();
  bool get_patternRec_seedtraj();
  bool get_patternRec_seed(State& seed);
  bool perform_least_squares(State& seed);
  bool perform_kalman_fit(State& seed);
  bool perform_pattern_rec(const State& seed);
  bool filter_close_measurements(measurement_vector& Fmeas,
				 const State& seed);

  //-------- get traj from event----------//
  bool readTrajectory(const bhep::particle& part);
  //get unfitted rec traj, i.e, measurements
  bool recTrajectory(const bhep::particle& part); 
  // Check traj passes cuts for fitting.
  bool check_valid_traj();
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
  RecpackManager _patRecman;

  vector<Surface*> virtual_planes;
  
  //counter for virtual planes
  size_t pnumber;
  
  bool userseed;
  
  //Parameters to define fitting method.
  bool refit; //Do second fit.
  bool patternRec; //Pattern recognition algorithm required?

  int min_seed_hits; //Minimum isolated hits required for Prec seed.
  int max_seed_hits; //Max. to be used.

  //Counters for fit fails and successes for various reasons.
  int totFitAttempts;
  int fitSucceed;
  int toomany;
  int toofew;
  int kink;
  int nonFid;
  int patFail;

  //counters to aid pattern rec.
  int iGroup;
  int _nConsecHoles;

  //bit to tell if reseed perfromed
  int reseed_ok;

  //------------------ Physics -----------------//
    
  double X0;

  int dim; //dimension of model state
  int meas_dim; //dimension of measurement
  
  State seedstate;   
  EVector qoverp;

  string model;  // fit model
  string kfitter; // kind of fit
    
  double chi2node_max;
  int max_outliers;
  double chi2fit_max;
  double facRef;

  double patRec_maxChi;
  int patRec_max_outliers;
  int max_consec_missed_planes;
    
  Trajectory _traj;
  Trajectory _traj2;
  measurement_vector _meas;
  measurement_vector _hadmeas;

  //Vectors to contain the relevant results of a pattern recognition run.
  //vector<bool> _patRecStat;
  EVector _recChi;

  //value set to identify where a traj failed:
  //0=too few hits. 1=too many hits. 2=outside fiducial. 3=no convergence with kink.
  //4=Couldn't find seed for patrec. 5=Failed in pat rec filtering
  int _failType;
  vector<double> _fitTracker;

  size_t nnodes;

  //Temporary fix (??). vector to store hadron unit dir vec.
  EVector _hadunit;
    
  //-------------- verbosity levels ------------//

  int vfit,vnav,vmod,vmat,vsim;
  
};

class reverseSorter{
public:
  bool operator()(const Measurement* p1, const Measurement* p2){
    if (p2->position()[2] < p1->position()[2]) return true;
    return false;
  }

};
  
#endif
