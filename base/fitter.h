/* -*- mode: c++ -*- */
#ifndef _fitter___
#define _fitter___

#include <recpack/RecpackManager.h>
#include <mind/MINDsetup.h>
#include <mind/Utilities.h>
#include <mind/MINDfitman.h>
#include <bhep/event.h>
#include <bhep/gstore.h>
#include <recpack/RayTool.h>
#include <recpack/KalmanFitter.h>
#include <recpack/LsqFitter.h>
#include <recpack/ParticleState.h>

#include <TH1F.h>
#include <TGraphErrors.h>
#include <TF1.h>

#include <mind/event_classif.h>

using namespace Recpack;

class fitter{
  
public:
  
  fitter(const bhep::gstore& pstore,
	    bhep::prlevel vlevel=bhep::NORMAL);
    
  virtual ~fitter();
    
  //------------------ main functions ---------------// 
  //init
  //bool initialize(const bhep::sstore&) ;
  bool initialize();
  //exe
  bool execute(bhep::particle& part,int evNo);
  //end
  bool finalize() ;
  //-------------------------------------------------// 
    
  const Trajectory& get_traj(){ 
    if (reseed_ok) return _traj2;
    else return _traj; }
  EVector& get_had_unit(){ return _hadunit; }
  double get_had_eng(){ return _hadEng; }
  int get_fail_type(){ return _failType; }
  bool check_reseed(){ return reseed_ok; }

  Measurement* get_meas(int num){return _meas[num];}
  int get_nMeas(){return (int)_meas.size();}

  double getChi2(){ return _traj.quality(); }
  
  Measurement* getMeasurement(bhep::hit& hit);
  
  int getQ();

  //calculate momentum from range.
  //void calculate_len_mom(double len, double *mom);
 
  //recpack manager
  
  RecpackManager& man(){
    return MINDfitman::instance().manager();}

  event_classif& get_classifier(){ return _classify; }
  
  //fit twice

  void setRefit(bool ok){refit=ok;}
   
protected:
  
  void resetVirtualPlanes(); 

  //read parameters from store
  void readParam();
    
  //seed for fit
  void computeSeed(int firsthit=0);
  void setSeed(EVector v, int firsthit=0);
  void mom_from_parabola(int nplanes, int firsthit, EVector& V);
  void set_de_dx(double mom);

  //seed error
  EMatrix setSeedCov(EMatrix C0, double factor);
 
  //fit trajectory
  bool fitTrajectory(State seed);
  bool reseed_traj();
  bool fitHadrons();
  double eng_scale(double visEng);

  //-------- get traj from event----------//
  bool readTrajectory(const bhep::particle& part);
  //get unfitted rec traj, i.e, measurements
  bool recTrajectory(const bhep::particle& part); 
  // Check traj passes cuts for fitting.
  bool check_valid_traj();
  string getPlaneName(bhep::hit);
  //--------------------------------------//
  
  bool checkQuality();
    
  void reset();

protected:

  bhep::gstore store;
  
  bhep::prlevel level;
    
  bhep::messenger m;
    
  MINDsetup geom;
  
  //counter for virtual planes
  size_t pnumber;
  
  //Parameters to define fitting method.
  bool refit; //Do second fit.
  bool patternRec; //Pattern recognition algorithm required?

  int min_seed_hits; //Minimum isolated hits required for Prec seed.
  double min_iso_prop;

  //bit to tell if reseed perfromed
  bool reseed_ok;

  //------------------ Physics -----------------//
    
  double X0;

  int dim; //dimension of model state
  int meas_dim; //dimension of measurement
  double _tolerance; //pos. resolution/plane tolerance
  
  State seedstate;   
  EVector qoverp;

  string model;  // fit model
  //string kfitter; // kind of fit
    
  //double chi2node_max;
  //int max_outliers;
  double chi2fit_max;
  double facRef;
    
  Trajectory _traj;
  Trajectory _traj2;
  measurement_vector _meas;
  measurement_vector _hadmeas;

  //value set to identify where a traj failed:
  //0=too few hits. 1=too many hits. 2=outside fiducial. 3=no convergence with kink.
  //4=Couldn't find seed for patrec. 5=Failed in pat rec filtering
  int _failType;

  size_t nnodes;

  //Temporary fix (??). vector to store hadron unit dir vec.
  EVector _hadunit;
  double _hadEng;
  
  
  // Stuff relevant for event classification, will be uncommented when needed.
  event_classif _classify;

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

class forwardSorter{
public:
  bool operator()(const Measurement* p1, const Measurement* p2){
    if (p2->position()[2] > p1->position()[2]) return true;
    return false;
  }

};
  
#endif
