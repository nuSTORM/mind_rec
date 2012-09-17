/* -*- mode: c++ -*- */
#ifndef _mind_plotter___
#define _mind_plotter___

#include <bhep/event.h>
#include <bhep/particle.h>

#include <TFile.h>
#include <TTree.h>

#include <recpack/RecpackManager.h>
//#include <recpack/Measurement.h>
#include <recpack/Ring.h>

#include <mind/fitter.h>
#include <mind/cluster.h>

using namespace std;
using namespace Recpack;
using namespace bhep;

/* A Class with functions to plot the results of fits to MIND data */

class MINDplotter{

 public:

  MINDplotter();

  virtual ~MINDplotter();

  //Main functions for∫best∫∫ output initialization.
  void initialize(string outFileName,  bool patRec, bool clust, bhep::prlevel vlevel=bhep::NORMAL);

  void execute(fitter& Fit, const bhep::event& evt);///
  void finalize();

  ////
  //To calculate expected vertex position give trajectory and vertex location.
  bool extrap_to_vertex(const Trajectory& traj, 
			const bhep::Point3D& vertexLoc, fitter& fitObj, State& ste, const int trajno);

  //To calculate pos, charge, direction and momentum without extrapolating to vertex //tapasi
  bool fill_kinematics(const Trajectory& traj, State& ste, const int trajno);

  /*Requires a vector of the node of interest, corresponding covariance Matrix
    and the event under study to calculate either position or momentum pulls */
  void position_pulls(const int trajno);
  void momentum_pulls(const int trajno);
  void direction_pulls(const int trajno);

  /* Function to find maximum local chi2 in a particular trajectory
     and enter that value in a histogram */
  void max_local_chi2(const Trajectory& traj, const int trajno);

  /* Function to plot stats about pattern recogntion, 1: no clust, 2: clust */
  void patternStats1(fitter& Fit,const Trajectory& traj, const int trajno ); ///
  void patternStats2(fitter& Fit, const Trajectory& traj, const int trajno);///

  /*Function to record quality of hadron fit*/
  void hadron_direction(fitter& fit);

  /*tru interaction type*/
  void get_tru_int_type(const bhep::event& evt);

protected:

  bhep::prlevel level;
    
  bhep::messenger m;

  TFile *outFile;

  TTree *statTree;
  bool _patR, _clu;

private:

  EVector vert;
  EMatrix vertMat;
  
  bhep::particle* _truPart;
  int _evNo;
  int _Fit[10];
  int _reFit[10];
  int _fail[10];
  int _failEvent;
  int _intType[10];
  double _vert[3];
  double _nuEng;
  int _charm[2];
  double _Q2;
  int _pdg[3];
  double _visEng;
  double _engTraj[10];
  double _engvar[2][10];
  //Information on hadrons.
  double _hadP[3][10];
  double _engTrans;
  int _nhad[2][10];
  double _chadP[4][10];
  double _nhadP[4][10];
  double _hadE[2][10];
  //double _haddot;
  double _X[6][2][10];
  double _Th[6][2][10];

  double _qP[6][10];
  
  double _leng[10];
  double _rangqP[3][10];
  int _Q[3][10];

  double _Chi[2][10];
  int _plns[2][10];
  // std::vector<int> _plns[2];
  
  int _nhits[10];
  int _hitType[5][10];
  

  //std::vector<double>_XPos[10];
  std::vector< vector<double> > _XPos;
  std::vector< vector<double> > _YPos;
  std::vector< vector<double> > _ZPos;
  std::vector< vector<double> > _Edep;

  //double _XPos[10][2500];
  //double _YPos[10][2500];
  //double _ZPos[10][2500];
  //double _Edep[10][2500];


  bool _mus[10][2500];
  bool _cand[10][2500];
  bool _node[10][2500];
  bool _had[10][2500];
  double _pChi[3];
  //TString _intName;
  int _truInt;
  int _reseed_count[10][0];///
  double _Xtent[10];///
  double _vertZ[10];///
  int _traj_hits;
  //For CA
  int _trajs_no; ///
  // std::vector<Trajectory*> _trajs;
 
  // double _chi2[50];
  double _chi2[50];
  int _hits[50];
  int _fitted_hits[50];
  int _traj_index;
  // Trajectory _best_traj;
  //int _nhit;
  //  vector<cluster*> hits ;
  //vector<cluster*> hadhits ;
 
  int _hType[4];
  

  void define_tree_branches();
  //1: Old method. 2: with clustering etc.
  bool extract_true_particle1(const bhep::event& evt, fitter& Fit, const Trajectory& traj, const int trajno);
  bool extract_true_particle2(const bhep::event& evt, fitter& Fit, const Trajectory& traj, const int trajno);
  void add_to_hads(const bhep::particle& part, const int trajno);
  //
  int _truMuHitIndex[2]; 
  int HitIndex; 
};

#endif
