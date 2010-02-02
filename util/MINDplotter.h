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

/* A Class with functions to plot the results of fits to MIND data */

class MINDplotter{

 public:

  MINDplotter();

  virtual ~MINDplotter();

  //Main functions for output initialization.
  void initialize(string outFileName, bool patRec, bool clust, bhep::prlevel vlevel=bhep::NORMAL);
  void execute(fitter& Fit, const bhep::event& evt, bool success);
  void finalize();
//   //
  //To calculate expected vertex position give trajectory and vertex location.
  bool extrap_to_vertex(const Trajectory& traj, 
			const bhep::Point3D& vertexLoc, fitter& fitObj, State& ste);

  /*Requires a vector of the node of interest, corresponding covariance Matrix
    and the event under study to calculate either position or momentum pulls */
  void position_pulls();
  void momentum_pulls();
  void direction_pulls();

  /* Function to find maximum local chi2 in a particular trajectory
     and enter that value in a histogram */
  void max_local_chi2(const Trajectory& traj);

  /* Function to plot stats about pattern recogntion, 1: no clust, 2: clust */
  void patternStats1(fitter& Fit);
  void patternStats2(fitter& Fit);

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
  bool _Fit;
  bool _reFit;
  int _fail;
  int _intType;
  double _nuEng;
  double _visEng;
  double _engTraj;
  double _engvar[2];
  double _hadP[3];
  double _hadE[2];
  //double _haddot;
  double _X[3][2];
  double _Th[3][2];
  double _qP[3];
  //double _leng;
  //double _rangP[2];
  int _Q[3];
  double _Chi[2];
  int _plns[2];
  int _nhits;
  int _hitType[4];
  double _XPos[2500];
  double _YPos[2500];
  double _ZPos[2500];
  double _Edep[2500];
  bool _mus[2500];
  bool _cand[2500];
  bool _node[2500];
  double _pChi[3];
  //TString _intName;
  int _truInt;

  void define_tree_branches();
  //1: Old method. 2: with clustering etc.
  bool extract_true_particle1(const bhep::event& evt, fitter& Fit);
  bool extract_true_particle2(const bhep::event& evt, fitter& Fit);
  void add_to_hads(const bhep::particle& part);
  //
};

#endif
