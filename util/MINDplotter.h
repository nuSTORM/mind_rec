/* -*- mode: c++ -*- */
#ifndef _mind_plotter___
#define _mind_plotter___

#include <bhep/event.h>
#include <bhep/particle.h>

#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TProfile.h>
#include <TTree.h>

#include <recpack/RecpackManager.h>
#include <recpack/Measurement.h>
#include <recpack/Ring.h>
#include <mind/fitter.h>

using namespace std;
using namespace Recpack;
using namespace bhep;

/* A Class with functions to plot the results of fits to MIND data */

class MINDplotter{

 public:

  MINDplotter();

  virtual ~MINDplotter();

  //Main functions for output initialization.
  bool initialize(TString outFileName, bhep::prlevel vlevel=bhep::NORMAL);
  bool execute(fitter& Fit, const bhep::event& evt, bool success, bool patRec);
  bool finalize();
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

  /* Function to plot stats about pattern recogntion */
  void patternStats(fitter& Fit);

  /*Function to record types of fit and quality of seed*/
  void track_fit(fitter& fit);

  /*Function to record quality of hadron fit*/
  void hadron_direction(fitter& fit);

protected:

  bhep::prlevel level;
    
  bhep::messenger m;

  TFile *outFile;

  TTree *statTree;

private:

  EVector vert;
  EMatrix vertMat;
  
  bhep::particle* _truPart;
  int _evNo;
  bool _Fit;
  int _fail;
  double _fitTracker[5];
  double _nuEng;
  double _hadP[3];
  double _haddot;
  double _X[3][2];
  double _Th[3][2];
  double _qP[3];
  int _Q[3];
  double _Chi[2];
  int _nhits;
  int _hitType[4];
  double _XPos[300];
  double _YPos[300];
  double _ZPos[300];
  bool _mus[300];
  bool _cand[300];
  bool _node[300];
  double _pChi[3];

  void define_tree_branches();

  bool extract_true_particle(const bhep::event& evt, fitter& Fit, bool patRec);

};

#endif
