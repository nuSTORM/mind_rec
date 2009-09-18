/*******************************************************
 *                                                     *
 * Singleton class containing the recpack manager used *
 * to fit and perform pattern recognition.             *
 *                                                     *
 * Andrew Laing, 09/09.                                *
 *                                                     *
 *******************************************************/

#ifndef __FITMAN__
#define __FITMAN__

#include <recpack/RecpackManager.h>
#include <recpack/KalmanFitter.h>
#include <recpack/CellularAutomaton.h>

#include <bhep/gstore.h>

class MINDfitman
{
 public:

  static MINDfitman& instance();

  RecpackManager& manager(){ return _man; }

  void set_man_parameters(const bhep::gstore params, Setup& det);

  void fit_mode();

  void rec_mode();

 private:
  //Constructors, destructor hidden
  MINDfitman();
  ~MINDfitman(){}
  MINDfitman(MINDfitman const&);
  MINDfitman& operator=(MINDfitman const&);

  void config_common_props();
  void config_geometry(Setup& det);
  void config_fitter();
  void config_navigator();
  void config_model();
  void config_matching();

 private:

  RecpackManager _man;

  bhep::gstore _store;

  //-------------- verbosity levels ------------//

  Messenger::Level l0, l1, l2, l3;

  string _model;
  string _fitter;

};

#endif
