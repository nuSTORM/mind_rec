/******************************************************
 *                                                    *
 * Class definition of cluster for hits formed        *
 * from many reconstructed hits in MIND               *
 *                                                    *
 * author: Andrew Laing, 11/2009                      *
 *                                                    *
 ******************************************************/

#ifndef _CLUSTER__
#define _CLUSTER__

#include <bhep/hit.h>

#include <recpack/Measurement.h>

class cluster: public Measurement
{
 public:
  //Constructor.
  cluster();
  //destructor.
  ~cluster();

  //add a contibuting hit.
  void add_hit(bhep::hit* dep);

  //Getters.
  double get_eng(){ return _eng; }
  int get_nhits(){ return _nhit; }
  size_t get_nVox(){ return _nVox; }
  double get_mu_prop(){ return _muProp / (double)_nhit; }
  std::vector<bhep::hit*> get_hits(){ return _voxes; }
  
 private:

  std::vector<bhep::hit*> _voxes;

  size_t _nVox;

  int _nhit;

  double _eng;

  double _muProp;

};

#endif
