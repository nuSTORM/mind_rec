/******************************************************
 *                                                    *
 * Class definition of the hit_reconstructor which    *
 * will create the rec_hit's from the deposits        *
 * from the simulation.                               *
 *                                                    *
 * author: Andrew Laing, 11/2009                      *
 *                                                    *
 ******************************************************/

#ifndef _HIT_CONSTRUCTOR__
#define _HIT_CONSTRUCTOR__

#include <recpack/Measurement.h>

#include <bhep/gstore.h>
#include <bhep/hit.h>

#include <mind/rec_hit.h>

#include <map>

class hit_constructor
{
 public:
  //constructor
  hit_constructor(const bhep::gstore& store);
  //destructor
  ~hit_constructor();

  //reconstruction.
  void execute(const std::vector<bhep::hit*>& hits, measurement_vector& meas);

 private:

  //reset
  void reset();

  //calculate z position of layers.
  void calculate_layerZ();
  //function which puts hits into the map.
  void parse_to_map(const std::vector<bhep::hit*> hits);
  //find plane of hit.
  double find_plane(bhep::hit& curHit);
  //find vox number of hit.
  int calculate_vox_no(bhep::hit& curHit);
  //Make the rec hits.
  void construct_hits(measurement_vector& meas);
  //make an individual rec hit.
  rec_hit* get_vhit(int vox, double z, const std::multimap<int,bhep::hit*>& map);
  

  //vector of plane z positions.
  std::vector<double> _zLayer;
  //iterator for z planes.
  std::vector<double>::iterator _zIt;

  //Parameters related to current MIND setup.
  double _detectorLength;
  double _detectorX;
  double _detectorY;
  double _passiveLength;
  double _activeLength;
  double _gapLength;
  int _nActive;
  
  EMatrix _cov;

  string _measType;

  //Voxel properties.
  double _voxXdim;
  double _voxYdim;
  int _nVoxX;
  int _nVox;

  //Container wchich will define voxels.
  std::map<double, std::multimap<int, bhep::hit*> > _voxels;
  
};

class forwardSort{
 public:
  bool operator()(const bhep::hit* p1, const bhep::hit* p2){
    if (p2->x()[2] > p1->x()[2]) return true;
    return false;
  }

};

#endif
