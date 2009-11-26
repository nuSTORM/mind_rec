#include <cluster.h>

cluster::cluster() : Measurement()
{
  _eng = 0;
  _muProp = 0;
  _nhit = 0;
}

cluster::~cluster()
{
}

void cluster::add_hit(bhep::hit* dep)
{

  _voxes.push_back( dep );

  _eng += dep->ddata( "TotalEng" );

  _muProp += dep->ddata( "MuonProportion" );
  
  _nVox = _voxes.size();

  _nhit += dep->idata( "NoPoints" );

}
