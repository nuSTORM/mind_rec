#include <hit_constructor.h>

hit_constructor::hit_constructor(const bhep::gstore& store)
{
  //Store detector geometry and calculate scint layer positions.

  _detectorLength = store.fetch_dstore("MIND_z") * m;
  _detectorX = store.fetch_dstore("MIND_x") * m;
  _detectorY = store.fetch_dstore("MIND_y") * m;
  _passiveLength = store.fetch_dstore("widthI") * cm;
  _activeLength = store.fetch_dstore("widthS") * cm;
  _gapLength = store.fetch_dstore("widthA") * cm;
  _nActive = store.fetch_istore("nplane");

  _voxXdim = store.fetch_dstore("rec_boxX") * cm;
  _voxYdim = store.fetch_dstore("rec_boxY") * cm;
  _nVoxX = (int)( _detectorX / _voxXdim );
  _nVox = _nVoxX * (int)( _detectorY / _voxYdim );

  double res = store.fetch_dstore("pos_res") * cm;
  _cov = EMatrix(2,2,0);
  for (int i = 0;i < 2;i++)
    _cov[i][i] = pow(res, 2);

  _measType = store.fetch_sstore("meas_type");

  calculate_layerZ();
  
}

hit_constructor::~hit_constructor()
{
}

void hit_constructor::reset()
{
  //Clear out map.correct??
  _voxels.clear();

}

void hit_constructor::execute(const std::vector<bhep::hit*>& hits,
			      std::vector<rec_hit*>& meas)
{
  //First clear out map.
  reset();
  
  //copy hits so they can be sorted in z.
  std::vector<bhep::hit*> sortedHits = hits;

  sort( sortedHits.begin(), sortedHits.end(), forwardSort() );
  
  //sort into voxels map.
  parse_to_map( sortedHits );
  
  //Make rec_hits from vox.
  construct_hits( meas );
  
}

void hit_constructor::calculate_layerZ()
{
  //Fill vector with all possible scint z positions.

  double pieceLength = _passiveLength + _nActive*_activeLength
    + (_nActive+1)*_gapLength;
  int npieces = (int)ceil( _detectorLength / pieceLength );

  //reset detector length to integer multiple of pieces.
  _detectorLength = npieces * pieceLength;

  int nLayers = npieces * _nActive;

  double z;

  for (int iLayer = 0;iLayer < nLayers;iLayer++){

    z = 0;

    for (int j = 0;j < _nActive;j++){

      z += _gapLength + _activeLength;

      if ( (iLayer-j) % _nActive == 0 ){
	int block_no = (iLayer-j) / _nActive;
	z += pieceLength*block_no + _passiveLength
	  - _activeLength/2;
	break;
      }

    }

    z -= _detectorLength/2;
    
    _zLayer.push_back( z );

  }
  
}

void hit_constructor::parse_to_map(const std::vector<bhep::hit*> hits)
{
  //Sort hits into voxel map.
  int voxNo;
  double zpos;

  //Set iterator to plane 1.
  _zIt = _zLayer.begin();

  std::vector<bhep::hit*>::const_iterator hitIt;

  for (hitIt = hits.begin();hitIt != hits.end();hitIt++){

    //find plane.
    zpos = find_plane( *(*hitIt) );
    
    //get Vox number;
    voxNo = calculate_vox_no( *(*hitIt) );

    //Add voxel to map (or hit to existing voxel).

    _voxels[zpos].insert( pair<int,bhep::hit*>(voxNo,(*hitIt)) );

  }

}

double hit_constructor::find_plane(bhep::hit& curHit)
{
  //find the appropriate z position by comparison to
  //Layer z.

  double modDiff = fabs( curHit.x()[2] - (*_zIt) );

  while ( modDiff > _activeLength/2 ){

    _zIt++;

    modDiff = fabs( curHit.x()[2] - (*_zIt) );

  }

  return (*_zIt);
}

int hit_constructor::calculate_vox_no(bhep::hit& curHit)
{
  //calculate the correct voxel number.

  int xbox = (int)( (curHit.x()[0] + _detectorX/2) / _voxXdim );

  int ybox = (int)( fabs( curHit.x()[1] - _detectorY/2 ) / _voxYdim );

  int vox_num = xbox + ybox*_nVoxX;

  return vox_num;
}

void hit_constructor::construct_hits(std::vector<rec_hit*>& meas)
{
  //takes the voxels which have been filled and make
  //rec_hit objects out of them.

  std::map<double,std::multimap<int,bhep::hit*> >::iterator vIt;

  for (vIt = _voxels.begin();vIt != _voxels.end();vIt++){

    for (int ivox = 0;ivox < _nVox;ivox++){

      if ( vIt->second.count(ivox) != 0 ){

	rec_hit* vhit = get_vhit( ivox, vIt->first, vIt->second );

	meas.push_back( vhit );

      }
    }
  }

}

rec_hit* hit_constructor::get_vhit(int vox, double z,
				   const std::multimap<int,bhep::hit*>& map)
{
  //Makes a rec_hit from the voxel position and adds the relevant points.
  
  int irow = vox / _nVoxX;
  int icol = vox % _nVoxX;
  //std::cout << "Vox to XY test: " << vox << ", " <<irow << ", " << icol << std::endl;
  double voxX = icol*_voxXdim + _voxXdim/2 - _detectorX/2;
  double voxY = _detectorY/2 - (irow*_voxYdim + _voxYdim/2);
  //std::cout << "Values: " << voxX << ", " << voxY << std::endl;
  EVector hit_pos(2,0);
  hit_pos[0] = voxX;
  hit_pos[1] = voxY;

  EVector meas_pos(3,0);
  meas_pos[0] = voxX;
  meas_pos[1] = voxY;
  meas_pos[2] = z;
  
  rec_hit* me = new rec_hit();
  me->set_name(_measType);
  me->set_hv(HyperVector(hit_pos,_cov));
  me->set_name("volume", "Detector");
  me->set_position( meas_pos );

  me->set_vox_no( vox );
  //add the hits to the rec_hit.
  std::multimap<int,bhep::hit*>::const_iterator hIt;
  for (hIt = map.equal_range(vox).first;hIt != map.equal_range(vox).second;hIt++)
    {
      me->add_point( (*hIt).second );
    }
  
  return me;
}
