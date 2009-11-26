#include <hit_clusterer.h>

#include <CLHEP/Random/RandGauss.h>

#include <map>

hit_clusterer::hit_clusterer(const bhep::gstore& store)
{
  //get relevant information from the store.
  long seed = (long)store.fetch_dstore("Gen_seed");
  _ranGen = RanluxEngine( seed, 4 );

  _sigMa = store.fetch_dstore("pos_sig") * cm;

  _vInX = (int)( (store.fetch_dstore("MIND_x")*m) / (store.fetch_dstore("rec_boxX")*cm) );

  _cov = EMatrix(2,2,0);
  for (int i = 0;i < 2;i++)
    _cov[i][i] = pow( _sigMa, 2 );

  _measType = store.fetch_sstore("meas_type");
  
}

hit_clusterer::~hit_clusterer()
{
}

void hit_clusterer::execute(const std::vector<bhep::hit*>& deps,
			    std::vector<cluster*>& clusts)
{
  //take a vector of voxel hits and make clusters.
  double zplane;
  std::vector<bhep::hit*>::const_iterator depIt;
  std::map<int,bhep::hit*> plane_hits;

  zplane = (*deps.begin())->x()[2];
  for (depIt = deps.begin();depIt != deps.end();depIt++){

    if ( (*depIt)->x()[2] != zplane ){

      //have all hits in this plane so add the clusters to vector.
      clusters_in_plane( zplane, plane_hits, clusts );

      plane_hits.clear();
      
      zplane = (*depIt)->x()[2];
      plane_hits.insert( pair<int,bhep::hit*>((*depIt)->idata("voxel"),(*depIt)) );

    } else plane_hits.insert( pair<int,bhep::hit*>((*depIt)->idata("voxel"),(*depIt)) );

  }
  //And the last plane.
  if ( plane_hits.size() != 0 ){
    clusters_in_plane( zplane, plane_hits, clusts );
    plane_hits.clear();
  }

}

void hit_clusterer::clusters_in_plane(double zpos, std::map<int,bhep::hit*>& deps,
				      std::vector<cluster*>& clusts)
{
  //if only a single voxel is hit in the plane make the 'cluster' directly.
  if ( deps.size() == 1 ){

    cluster* clst = make_cluster( zpos, deps.begin()->second );

    clusts.push_back( clst );

  } else {
    //Need to find out how many clusters there are and there positions.
    form_clusters( zpos, deps, clusts );
  }

}

void hit_clusterer::form_clusters(double zpos, std::map<int,bhep::hit*>& deps,
				  std::vector<cluster*>& clusts)
{
  //Have to group voxels.
  int vox1 = deps.begin()->first;
  int max_vox;

  std::map<int,bhep::hit*>::iterator mIt;

  std::map<int,double> clhit;
  clhit.insert( pair<int,double>(vox1, deps.begin()->second->fetch_dproperty("TotalEng")) );
  std::map<int,double>::iterator dIt;

  int isearch[] = {_vInX+1,_vInX,_vInX-1,1,-1,-(_vInX-1),-_vInX,-(_vInX+1)};
  
  while ( deps.size() != 0 ) {
 
    for (int ineigh = 0;ineigh < 8;ineigh++){

      mIt = deps.find( vox1 - isearch[ineigh] );

      if ( mIt != deps.end() )
	clhit.insert( pair<int,double>(vox1 - isearch[ineigh],
				       (*mIt).second->fetch_dproperty("TotalEng")) );

    }
    if ( clhit.size() == 1){
      cluster* cl1 = make_cluster( zpos, deps[clhit.begin()->first] );
      clusts.push_back( cl1 );
      deps.erase( clhit.begin()->first );
      clhit.clear();

      if ( deps.size() != 0 ){
	vox1 = deps.begin()->first;
	clhit.insert( pair<int,double>(vox1,
				       deps.begin()->second->fetch_dproperty("TotalEng")) );
      }

    } else {
      max_vox = get_max_vox( clhit );
      if ( max_vox == vox1 ){
	EVector pos(3,0);
	pos[2] = zpos;
	std::vector<bhep::hit*> hits;
	for (dIt = clhit.begin();dIt != clhit.end();dIt++)
	  hits.push_back( deps[dIt->first] );
	calculate_clust_pos( hits, pos );
	cluster* cl2 = make_cluster( pos, hits );
	clusts.push_back( cl2 );
	for (dIt = clhit.begin();dIt != clhit.end();dIt++)
	  deps.erase( dIt->first );
	clhit.clear();

	if ( deps.size() != 0 ){
	  vox1 = deps.begin()->first;
	  clhit.insert( pair<int,double>(vox1,
					 deps.begin()->second->fetch_dproperty("TotalEng")) );
	}
      } else {
	vox1 = max_vox;
	clhit.clear();
	clhit.insert( pair<int,double>( max_vox, 
					deps[max_vox]->fetch_dproperty("TotalEng")) );
      }
    }
      
  }

}

int hit_clusterer::get_max_vox(const std::map<int,double>& voxes)
{
  int max;
  double max_val;
  std::map<int,double>::const_iterator vIt2;

  max = voxes.begin()->first;
  max_val = voxes.begin()->second;

  for (vIt2 = voxes.begin();vIt2 != voxes.end();vIt2++){

    if ( vIt2->second > max_val ){
      max = vIt2->first;
      max_val = vIt2->second;
    }

  }
  
  return max;
}

cluster* hit_clusterer::make_cluster(double zpos, bhep::hit* dep)
{

  EVector hit_pos(2,0);
  hit_pos[0] = dep->x()[0];
  hit_pos[1] = dep->x()[1];
  
  EVector meas_pos(3,0);
  meas_pos[0] = hit_pos[0];
  meas_pos[1] = hit_pos[1];
  meas_pos[2] = zpos;

  cluster* me = new cluster();
  me->set_name(_measType);
  me->set_hv(HyperVector(hit_pos,_cov));
  me->set_name("volume", "Detector");
  me->set_position( meas_pos );

  me->add_hit( dep );

  return me;
}

cluster* hit_clusterer::make_cluster(const EVector& vec,
				     const std::vector<bhep::hit*>& deps)
{

  EVector hit_pos(2,0);
  hit_pos[0] = vec[0];
  hit_pos[1] = vec[1];

  cluster* me = new cluster();
  me->set_name(_measType);
  me->set_hv(HyperVector(hit_pos,_cov));
  me->set_name("volume", "Detector");
  me->set_position( vec );

  std::vector<bhep::hit*>::const_iterator depIt2;
  for (depIt2 = deps.begin();depIt2 != deps.end();depIt2++)
    me->add_hit( (*depIt2) );

  return me;
}

void hit_clusterer::calculate_clust_pos(const std::vector<bhep::hit*>& hits, EVector& vec)
{
  //Weighted mean for x and y position.
  double X = 0, Y = 0, Q = 0, Q1;

  std::vector<bhep::hit*>::const_iterator hitIt;
  for (hitIt = hits.begin();hitIt != hits.end();hitIt++){
    
    Q1 = (*hitIt)->ddata("TotalEng");

    X += Q1 * (*hitIt)->x()[0];
    Y += Q1 * (*hitIt)->x()[1];
    Q += Q1;

  }
    
  vec[0] = X / Q + RandGauss::shoot(&_ranGen, 0, _sigMa);
  vec[1] = Y / Q + RandGauss::shoot(&_ranGen, 0, _sigMa);

}