#include "atlas/library/Library.h"
#include "atlas/util/Config.h"
#include "atlas/runtime/Log.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/array.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/functionspace.h"


#include "fa.h"


#include <iostream>
#include <algorithm>
#include <vector>

using namespace atlas;
using namespace atlas::grid;
using namespace atlas::mesh;
using namespace atlas::util;
using namespace atlas::array;
using namespace atlas::meshgenerator;

StructuredGrid forgeGrid ()
{
#include "t32.h"

  std::vector<Spacing> spacings (N);

  for (int i = 0; i < N; i++)
    spacings[i] = Spacing (Config ("type", "linear") | Config ("N", pl[i])
                         | Config ("start", 0) | Config ("end", 360));

  StructuredGrid::XSpace xspace (spacings);
  StructuredGrid::YSpace yspace (Config ("type", "gaussian") | Config ("N", N));
  Projection proj (Config ("type", "rotated_schmidt") | Config ("stretching_factor", stretchingFactor)
                 | Config ("north_pole", std::vector<double>{longitudeOfStretchingPoleInDegrees, latitudeOfStretchingPoleInDegrees}));

  return StructuredGrid (xspace, yspace, proj);
}

static
void distributeGrid (const Grid & grid)
{
  int nproc = atlas::mpi::comm ().size (); 
  int iproc = atlas::mpi::comm ().rank (); 

  grid::Distribution dist (grid, grid::Partitioner (Config ("type", "equal_regions")));

  const std::vector<idx_t> & nb_pts = dist.nb_pts (); 

  Mesh mesh = MeshGenerator ("structured").generate (grid, dist);
  int lsize = mesh.nodes ().size (); 
  int gsize = grid.size (); 


  std::cout << " partition = " << mesh.partition () << std::endl;
  std::cout << " nproc = " << nproc << std::endl;
  std::cout << " iproc = " << iproc << std::endl;
  std::cout << " lsize = " << lsize << std::endl;
  std::cout << " gsize = " << gsize << std::endl;
  
  functionspace::NodeColumns nodes_fs = functionspace::NodeColumns (mesh);

  {
    printf ("%s:%d\n", __FILE__, __LINE__);
    const auto & remote_index = array::make_view<idx_t, 1>( mesh.nodes().remote_index() );  
    printf (" remote_index \n");
    for (int i = 0; i < 10; i++)
      printf (" %8d > %8d\n", i, remote_index (i));
  }

  Field f = nodes_fs.createField<double> (option::name ("field"));

  FieldSet lfs;
  lfs.add (f);

  {
    auto fv = array::make_view<double,1>(f);
    for (int i = 0; i < fv.size (); i++)
      fv (i) = double (iproc);
  }

  FieldSet gfs;
  
  const int iproc0 = 1; // Gather on task #1
  Field g ("global", array::DataType ("real64"), make_shape (iproc == iproc0 ? grid.size () : 0));
  g.metadata ().set ("owner", iproc0);

  gfs.add (g);

  auto gv = array::make_view<double,1>(g);
  for (int i = 0; i < gv.size (); i++)
    gv (i) = -1; 

  nodes_fs.gather (lfs, gfs);


  std::cout << "gv" << std::endl;
  for (int i = 0; i < gv.size (); i++)
    std::cout << i << " " << gv (i) << std::endl;

}

int main (int argc, char * argv[])
{
  
  character * CLNOMF = argv[1];
  
  integer64 IREP, INUMER = 0, INIMES = 2, INBARP = 0, INBARI = 0;
  logical LLNOMM = fort_TRUE, LLERFA = fort_TRUE, LLIMST = fort_TRUE;
  character CLNOMC[2] = "c";

  faitou64_ (&IREP, &INUMER, &LLNOMM, CLNOMF, (character*)"OLD", &LLERFA, &LLIMST, &INIMES, &INBARP, &INBARI,      
             CLNOMC, strlen (CLNOMF), 3, strlen (CLNOMC));

  const int JPLATMAX = 2000;
  real64 ZSLAPO, ZCLOPO, ZSLOPO, ZCODIL, ZREFER;
  integer64 ITYPTR, ITRONC, INLATI, INXLON, INIVER;
  integer64 INLOPA[JPLATMAX], INOZPA[JPLATMAX];
  real64 ZSINLA[JPLATMAX], ZAHYBR[JPLATMAX], ZBHYBR[JPLATMAX];
  logical LLGARD = fort_FALSE;

  facies64_ (CLNOMC, &ITYPTR, &ZSLAPO, &ZCLOPO, &ZSLOPO,  &ZCODIL, &ITRONC, &INLATI, &INXLON, &INLOPA[0],       
             &INOZPA[0], &ZSINLA[0], &INIVER, &ZREFER, &ZAHYBR[0], &ZBHYBR[0], &LLGARD, strlen (CLNOMC));

  for (int i = 0; i < INLATI / 2; i++)
    printf (" %8d > %8ld\n", i, INLOPA[i]);

  fairme64_ (&IREP, &INUMER, "KEEP", 4);

  return 0;

  atlas::Library::instance ().initialise (argc, argv);

  StructuredGrid grid = forgeGrid (); 

  distributeGrid (grid);
  
  atlas::Library::instance ().finalise (); 
  return 0;
}