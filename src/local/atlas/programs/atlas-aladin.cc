#include "atlas/library/Library.h"
#include "atlas/util/Config.h"
#include "atlas/runtime/Log.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/array.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/functionspace.h"


#include <iostream>
#include <algorithm>
#include <vector>

using namespace atlas;
using namespace atlas::grid;
using namespace atlas::mesh;
using namespace atlas::util;
using namespace atlas::array;
using namespace atlas::meshgenerator;


static
void showGrid (const StructuredGrid & grid)
{
  std::cout << " ny = " << grid.xspace ().ny () << std::endl;
  std::cout << " nxmin = " << grid.xspace ().nxmin () << std::endl;
  std::cout << " nxmax = " << grid.xspace ().nxmax () << std::endl;
  
  std::cout << " xspace = " << grid.xspace ().type () << std::endl;
  std::cout << " nx = ";
  std::for_each (grid.xspace ().nx ().begin (), grid.xspace ().nx ().end (), 
                 [&] (idx_t k) { std::cout << k << ", "; });
  std::cout << std::endl;
  
  std::cout << " xmin = " << grid.xspace ().xmin ().front () << std::endl;
  
  std::cout << " xmax = " << grid.xspace ().xmax ().front () << std::endl;
  
  std::cout << " yspace = " << grid.yspace ().type () << std::endl;
  std::cout << " ymin = " << grid.yspace ().min () << std::endl;
  std::cout << " ymax = " << grid.yspace ().max () << std::endl;
  std::cout << " size = " << grid.yspace ().size () << std::endl;

  std::cout << " y = ";
  std::for_each (grid.yspace ().begin (), grid.yspace ().end (), [&] (float x) { std::cout << x << ", "; });
  std::cout << std::endl;

  int k = 0;
  for (int iy = 0; iy < grid.ny (); iy++)
    {
      for (int ix = 0; ix < grid.nx (iy); ix++)
        {
          PointXY p = grid.xy (ix, iy);
          PointLonLat ll = grid.lonlat (ix, iy);
          printf (" %8d > %20.10f, %20.10f | %20.10f, %20.10f\n", k++, ll.lon (), ll.lat (), p.x (), p.y ());
        }
    }
  printf ("...\n");

}

static
void showGrid (const Grid & grid)
{

  const Projection & proj = grid.projection ();

  std::cout << " size = " << grid.size () << std::endl;

  int k = 0;
  for (const auto & xy : grid.xy ())
    {
      const PointLonLat ll = proj.lonlat (xy);
      printf (" %20.10f, %20.10f\n", ll.lon (), ll.lat ());
      k++;
    }

  const Config config = grid.spec ();

  std::cout << config << std::endl;

}

static
StructuredGrid forgeStructuredGrid ()
{
#include "aro0064x0064.h"
//#include "aro1440x1536.h"

  std::vector<Spacing> spacings (Ny);

  for (int i = 0; i < Ny; i++)
    spacings[i] = Spacing (Config ("type", "linear") | Config ("N", Nx)
                         | Config ("start", -Nux / 2 * DxInMetres) | Config ("end", (Nx - 1 - Nux / 2) * DxInMetres));

  StructuredGrid::XSpace xspace (spacings);
  StructuredGrid::YSpace yspace (Config ("type", "linear") | Config ("N", Ny) | Config ("start", -Nuy / 2 * DyInMetres) | Config ("end", (Ny - 1 - Nuy / 2) * DyInMetres));
  Projection proj (Config ("type", "lambert_conformal_conic") | Config ("longitude0", LoVInDegrees)
                 | Config ("latitude0", LaDInDegrees) | Config ("latitude1", Latin1InDegrees) 
                 | Config ("latitude2", Latin2InDegrees));

  return StructuredGrid (xspace, yspace, proj);
}

static
Grid forgeRegionalGrid ()
{
#include "aro0064x0064.h"
//#include "aro1440x1536.h"

  Config c;

  c.set ("type", "regional");
  c.set ("nx", Nx);
  c.set ("ny", Ny);
  c.set ("dx", DxInMetres);
  c.set ("dy", DyInMetres); 
  c.set ("xmin", -Nux / 2 * DxInMetres); c.set ("xmax", (Nx - 1 - Nux / 2) * DxInMetres);
  c.set ("ymax", -Nuy / 2 * DyInMetres); c.set ("ymin", (Ny - 1 - Nuy / 2) * DyInMetres);

  c.set ("projection", Config ("type", "lambert_conformal_conic") | Config ("longitude0", LoVInDegrees)
                     | Config ("latitude0", LaDInDegrees) | Config ("latitude1", Latin1InDegrees) 
                     | Config ("latitude2", Latin2InDegrees));
  
  return Grid (c);
}

static
void distributeGrid (const Grid & grid)
{
  int nproc = atlas::mpi::comm ().size ();
  int iproc = atlas::mpi::comm ().rank ();

  grid::Distribution dist (grid, grid::Partitioner (Config ("type", "checkerboard")));

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

  Field f = nodes_fs.createField<double> (option::name ("field"));

  {
    auto fv = array::make_view<double,1>(f);
    for (int i = 0; i < fv.size (); i++)
      fv (i) = double (iproc);
  }

  FieldSet lfs;
  lfs.add (f);

  FieldSet gfs;
  
  const int iproc0 = 1; // Gather on task #1
  Field g ("regional", array::DataType ("real64"), make_shape (iproc == iproc0 ? grid.size () : 0));
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
  atlas::Library::instance ().initialise (argc, argv);

  StructuredGrid grid = forgeStructuredGrid ();
  
  showGrid (static_cast<Grid> (grid));

  Grid grid1 = forgeRegionalGrid ();

//showGrid (grid1);

//distributeGrid (grid);
  
  atlas::Library::instance ().finalise ();
  return 0;
}
