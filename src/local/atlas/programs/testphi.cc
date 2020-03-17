#include "atlas/library/Library.h"
#include "atlas/util/Config.h"
#include "atlas/runtime/Log.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/array.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/functionspace.h"
#include "atlas/output/Gmsh.h"
#include "atlas/runtime/Trace.h"

#include "atlas/util/CoordinateEnums.h" //to have LON LAT
#include "atlas/interpolation.h"


#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>
#include <string>  

#include <sys/time.h>
#include <stdio.h>

using namespace atlas;
using namespace atlas::grid;
using namespace atlas::mesh;
using namespace atlas::util;
using namespace atlas::array;
using namespace atlas::meshgenerator;
using namespace atlas::output;

static 
void pt (const char * label)
{
  int myproc = atlas::mpi::comm ().rank (); 

  static double t0 = -1;
  struct timeval tv;
  struct timezone tz;
  gettimeofday (&tv, &tz);
  double t = tv.tv_sec + tv.tv_usec * 1E-6;
  if ((t0 > 0) && (myproc == 0))
    printf ("%-30s %12.4f\n", label, t-t0);
  t0 = t;
}



int main (int argc, char * argv[]) 
{
  if (argc < 4)
    {
      printf ("Usage %s: verbose NX NY\n", argv[0]);
      return 0;
    }


  bool verbose = atoi (argv[1]);
  const int NX = atoi (argv[2]), NY = atoi (argv[3]);

  atlas::Library::instance ().initialise (argc, argv);

  int myproc = atlas::mpi::comm ().rank (); 
  
  pt (nullptr);

  StructuredGrid grid (std::string ("L") + std::to_string (NX) + "x" + std::to_string (NY));

  pt ("grid");

  grid::Distribution dist (grid, Config ("light", true) | Config ("blocksize", NX));

  pt ("dist");

if(verbose)
  for (int j = 0, g = 0; j < NY; j++)
  for (int i = 0; i < NX; i++, g++)
    {
      int p = dist.partition (g);
      printf (" %8d -> %8d\n", g, p);
    }
  
  functionspace::StructuredColumns fs {grid, dist, Config ("halo", 1)};

  pt ("fs");

  int k = 0;

if(verbose)
  for (const auto & xy : grid.xy ())
    {
      const PointLonLat ll = grid.projection ().lonlat (xy);
      printf (" %8d > %20.10f, %20.10f | %20.10f, %20.10f\n", k, ll.lon (), ll.lat (), xy.x (), xy.y ());
      k++;
    }

  auto lonlat       = array::make_view<double, 2> (fs.xy ());
  auto global_index = array::make_view<  long, 1> (fs.global_index ());
  auto partition    = array::make_view<   int, 1> (fs.partition ());

  Field field = fs.createField<double> (option::name ("myproc"));

  pt ("field");

if(verbose)
  printf (" fs.size () = %d\n", fs.size ());

if(verbose)
  for (int i = 0; i < fs.size (); i++)
    printf (" %8d > %8d | %c | %20.10f, %20.10f\n", 
            global_index (i)-1, partition (i),  
            partition (i) == myproc ? ' ' : 'G', 
            lonlat (i, LON), lonlat (i, LAT));

  auto view = array::make_view<double, 1> (field);

  for (int i = 0; i < fs.size (); i++)
    if (partition (i) == myproc)
      view (i) = double (myproc);
    else
      view (i) = -1.0;

  field.haloExchange ();
  pt ("halo");

if(verbose)
  printf (" -- field --\n");

if(verbose)
  for (int i = 0; i < fs.size (); i++)
    printf (" %8d > %8d | %c | %20.10f\n", 
            global_index (i)-1, partition (i),  
            partition (i) == myproc ? ' ' : 'G', 
            view (i));

  atlas::Library::instance ().finalise ();

  return 0;
}

