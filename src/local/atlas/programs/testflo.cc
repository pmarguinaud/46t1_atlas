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

using namespace atlas;
using namespace atlas::grid;
using namespace atlas::mesh;
using namespace atlas::util;
using namespace atlas::array;
using namespace atlas::meshgenerator;
using namespace atlas::output;


static
void showGrid (const StructuredGrid & grid, bool fullPrint)
{
  std::cout << " ny = " << grid.xspace ().ny () << std::endl;
  std::cout << " nxmin = " << grid.xspace ().nxmin () << std::endl;
  std::cout << " nxmax = " << grid.xspace ().nxmax () << std::endl;
  
  std::cout << " xspace = " << grid.xspace ().type () << std::endl;
  std::cout << " nx = ";
//  std::for_each (grid.xspace ().nx ().begin (), grid.xspace ().nx ().end (), 
//                 [&] (idx_t k) { std::cout << k << ", "; });
  std::cout << std::endl;
  
  std::cout << " xmin = " << grid.xspace ().xmin ().front () << std::endl;
  
  std::cout << " xmax = " << grid.xspace ().xmax ().front () << std::endl;
  
  std::cout << " yspace = " << grid.yspace ().type () << std::endl;
  std::cout << " ymin = " << grid.yspace ().min () << std::endl;
  std::cout << " ymax = " << grid.yspace ().max () << std::endl;
  std::cout << " size = " << grid.yspace ().size () << std::endl;

  std::cout << " y = ";
//  std::for_each (grid.yspace ().begin (), grid.yspace ().end (), [&] (float x) { std::cout << x << ", "; });
  std::cout << std::endl;
  if(fullPrint){
    int k = 0, l = 0;
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
      printf (" %8d > %20.10f, %20.10f | %20.10f, %20.10f\n", k, ll.lon (), ll.lat (), xy.x (), xy.y ());
      k++;
    }

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
void swap (short &b)
{
  char * x = (char *)&b;
  std::swap (x[0], x[1]);
}

struct Supergrid {
  StructuredGrid grid;
  Mesh mesh;
  functionspace::StructuredColumns columns;
  functionspace::NodeColumns nodes;
  Field field;
};

static
Supergrid makeSubSuperGrid(const Supergrid supergrid,const int reductionFactor){
  int nx,ny;
  double dx,dy,xmin,xmax,ymin,ymax;
  nx=supergrid.grid.nx ()[0];
  ny=supergrid.grid.ny ();
  dx=supergrid.grid.xspace ().dx()[0];
  dy= (supergrid.grid.yspace ().max()-supergrid.grid.yspace ().min())/supergrid.grid.yspace ().size();
  xmin=supergrid.grid.xspace ().min ();
  xmax=supergrid.grid.xspace ().max ();
  ymin=supergrid.grid.yspace ().min ();
  ymax=supergrid.grid.yspace ().max ();
  if(!(nx%reductionFactor==0)&&!(ny%reductionFactor==0)){
    std::cout <<"Unable to generate subGrid : "<<nx<<" "<<ny<<" "<<reductionFactor<< std::endl;
    exit;
  }


  Config c;
  c.set ("type", "regional");
  c.set ("nx", nx/reductionFactor);
  c.set ("ny", ny/reductionFactor);
  c.set ("dx", dx/reductionFactor);
  c.set ("dy", dy/reductionFactor);
  c.set ("xmin", xmin);
  c.set ("xmax", xmax);
  c.set ("ymax", ymax);
  c.set ("ymin", ymin);
  Grid g = Grid(c);
  showGrid(g,false);
  
  Supergrid  reducedSupergrid;
  Mesh m = MeshGenerator ("structured").generate (g);
  
  
  auto sourceView = array::make_view<float, 1>( supergrid.field );

  Config interpConfig;
//  interpConfig.set( "type", "nearest-neighbour" );
  interpConfig.set( "type", "structured-bilinear" );
  interpConfig.set( "halo", 1 );
  functionspace::StructuredColumns reducedColumns = functionspace::StructuredColumns(g);
  functionspace::NodeColumns reducedNodes = functionspace::NodeColumns(m);
  Field reducedField = reducedNodes.createField<double>(interpConfig );



  std::cout <<"Prepare interp"<< std::endl;
  std::cout <<supergrid.columns.size()<< std::endl;
  std::cout <<reducedColumns.size()<< std::endl;
  Interpolation interpolation_fwd(interpConfig, supergrid.columns, reducedColumns );
  std::cout <<"exec interp"<< std::endl;
  interpolation_fwd.execute( supergrid.field, reducedField );
  reducedField.haloExchange();
  std::cout <<"exchange ok"<< std::endl;


  reducedSupergrid.grid=g;
  reducedSupergrid.mesh=m;
  reducedSupergrid.nodes=reducedNodes;
  reducedSupergrid.columns=reducedColumns;
  reducedSupergrid.field=reducedField;
  


  return reducedSupergrid;
}	

static
Supergrid dirHdrToAtlas (const std::string &dirHdrPath)
{
  std::ifstream ifs_hdr( dirHdrPath+".hdr",std::ifstream::in );
  if ( ! ifs_hdr.is_open() ) {                 
    std::cout <<" Failed to open hdr file: "<<dirHdrPath+".hdr"<< std::endl;
    exit;
  }
  else {
    std::cout <<"Opened OK" << std::endl;
  }

  std::map<std::string, double> doubleMap;  
  std::string line;
  while (std::getline(ifs_hdr, line)) {
    std::size_t pos = line.find(":");	  
    if(pos != std::string::npos){
        std::string key = line.substr (0,pos);
	std::string val = line.substr (pos+1);
	doubleMap.insert ( std::pair<std::string,double>(key,atof(val.c_str())));

    }
  }
  
  Config c;

  c.set ("type", "regional");
  c.set ("nx", int(doubleMap.find("cols")->second));
  c.set ("ny", int(doubleMap.find("rows")->second));
  c.set ("dx", (doubleMap.find("east")->second - doubleMap.find("west")->second)/doubleMap.find("cols")->second);
  c.set ("dy", (doubleMap.find("north")->second - doubleMap.find("south")->second)/doubleMap.find("rows")->second);
  c.set ("xmin", doubleMap.find("west")->second); 
  c.set ("xmax", doubleMap.find("east")->second);
  c.set ("ymax", doubleMap.find("north")->second); 
  c.set ("ymin", doubleMap.find("south")->second);

  c.set ("projection", Config ("type", "rotated_lonlat"));
  StructuredGrid grid = Grid (c);
  ifs_hdr.close();
  std::cout <<"Grid initialized... Reading dir" << std::endl;
  
  std::ifstream ifs( dirHdrPath+".dir",std::ifstream::in );
  if ( ! ifs.is_open() ) {
    std::cout <<" Failed to open dir: "<<dirHdrPath+".dir" << std::endl;
    exit;
  }
  else {
    std::cout <<"Opened OK" << std::endl;
  }
  Mesh mesh = MeshGenerator ("structured").generate (grid);
  auto lonlat = array::make_view<double, 2>( mesh.nodes().xy() );
  int nproc = atlas::mpi::comm ().size ();
  int iproc = atlas::mpi::comm ().rank ();
  int gsize = grid.size ();

  std::cout <<"nproc: "<<nproc << std::endl;
  std::cout <<"iproc: "<<iproc << std::endl;
  std::cout <<"lonlat size: "<<lonlat.size() << std::endl;
  std::cout <<"grid size: "<<gsize << std::endl;
  std::cout <<"mesh size: "<<mesh.partitionGraph().size() << std::endl;

  Config fsConf;
  fsConf.set("halo",2); //needed for pole inter?
  fsConf.set("name","toto");

  functionspace::NodeColumns nodes = functionspace::NodeColumns (mesh,fsConf);
  Field f = nodes.createField<float> ();
  
  const int nx = grid.nx ()[0];
  const int ny = grid.ny ();
  short num[nx];
  int i;
  auto fv = array::make_view<float,1>(f);
  int j;
  std::cout <<"startLoop" << std::endl;

  for( j = 0; j < ny; ++j ) {
    ifs.read((char*)&num, sizeof(num));
    for(i=0;i<nx;++i){
      swap (num[i]);
      fv (j*nx+i) = float(num[i]);
    }
  }

  f.haloExchange();

  functionspace::StructuredColumns columns = functionspace::StructuredColumns (grid,fsConf);

  Supergrid  supergrid;
  supergrid.grid=grid;
  supergrid.mesh=mesh;
  supergrid.nodes=nodes;
  supergrid.columns=columns;
  supergrid.field=f;
  return supergrid;
  
}

double vortex_rollup( double lon, double lat, double t ) {
    // lon and lat in degrees!

    // Formula found in "A Lagrangian Particle Method with Remeshing for Tracer Transport on the Sphere"
    // by Peter Bosler, James Kent, Robert Krasny, Christiane Jablonowski, JCP 2015

    lon *= M_PI / 180.;
    lat *= M_PI / 180.;

    auto sqr           = []( const double x ) { return x * x; };
    auto sech          = []( const double x ) { return 1. / std::cosh( x ); };
    const double T     = 1.;
    const double Omega = 2. * M_PI / T;
    t *= T;
    const double lambda_prime = std::atan2( -std::cos( lon - Omega * t ), std::tan( lat ) );
    const double rho          = 3. * std::sqrt( 1. - sqr( std::cos( lat ) ) * sqr( std::sin( lon - Omega * t ) ) );
    double omega              = 0.;
    double a                  = util::Earth::radius();
    if ( rho != 0. ) {
        omega = 0.5 * 3 * std::sqrt( 3 ) * a * Omega * sqr( sech( rho ) ) * std::tanh( rho ) / rho;
    }
    double q = 1. - std::tanh( 0.2 * rho * std::sin( lambda_prime - omega / a * t ) );
    return q;
};


static 
void pt (const char * label)
{
  static double t0 = -1;
  struct timeval tv;
  struct timezone tz;
  gettimeofday (&tv, &tz);
  double t = tv.tv_sec + tv.tv_usec * 1E-6;
  if (t0 > 0)
    printf ("%-30s %12.4f\n", label, t-t0);
  t0 = t;
}



int main (int argc, char * argv[]) 
{

  pt (nullptr);

  atlas::Library::instance ().initialise (argc, argv);
    Config config;
    config.set( "type", "structured-linear2D" );
    config.set( "halo", 1 );
    int nlev=0;
    std::cout << "sourceGrid creation"<< std::endl;
//    StructuredGrid src_grid("L2160x4320");
    StructuredGrid src_grid("L21600x43200");
    pt ("src_grid");

//    std::cout << "sourceGrid created : footprint:"<<src_grid.footprint()<< std::endl;

    std::cout << "sourceDist creation"<< std::endl;
    grid::Distribution src_dist (src_grid, grid::Partitioner (Config ("type", "checkerboard")));
    pt ("src_dist");
//    std::cout << "sourceDist created : footprint:"<<src_dist.footprint()<< std::endl;
//
        for ( long p = 0; p < mpi::comm().size(); ++p ) {
          std::cout << src_dist.nb_pts()[p]<< std::endl;
        }

//
    std::cout << "targetGrid creation"<< std::endl;
//    StructuredGrid tgt_grid("L5400x10800");
    StructuredGrid tgt_grid("L100x200");
    pt ("tgt_grid");
//    std::cout << "targetGrid created : footprint:"<<tgt_grid.footprint()<< std::endl;


    std::cout << "targetDist creation"<< std::endl;
    grid::Distribution tgt_dist (tgt_grid, grid::Partitioner (Config ("type", "checkerboard")));
    pt ("tgt_dist");
//    std::cout << "targetDist created : footprint:"<<tgt_dist.footprint()<< std::endl;

    functionspace::StructuredColumns src_fs;
    functionspace::StructuredColumns tgt_fs;
  
    std::cout << "sourceFunctionspace creation"<< std::endl;
    src_fs = functionspace::StructuredColumns{src_grid,src_dist, config | option::levels( nlev )};
    std::cout << "sourceFunctionspace created : footprint:"<<src_fs.footprint()<< std::endl;
    pt ("src_fs");

    std::cout << "targetFunctionspace creation"<< std::endl;
    tgt_fs = functionspace::StructuredColumns{tgt_grid,tgt_dist, config | option::levels( nlev )};
    std::cout << "targetFunctionspace created : footprint:"<<tgt_fs.footprint()<< std::endl;
    pt ("tgt_fs");

    std::cout << "sourceField creation"<< std::endl;
    Field src_field = src_fs.createField<double>( option::name( "source" ) );
    std::cout << "sourceField created : footprint:"<<src_field.footprint()<< std::endl;

    std::cout << "targetField creation"<< std::endl;
    Field tgt_field = tgt_fs.createField<double>( option::name( "target" ) );
    std::cout << "tagetField created : footprint:"<<tgt_field.footprint()<< std::endl;


    std::cout << "lonlat creation"<< std::endl;
    auto lonlat = array::make_view<double, 2>( src_fs.xy() );
//    std::cout << "lonlat created : footprint:"<<lonlat.footprint()<< std::endl;

    std::cout << "source creation"<< std::endl;
    auto source = array::make_view<double, 1>( src_field );
//    std::cout << "source created : footprint:"<<source.footprint()<< std::endl;
    idx_t size  = src_fs.size();

    std::cout << " src_fs.size() = " <<  src_fs.size() << std::endl;

    atlas_omp_parallel_for( idx_t n = 0; n < size; ++n ) {
//       source( n ) = vortex_rollup( lonlat( n, LON ), lonlat( n, LAT ),mpi::comm().rank()*100  );
       source( n ) = mpi::comm().rank()*100;
    }




    std::cout << "halo exchange"<< std::endl;
    src_field.haloExchange();


    return 0;

    std::cout << "interp"<< std::endl;
    Interpolation interpolation_fwd(config, src_fs, tgt_fs );
//    std::cout << "interp created : footprint:"<<interpolation_fwd.footprint()<< std::endl;

    interpolation_fwd.execute( src_field, tgt_field );
    std::cout << "interp ended"<< std::endl;
    tgt_field.haloExchange();
    std::cout << "interp halo ended"<< std::endl;
    idx_t tgt_size  = tgt_fs.size();

//    auto dest = array::make_view<double, 1>( tgt_field );
    Mesh mesh = MeshGenerator ("structured").generate (tgt_grid, tgt_dist);

    mpi::comm().barrier();
    Gmsh gmsh( "tgt_.msh" );
    gmsh.write( mesh); 
    gmsh.write( tgt_field); 
  
//    const std::string path="/home/suzat/dirHdr/medium/orographyCompleted"; 
/*    const std::string path="/home/suzat/dirHdr/light/orographyCompleted"; 
    Supergrid supergrid;
    supergrid=dirHdrToAtlas(path);
//  StructuredGrid grid = forgeStructuredGrid ();
//  bool details=false; 
//
/*   
    Config config;
    config.set( "type", "structured-linear2D" );
    config.set( "halo", 1 );
    
    StructuredGrid src_grid("L1080x2160");
    StructuredGrid tgt_grid("L540x1080");
    functionspace::StructuredColumns src_fs;
    functionspace::StructuredColumns tgt_fs;

    grid::Distribution src_dist (src_grid, grid::Partitioner (Config ("type", "checkerboard")));
    grid::Distribution tgt_dist (tgt_grid, grid::Partitioner (Config ("type", "checkerboard")));

    Mesh src_mesh = MeshGenerator ("structured").generate (src_grid, src_dist);
    Mesh tgt_mesh = MeshGenerator ("structured").generate (tgt_grid, tgt_dist);


    int nlev=0;
    src_fs = functionspace::StructuredColumns{src_grid, config | option::levels( nlev )};
    auto tgt_partitioner = grid::MatchingPartitioner(src_fs);
    tgt_fs = functionspace::StructuredColumns{tgt_grid, config | option::levels( nlev )};
    Field src_field = src_fs.createField<double>( option::name( "source" ) );
    Field tgt_field = tgt_fs.createField<double>( option::name( "target" ) );
    auto lonlat = array::make_view<double, 2>( src_fs.xy() );
    auto source = array::make_view<double, 1>( src_field );
    idx_t size  = src_fs.size();
    atlas_omp_parallel_for( idx_t n = 0; n < size; ++n ) {
       source( n ) = vortex_rollup( lonlat( n, LON ), lonlat( n, LAT ), (0.5 + n) / 2 );
    }
    src_field.haloExchange();
    FieldSet src_fset;
    src_fset.add(src_field);
    output::Output src_gmsh;
    output::Output tgt_gmsh;
    src_gmsh = output::Gmsh( "src_field.msh" );
    src_gmsh.write( src_mesh );
    src_gmsh.write( src_fset );
    Interpolation interpolation_fwd(config, src_fs, tgt_fs );
    interpolation_fwd.execute( src_field, tgt_field );
    tgt_field.haloExchange();
    tgt_gmsh = output::Gmsh( "tgt_field.msh" );
    tgt_gmsh.write( tgt_mesh );
    tgt_gmsh.write( tgt_field );
    idx_t tgt_size  = tgt_fs.size();

    auto dest = array::make_view<double, 1>( tgt_field );
    atlas_omp_parallel_for( idx_t n = 0; n < tgt_size; ++n ) {
       std::cout << dest(n) << std::endl;
    }

    mpi::comm().barrier();
*/
  
//    const std::string path="/home/suzat/dirHdr/medium/orographyCompleted"; 
/*    const std::string path="/home/suzat/dirHdr/light/orographyCompleted"; 
    Supergrid supergrid;
    supergrid=dirHdrToAtlas(path);

    showGrid (supergrid.grid,false);
    int reducFactor=4;
    Supergrid subSuperGrid;
    subSuperGrid=makeSubSuperGrid(supergrid,reducFactor);
//    grid::Distribution dist (grid, grid::Partitioner (Config ("type", "checkerboard")));
//    Mesh mesh = MeshGenerator ("structured").generate (grid, dist);

     
    Gmsh gmsh( "src_.msh" );
    gmsh.write( supergrid.mesh ); 
    gmsh.write( supergrid.field); 


  StructuredGrid tgt_grid("L540x1080");
//  Grid grid1 = forgeRegionalGrid ();

//showGrid (grid1);

//  distributeGrid (grid);
*/ 
  atlas::Library::instance ().finalise ();
  return 0;
}
