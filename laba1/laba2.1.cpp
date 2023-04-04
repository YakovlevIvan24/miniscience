// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 13
//
//  Remeshing an STL file without an underlying CAD model
//
// -----------------------------------------------------------------------------
#define _USE_MATH_DEFINES
#include <set>
#include <cmath>
#include <gmsh.h>

int main(int argc, char **argv)
{
  gmsh::initialize();

  gmsh::model::add("ROCKET");

  // Let's merge an STL mesh that we would like to remesh (from the parent
  // directory):
  try {
    gmsh::merge("ROCKET4.stl");
  } catch(...) {
    gmsh::logger::write("Could not load STL mesh: bye!");
    gmsh::finalize();
    return 0;
  }
    gmsh::option::setNumber("Mesh.Algorithm3D",1);
  // We first classify ("color") the surfaces by splitting the original surface
  // along sharp geometrical features. This will create new discrete surfaces,
  // curves and points.

    gmsh::model::mesh::classifySurfaces(0 * M_PI / 180., true, true);

    gmsh::model::mesh::createGeometry();
  // Angle between two triangles above which an edge is considered as sharp:


  // For complex geometries, patches can be too complex, too elongated or too
  // large to be parametrized; setting the following option will force the
  // creation of patches that are amenable to reparametrization

  // Create a geometry for all the discrete curves and surfaces in the mesh, b

  // Create a volume from all the surfaces
  std::vector<std::pair<int, int> > s;
  gmsh::model::getEntities(s, 2);
  std::vector<int> sl;
  for(auto surf : s) sl.push_back(surf.second);
  int l = gmsh::model::geo::addSurfaceLoop(sl);
  gmsh::model::geo::addVolume({l});

  gmsh::model::geo::synchronize();

  // We specify element sizes imposed by a size field, just because we can :-)
  gmsh::model::mesh::generate(3);
    gmsh::write("FINAL_ROCKET.msh");
  // Launch the GUI to see the results:
  std::set<std::string> args(argv, argv + argc);
  if(!args.count("-nopopup")) gmsh::fltk::run();

  gmsh::finalize();
  return 0;
}
