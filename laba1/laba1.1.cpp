
// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 19
//
//  Thrusections, fillets, pipes, mesh size from curvature
//
// -----------------------------------------------------------------------------

// The OpenCASCADE geometry kernel supports several useful features for solid
// modelling.

#include <set>
#include <cmath>
#include <cstdlib>
#include <gmsh.h>

int main(int argc, char **argv)
{
    gmsh::initialize(argc, argv);

    gmsh::model::add("t19");
     gmsh::model::occ::addTorus(0,0,0,5,3,1);
     gmsh::model::occ::addTorus(0,0,0,5,1.5,2);

    std::vector<std::pair<int, int> > ov;
    std::vector<std::vector<std::pair<int, int> > > ovv;
    gmsh::model::occ::cut({{3, 1}}, {{3, 2}}, ov, ovv, 3);



    gmsh::model::occ::synchronize();


    gmsh::option::setNumber("Mesh.MeshSizeMin", 0.3);
    gmsh::option::setNumber("Mesh.MeshSizeMax", 1);

    gmsh::model::mesh::generate(3);
    gmsh::write("t19.msh");

    // Launch the GUI to see the results:
    std::set<std::string> args(argv, argv + argc);
    if(!args.count("-nopopup")) gmsh::fltk::run();

    gmsh::finalize();
    return 0;
}