#include <dolfin.h>
#include "Poisson.h"


using namespace dolfin;

// Source term (right-hand side)
class Source : public Expression
{
  void eval(Array<double>& values, const Array<double>& x) const
  {
    double dx = x[0] - 0.5;
    double dy = x[1] - 0.5;
    values[0] = 10*exp(-(dx*dx + dy*dy) / 0.02);
  }
};

// Normal derivative (Neumann boundary condition)
class dUdN : public Expression
{
  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = sin(5*x[0]) + cos(5*x[1]);
  }
};

class DirichletValue : public Expression
{
    void eval(Array<double>& values, const Array<double>& x) const
    {
        values[0] = sin(x[0]);
    }
};

// Sub domain for Dirichlet boundary condition
class DirichletBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return x[0] < DOLFIN_EPS  or x[0] > 0.5 - DOLFIN_EPS;
  }
};

int main()
{


//    MeshEditor editor;
//    editor.init_vertices(2500);
//    editor.init_cells(1250);
//    for(int i = 0; i < 50;i++)
//    {
//        for(int j = 0; j < 50;j++)
//        {
//            editor.add_vertex(i+j*50,float(i)/50,float(j)/50);
//        }
//    }
//    for(int i = 0; i < 1250;i++)
//    {
//        editor.add_cell(i,)
//    }

    // Create mesh and function space

    auto mesh = std::make_shared<Mesh>("circle.xml");
    auto V = std::make_shared<Poisson::FunctionSpace>(mesh);

  // Define boundary condition
    auto u0 = std::make_shared<DirichletValue>();
    auto boundary =  std::make_shared<DirichletBoundary>();
    DirichletBC bc(V,u0,boundary);
  // Define variational forms
    Poisson::BilinearForm a(V, V);
    Poisson::LinearForm L(V);
    auto f = std::make_shared<Source>();
    auto g = std::make_shared<dUdN>();
    L.f = f;
    L.g = g;

  // Compute solution
    Function u(V);
    solve(a == L, u, bc);

  // Save solution in VTK format
    File file("poisson.pvd");
    file << u;

    return 0;
}
