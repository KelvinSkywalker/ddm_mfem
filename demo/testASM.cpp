#include "ddmsolver.hpp"
#include <fstream>
#include <iostream>
#define pointsource
//#define randomrhs
using namespace std;
using namespace mfem;
void E_exact(const Vector &, Vector &);
void f_exact(const Vector &, Vector &);
double freq = 3.0, kappa=freq*M_PI;
double sigmaconst=kappa*kappa;
int dim;
int main(int argc, char *argv[])
{
    //1. Read the mesh & set options
   //const char *mesh_file = "../../mfem/mfem-4.7//data/beam-tet.mesh";
   const char *mesh_file = "../cpp/demo_mesh/cube.msh";
   int order = 1;
   bool static_cond = false;
   bool pa = false;
   const char *device_config = "cpu";
   bool visualization = true;

   OptionsParser args(argc, argv);
   // args.AddOption(&mesh_file, "-m", "--mesh",
   //                "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&freq, "-f", "--frequency", "Set the frequency for the exact"
                  " solution.");
   args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
                  "--no-static-condensation", "Enable static condensation.");
   args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
                  "--no-partial-assembly", "Enable Partial Assembly.");
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);
   Device device(device_config);
   device.Print();
   mfem::Mesh *mesh=new mfem::Mesh(mesh_file, 1, 1);
   mesh->Finalize();
   int dim = mesh->Dimension();
   int sdim = mesh->SpaceDimension();
   int ref_levels =4;
   int nsubdomain=512;
   int delta=1;
   for(int i=0;i<mesh->GetNBE();i++)
   {
      mesh->SetBdrAttribute(i,2);
   }
   for(int i=0;i<mesh->GetNE();i++)
   {
      mesh->SetAttribute(i,1);
   }
         // (int)floor(log(50000./mesh.GetNE())/log(2.)/dim);
   for (int l = 0; l < ref_levels; l++)
      {
         mesh->UniformRefinement();
      }
   ofstream outfile;
	outfile.open("mesh.vtk");
	mesh->PrintVTK(outfile);
	outfile.close();
   //2. Define FEM Space,here we use Nedlec FEM
   FiniteElementCollection *fec = new ND_FECollection(order, dim);
   FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);
   cout<<"the number of Dof is "<<fespace->GetNDofs()<<endl;
   //3. Classify the boundary dof
   Array<int> ess_tdof_list;
   if (mesh->bdr_attributes.Size())
   {
      Array<int> ess_bdr(mesh->bdr_attributes.Max());
      ess_bdr = 1;
      fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }
   //4.define linear form of RHS
   VectorFunctionCoefficient f(sdim, f_exact);
   LinearForm *b = new LinearForm(fespace);
   b->AddDomainIntegrator(new VectorFEDomainLFIntegrator(f));
   b->Assemble();
   GridFunction x(fespace);
   VectorFunctionCoefficient E(sdim, E_exact);
   x.ProjectCoefficient(E);
   cout<<"the interpolation error is "<< x.ComputeL2Error(E)<<endl;
   //5. define bilinear form
   Coefficient *muinv = new ConstantCoefficient(1.0);
   Coefficient *sigma = new ConstantCoefficient(sigmaconst);
   BilinearForm *a = new BilinearForm(fespace);
   a->SetDiagonalPolicy(mfem::AbstractSparseMatrix::DiagonalPolicy::DIAG_ONE);
   a->AddDomainIntegrator(new CurlCurlIntegrator(*muinv));
   a->AddDomainIntegrator(new VectorFEMassIntegrator(*sigma));
   //6. Assemble the linear system
   a->Assemble();
   SparseMatrix A;
   Vector B, X;
   a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);
   cout << "Size of linear system: " << A.Height() << endl;
   cout<<"norm of B is "<<B.Norml2()<<endl;
   //7. solve by DDM
   DDMSolver::ASMSolver ASMtest(nsubdomain, A.Height());
   ASMtest.Setup(*mesh,A,*fespace,delta);
   //ASMtest.Setup(*mesh,*fespace,delta,a);
   ASMtest.setUnityVector(fespace->GetNDofs());
   X=0.0;
   #ifdef randomrhs
   B.Randomize();
   #endif
   // UMFPackSolver umfsol(true);
   // umfsol.SetOperator(A);
   // umfsol.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
   // umfsol.Mult(B,X);
   GMRES(A, ASMtest, B, X, 1, 2000, 2000, 1e-12, 0.0);
   a->RecoverFEMSolution(X, *b, x);
   cout << "\n|| E_h - E ||_{L^2} = " << x.ComputeL2Error(E) << '\n' << endl;
   //6. Visualize the solution and mesh
   if (visualization)
   {
      ofstream mesh_ofs("refined.mesh");
      mesh_ofs.precision(8);
      mesh->Print(mesh_ofs);
      ofstream sol_ofs("sol.gf");
      sol_ofs.precision(8);
      x.Save(sol_ofs);
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock.precision(8);
      sol_sock << "solution\n" << mesh << x << flush;
   }
   return 0;
}
#ifdef pointsource
void E_exact(const Vector &x, Vector &E)
{
   E=0.0;
}
void f_exact(const Vector &x, Vector &f)
{
   f(0) = 0;
   f(1) = 0;
   f(2) = exp(-60*(x(0)-0.5)*(x(0)-0.5)+(x(1)-0.5)*(x(1)-0.5)+(x(2)-0.5)*(x(2)-0.5));
}
#else
void E_exact(const Vector &x, Vector &E)
{
   if (dim == 3)
   {
      E(0) = sin(kappa * x(1));
      E(1) = sin(kappa * x(2));
      E(2) = sin(kappa * x(0));
   }
   else
   {
      E(0) = sin(kappa * x(1));
      E(1) = sin(kappa * x(0));
      if (x.Size() == 3) { E(2) = 0.0; }
   }
}

void f_exact(const Vector &x, Vector &f)
{
   if (dim == 3)
   {
      f(0) = (sigmaconst + kappa * kappa) * sin(kappa * x(1));
      f(1) = (sigmaconst + kappa * kappa) * sin(kappa * x(2));
      f(2) = (sigmaconst + kappa * kappa) * sin(kappa * x(0));
   }
   else
   {
      f(0) = (sigmaconst + kappa * kappa) * sin(kappa * x(1));
      f(1) = (sigmaconst + kappa * kappa) * sin(kappa * x(0));
      if (x.Size() == 3) { f(2) = 0.0; }
   }
}
#endif
