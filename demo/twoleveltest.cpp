#include "ddmsolver.hpp"
#include <fstream>
#include <iostream>
//#define pointsource
#define randomrhs
using namespace std;
using namespace mfem;
void E_exact(const Vector &, Vector &);
void f_exact(const Vector &, Vector &);
void testu(const Vector &, Vector &);
double freq = 6.0, kappa=freq*M_PI;
double sigmaconst=kappa*kappa;
int dim;
int main(int argc, char *argv[])
{
    //1. Read the mesh & set options
   //const char *mesh_file = "../../mfem/mfem-4.7//data/beam-tet.mesh";
   const char *mesh_file = " cubemesh.msh";
   int order = 1;
   bool static_cond = false;
   bool pa = false;
   const char *device_config = "cpu";
   bool visualization = 1;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
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
   int dim = mesh->Dimension();
   int sdim = mesh->SpaceDimension();
   int nsubdomain=4;
   int ref_levels =0;
   int delta=1;
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
   // cout<<"boundary dofs are "<<endl;
   // for(int i=0;i<ess_tdof_list.Size();i++)
   // {
   //    cout<<ess_tdof_list[i]<<"  ";
   // }

   //4.define linear form of RHS
   VectorFunctionCoefficient f(sdim, f_exact);
   LinearForm *b = new LinearForm(fespace);
   b->AddDomainIntegrator(new VectorFEDomainLFIntegrator(f));
   b->Assemble();
   GridFunction x(fespace);
   VectorFunctionCoefficient E(sdim, E_exact);
   x.ProjectCoefficient(E);
   //5. define bilinear form
   Coefficient *muinv = new ConstantCoefficient(1.0);
   Coefficient *sigma = new ConstantCoefficient(sigmaconst);
   BilinearForm *a = new BilinearForm(fespace);
   a->AddDomainIntegrator(new CurlCurlIntegrator(*muinv));
   a->AddDomainIntegrator(new VectorFEMassIntegrator(*sigma));
   //6. Assemble the linear system
   a->Assemble();
   SparseMatrix A_c;
   Vector B_c, X_c;
   a->FormLinearSystem(ess_tdof_list, x, *b, A_c, X_c, B_c);
   cout << "Size of linear system: " << A_c.Height() << endl;
   DDMSolver::ASMSolver *ASM_c=new DDMSolver::ASMSolver(nsubdomain, A_c.Height());
   ASM_c->Setup(*mesh,A_c,*fespace,0);
   //7. assemble in fine mesh
   //set the FEMHierarchy
   FiniteElementSpaceHierarchy FEMH(mesh,fespace,false,false);
   FEMH.AddUniformlyRefinedLevel();
   Operator *P=FEMH.GetProlongationAtLevel(0);
   cout<<"multigrid has level "<<FEMH.GetNumLevels()<<endl;
   // Operator* P=FEMH.GetProlongationAtLevel(0);
   cout<<"the number of dof in fine mesh is "<<(FEMH.GetFinestFESpace()).GetNDofs()<<endl;
   cout<<"the prologation operator is "<<P->Width()<<" times "<<P->Height()<<endl;
   //set boundary
   FiniteElementSpace &finefespace=FEMH.GetFinestFESpace();
   Array<int> ess_bdr(mesh->bdr_attributes.Max()),f_ess_tdof_list;
   ess_bdr=1;
   finefespace.GetEssentialTrueDofs(ess_bdr,f_ess_tdof_list);
   //set bilinearform
   BilinearForm *a_f = new BilinearForm(&finefespace);
   a_f->AddDomainIntegrator(new CurlCurlIntegrator(*muinv));
   a_f->AddDomainIntegrator(new VectorFEMassIntegrator(*sigma));
   a_f->Assemble();
   //set linearform
   LinearForm *b_f = new LinearForm(&finefespace);
   b_f->AddDomainIntegrator(new VectorFEDomainLFIntegrator(f));
   b_f->Assemble();
   //form linear system
   SparseMatrix A_f;
   Vector B_f, X_f;
   GridFunction x_f(&finefespace);
   x_f.ProjectCoefficient(E);
   cout<<"Discretized Error is "<<x_f.ComputeL2Error(E)<<endl;
   a_f->FormLinearSystem(f_ess_tdof_list, x_f, *b_f, A_f, X_f, B_f);
   Mesh *mesh_f=finefespace.GetMesh();
   DDMSolver::ASMSolver *ASM_f=new DDMSolver::ASMSolver(8*nsubdomain,A_f.Height());
   ASM_f->Setup(*mesh_f,A_f,finefespace,delta);
   ASM_f->setUnityVector(A_f.Height());
   A_c.Finalize();
   ASM_f->setCoarseSpace(P,A_c);
   #ifdef randomrhs
   B_f.Randomize();
   #endif
   cout<<"norm of B is "<<B_f.Norml2()<<endl;
   GMRES(A_f, *ASM_f, B_f, X_f, 1, 1000, 1000, 1e-12, 0.0);
   a_f->RecoverFEMSolution(X_f, *b_f, x_f);
   #ifndef pointsource
   cout << "\n|| E_h - E ||_{L^2} = " << x_f.ComputeL2Error(E) << '\n' << endl;
   #endif
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
   f(2) = exp(-60*(x(0)-4)*(x(0)-4)+(x(1)-0.5)*(x(1)-0.5)+(x(2)-0.5)*(x(2)-0.5));
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
void testu(const Vector &x, Vector &f)
{
   if (dim == 3)
   {
      f(0) =  x(2);
      f(1) = x(2);
      f(2) = 1-x(1)-x(0);
   }
   else
   {
      f(0) = (sigmaconst + kappa * kappa) * sin(kappa * x(1));
      f(1) = (sigmaconst + kappa * kappa) * sin(kappa * x(0));
      if (x.Size() == 3) { f(2) = 0.0; }
   }
}
#endif
