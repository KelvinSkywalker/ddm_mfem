#include <iostream>
#include "mfem.hpp"
#include "ddmsolver.hpp"
using namespace std;
using namespace mfem;

void TestSplitMesh()
{
	const char *mesh_file = "./demo_mesh/cube.msh";
	Mesh mesh(mesh_file,1,1);

	int num_node = mesh.GetNV();
	int num_domain = 2;
	int ncommon = 3;

	DDMSolver::ASMSolver asmsolver(num_domain,num_node);
	asmsolver.SplitMesh(mesh,ncommon,true);
}

void TestUMFPack()
{
   const char *mesh_file = "../../mfem/mfem-4.6/data/star.mesh";
   int order = 1;
   bool static_cond = false;
   bool pa = false;
   bool fa = false;

   // 3. Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
   //    the same code.
   Mesh mesh(mesh_file, 1, 1);
   int dim = mesh.Dimension();

   // 4. Refine the mesh to increase the resolution. In this example we do
   //    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
   //    largest number that gives a final mesh with no more than 50,000
   //    elements.
   {
      int ref_levels =
         (int)floor(log(50000./mesh.GetNE())/log(2.)/dim);
      for (int l = 0; l < ref_levels; l++)
      {
         mesh.UniformRefinement();
      }
   }

   // 5. Define a finite element space on the mesh. Here we use continuous
   //    Lagrange finite elements of the specified order. If order < 1, we
   //    instead use an isoparametric/isogeometric space.
   FiniteElementCollection *fec;
   bool delete_fec;
   if (order > 0)
   {
      fec = new H1_FECollection(order, dim);
      delete_fec = true;
   }
   else if (mesh.GetNodes())
   {
      fec = mesh.GetNodes()->OwnFEC();
      delete_fec = false;
      cout << "Using isoparametric FEs: " << fec->Name() << endl;
   }
   else
   {
      fec = new H1_FECollection(order = 1, dim);
      delete_fec = true;
   }
   FiniteElementSpace fespace(&mesh, fec);
   cout << "Number of finite element unknowns: "
        << fespace.GetTrueVSize() << endl;

   // 6. Determine the list of true (i.e. conforming) essential boundary dofs.
   //    In this example, the boundary conditions are defined by marking all
   //    the boundary attributes from the mesh as essential (Dirichlet) and
   //    converting them to a list of true dofs.
   Array<int> ess_tdof_list;
   if (mesh.bdr_attributes.Size())
   {
      Array<int> ess_bdr(mesh.bdr_attributes.Max());
      ess_bdr = 1;
      fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }

   // 7. Set up the linear form b(.) which corresponds to the right-hand side of
   //    the FEM linear system, which in this case is (1,phi_i) where phi_i are
   //    the basis functions in the finite element fespace.
   LinearForm b(&fespace);
   ConstantCoefficient one(1.0);
   b.AddDomainIntegrator(new DomainLFIntegrator(one));
   b.Assemble();

   // 8. Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   GridFunction x(&fespace);
   x = 0.0;

   // 9. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the Laplacian operator -Delta, by adding the Diffusion
   //    domain integrator.
   BilinearForm a(&fespace);
   if (pa) { a.SetAssemblyLevel(AssemblyLevel::PARTIAL); }
   if (fa)
   {
      a.SetAssemblyLevel(AssemblyLevel::FULL);
      // Sort the matrix column indices when running on GPU or with OpenMP (i.e.
      // when Device::IsEnabled() returns true). This makes the results
      // bit-for-bit deterministic at the cost of somewhat longer run time.
      a.EnableSparseMatrixSorting(Device::IsEnabled());
   }
   a.AddDomainIntegrator(new DiffusionIntegrator(one));

   // 10. Assemble the bilinear form and the corresponding linear system,
   //     applying any necessary transformations such as: eliminating boundary
   //     conditions, applying conforming constraints for non-conforming AMR,
   //     static condensation, etc.
   if (static_cond) { a.EnableStaticCondensation(); }
   a.Assemble();

   OperatorPtr A;
   Vector B, X;
   a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

   cout << "Size of linear system: " << A->Height() << endl;

   //=====================================================================
   // int nrow = X.Size(); 
   int nrow = 100;
   cout << "The value of X: ";
   for (int i=0;i<nrow;i++)
	   cout << X[i] << " " ;
   cout << endl;

   cout << "begin to solve " << endl;
   UMFPackSolver umf_solver;
   umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
   umf_solver.SetOperator(*A);
   umf_solver.Mult(B, X);

   cout << "After solving, the value of X: ";
   for (int i=0;i<nrow;i++)
	   cout << X[i] << " " ;
   cout << endl;

   if (delete_fec)
   {
      delete fec;
   }

}

int main()
{
	/* TestSplitMesh(); */

	TestUMFPack();

	return 0;
}
