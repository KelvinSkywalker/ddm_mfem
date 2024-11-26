#include "ddmsolver.hpp"
#include <fstream>
#include <iostream>
#include <memory>
//#define usetwolevel
//#define pointsource
//#define randomrhs
//#define ASMDir
using namespace std;
using namespace mfem;
void E_exact(const Vector &, Vector &);
void minusE_exact(const Vector &, Vector &);
void ietaE_exact(const Vector &, Vector &);
void f_exact(const Vector &, Vector &);
void UnitNormal(const Vector &,Vector &);
int dim;
double freq=10.0;
double kappa = freq*M_PI;
double sigmaconst=-kappa*kappa;
double eta=kappa;
int main(int argc, char *argv[]){
	//1. set the mesh
	int order = 1;
	bool herm_conv = true;
	int ref_levels =2; // if use twolevel, ref_level > = 1
	int num_subdomain=8;
	int delta=1;
	#ifdef usetwolevel
	Mesh *mesh_c = new Mesh("../../../mfem/mfem-4.7//data/beam-tet.mesh", 1, 1);
	mesh_c->FinalizeMesh();
	dim=mesh_c->Dimension();
	for (int l = 0; l < ref_levels-1; l++)
	{
		mesh_c->UniformRefinement();
	}
	FiniteElementCollection *fec_c = new ND_FECollection(order, dim);
   	FiniteElementSpace *fespace_c = new FiniteElementSpace(mesh_c, fec_c);
	FiniteElementSpaceHierarchy FEMH(mesh_c,fespace_c,false,false);
	FEMH.AddUniformlyRefinedLevel();
	Operator *P=FEMH.GetProlongationAtLevel(0);
	cout<<"the prologation operator is "<<P->Width()<<" times "<<P->Height()<<endl;
   	FiniteElementSpace &fespacetemplate=FEMH.GetFinestFESpace();
	FiniteElementSpace *fespace=new FiniteElementSpace(fespacetemplate);
	Mesh *mesh=fespace->GetMesh();
	#else
	//Mesh *mesh = new Mesh("../../mfem/mfem-4.7//data/beam-tet.mesh", 1, 1);
	Mesh *mesh= new Mesh("../demo_mesh/cube.msh",1,1);
    mesh->FinalizeMesh();
	dim=mesh->Dimension();
	for(int i=0;i<mesh->GetNBE();i++)
	{
		mesh->SetBdrAttribute(i,1);
	}
   for(int i=0;i<mesh->GetNE();i++)
	{
		mesh->SetAttribute(i,1);
	}
	for (int l = 0; l < ref_levels; l++)
	{
		mesh->UniformRefinement();
	}
	//2. set the Finite Element Space
	FiniteElementCollection *fec;
   	FiniteElementSpace *fespace;
	fec = new ND_FECollection(order, dim);
    fespace = new FiniteElementSpace(mesh, fec);
	cout<<"number of dof is "<<	fespace->GetNDofs()<<endl;
	#endif
	ofstream outfile;
	outfile.open("mesh.vtk");
	mesh->PrintVTK(outfile);
	outfile.close();
	//3. set up boundary condition
	Array<int> ess_tdof_list;
	Array<int> ess_bdr;
	Array<int> rbmarker;
	rbmarker.SetSize(mesh->bdr_attributes.Max());
	rbmarker = 1;
	rbmarker[0] = 0;
    ess_bdr.SetSize(mesh->bdr_attributes.Max());
	ess_bdr = 0;
	fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
	cout<<"essential dof is "<<ess_tdof_list.Size()<<endl;
 	// 4. Setup Complex Operator convention
 	ComplexOperator::Convention conv =
 		herm_conv ? ComplexOperator::HERMITIAN : ComplexOperator::BLOCK_SYMMETRIC;
	//5. Define the solution vector x and rhs f as a complex finite element grid function
    //     corresponding to fespace. Define some vector function used by BC
    ComplexGridFunction x(fespace);
	GridFunction *minusu=new GridFunction(fespace);
	Vector zero_vector(3);
	zero_vector=0.0;
	ConstantCoefficient *minusone=new ConstantCoefficient(-1.0);
	VectorConstantCoefficient zero_v(zero_vector);
	VectorFunctionCoefficient f(dim, f_exact);
	VectorFunctionCoefficient E(dim,E_exact);
	#ifndef pointsource
	VectorFunctionCoefficient minusE(dim,E_exact,minusone);
	VectorFunctionCoefficient ietaE(dim,ietaE_exact);
	VectorFunctionCoefficient n(dim,UnitNormal);
	VectorCrossProductCoefficient ietaEtimesn(n,ietaE);
	minusu->ProjectCoefficient(minusE);
	CurlGridFunctionCoefficient curlE(minusu);
	x.real().ProjectCoefficient(E);
   	cout<<"the interpolation error is "<< x.real().ComputeL2Error(E)<<endl;
	#endif
	//6. define the linear form
 	ComplexLinearForm *b = new ComplexLinearForm(fespace, conv);
	b->AddDomainIntegrator(new VectorFEDomainLFIntegrator(f),new VectorFEDomainLFIntegrator(zero_v));
 	#ifndef pointsource
	b->AddBoundaryIntegrator(new VectorFEBoundaryTangentLFIntegrator(curlE),new VectorFEBoundaryTangentLFIntegrator(ietaEtimesn));
	#endif
	b->Assemble();
	//7. define the bilinearform
	// Integrators inside the computational domain
	ConstantCoefficient zero(0.0);
	ConstantCoefficient one(1.0);
	ConstantCoefficient nksqr(sigmaconst);
	ConstantCoefficient ksi(0.0);
	ConstantCoefficient ik(-eta); 
	SesquilinearForm *a = new SesquilinearForm(fespace, conv);
	a->AddDomainIntegrator(new CurlCurlIntegrator(one), new CurlCurlIntegrator(zero));
	a->AddDomainIntegrator(new VectorFEMassIntegrator(nksqr),new VectorFEMassIntegrator(ksi));
 	a->AddBoundaryIntegrator(new VectorFEMassIntegrator(zero), new VectorFEMassIntegrator(ik));
 	a->Assemble(0);
	//8. Assemble the linear system
 	OperatorPtr A;
 	Vector B, X;
 	a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);
	ComplexSparseMatrix sparseA=*A.As<ComplexSparseMatrix>();
	sparseA.real().Finalize();
	sparseA.imag().Finalize();
	//9.set up the ASM preconditioner
	#ifdef ASMDir
	DDMSolver::ASMSolver_ComplexDir ASM_Dir(num_subdomain,fespace->GetNDofs());
	ASM_Dir.Setup(*mesh,*fespace,delta,&a->real());
	#else
	DDMSolver::ASMSolver_Impedence ASM_Imp(num_subdomain,2*fespace->GetNDofs());
	ASM_Imp.Setup(*mesh,*fespace,delta,a,kappa,eta);
	ASM_Imp.setUnityVector(fespace->GetNDofs());
	#endif
	#ifdef usetwolevel
	//9.5 Assemble the coarsematrix
	SesquilinearForm *a_c = new SesquilinearForm(fespace_c, conv);
	a_c->AddDomainIntegrator(new CurlCurlIntegrator(one), new CurlCurlIntegrator(zero));
	a_c->AddDomainIntegrator(new VectorFEMassIntegrator(nksqr),new VectorFEMassIntegrator(ksi));
 	a_c->AddBoundaryIntegrator(new VectorFEMassIntegrator(zero), new VectorFEMassIntegrator(ik));
 	a_c->Assemble(0);
	OperatorPtr A_c;
	a_c->FormSystemMatrix(ess_tdof_list,A_c);
	ComplexSparseMatrix sparseA_c=*A_c.As<ComplexSparseMatrix>();
	#ifdef ASMDir
	ASM_Dir.setCoarseSpace(P,sparseA_c);
	#else
	ASM_Imp.setCoarseSpace(P,sparseA_c);
	#endif
	#endif
	#ifdef randomrhs
	B.Randomize();
	#endif
	#ifdef ASMDir
	GMRES(*A, ASM_Dir, B, X, 1, 1000, 1000, 1e-12, 0.0);
	#else
	GMRES(*A, ASM_Imp, B, X, 1, 1000, 1000, 1e-12, 0.0);
	#endif
	//10.solve the linear system
#ifdef MFEM_USE_SUITESPARSE
    // ComplexUMFPackSolver csolver(sparseA);
	// csolver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
    // csolver.SetPrintLevel(1);
    // csolver.Mult(B, X);
#else
	Array<int> offsets(3);
      	offsets[0] = 0;
      	offsets[1] = fespace->GetTrueVSize();
      	offsets[2] = fespace->GetTrueVSize();
      	offsets.PartialSum();
	std::unique_ptr<Operator> pc_r;
      	std::unique_ptr<Operator> pc_i;
      	double s = (conv == ComplexOperator::HERMITIAN) ? -1.0 : 1.0;

	BilinearForm prec(fespace);
	ConstantCoefficient ksqr(kappa * kappa);
	prec.AddDomainIntegrator(new CurlCurlIntegrator(one));
	prec.AddDomainIntegrator(new VectorFEMassIntegrator(ksqr));
	prec.Assemble();
        OperatorPtr PCOpAh;
        prec.SetDiagonalPolicy(mfem::Operator::DIAG_ONE);
        prec.FormSystemMatrix(ess_tdof_list, PCOpAh);
        // Gauss-Seidel Smoother
        pc_r.reset(new GSSmoother(*PCOpAh.As<SparseMatrix>()));
        pc_i.reset(new ScaledOperator(pc_r.get(), s));
	cout << "set up block diagnal"<<endl;
	BlockDiagonalPreconditioner BlockDP(offsets);
	BlockDP.SetDiagonalBlock(0, pc_r.get());
	BlockDP.SetDiagonalBlock(1, pc_i.get());
	cout<<"BlockDP has size"<<BlockDP.Height()<<endl;
 	cout << "start gmres"<<endl;	
	GMRESSolver gmres;
	gmres.SetPrintLevel(1);
	gmres.SetKDim(200);
	gmres.SetMaxIter(2000);
	gmres.SetRelTol(1e-4);
	gmres.SetAbsTol(1e-4);
	gmres.SetOperator(*A);
	gmres.SetPreconditioner(BlockDP);
	
	cout << "gmres iterations: "<<endl;
	gmres.Mult(B, X);
	cout << "gmres stopped. " << endl;
#endif
	// 11. Recover the solution as a finite element grid function and compute error. 
	//the solutioon is store and can be viewed using GLVis: 
	//"glvis -m cavity.mesh -g cavity_sol.gf".
	a->RecoverFEMSolution(X, *b, x);
	cout << "\n|| E_h - E ||_{L^2} = " << x.real().ComputeL2Error(E)<< '\n' << endl;
	{
		ofstream mesh_ofs("cavity.mesh");
		mesh_ofs.precision(8);
		mesh->Print(mesh_ofs);

		ofstream sol_r_ofs("cavity_sol.gf");
		sol_r_ofs.precision(8);
		x.real().Save(sol_r_ofs);
	}

	// 14. Send the solution by socket to a GLVis server.
	char vishost[] = "localhost";
	int  visport   = 19916;
	socketstream sol_sock(vishost, visport);
	sol_sock.precision(8);
	sol_sock << "solution\n" << *mesh << x.real() << "keys amrRljcUUuuu\n" << flush;
	#ifdef usetwolevel
	if(fec_c)
	{
		delete fec_c;
		delete fespace_c;
		delete fespace;
	}
	#else
	if (fec)
   	{
	delete fespace;
	delete fec;
	delete mesh;
	}
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
void ietaE_exact(const Vector &x, Vector &ikE)
{
    ikE(0) = eta*sin(kappa * x(1));
    ikE(1) = eta*sin(kappa * x(2));
    ikE(2) = eta*sin(kappa * x(0));
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
void UnitNormal(const Vector &x,Vector &f)
{
	if(x(1)==0)
	{
		f(0)=0;
		f(1)=-1;
		f(2)=0;
	}
	else if(x(1)==1)
	{
		f(0)=0;
		f(1)=1;
		f(2)=0;
	}
	else if(x(0)==0)
	{
		f(0)=-1;
		f(1)=0;
		f(2)=0;
	}
	else if(x(0)==1)
	{
		f(0)=1;
		f(1)=0;
		f(2)=0;
	}
	else if(x(2)==0)
	{
		f(0)=0;
		f(1)=0;
		f(2)=-1;
	}
	else if(x(2)==1)
	{
		f(0)=0;
		f(1)=0;
		f(2)=1;
	}
	else
	{
		f=0.0;
	}
}