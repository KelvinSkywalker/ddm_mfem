#include "ddmsolver.hpp"
#include <fstream>
#include <iostream>
#include <memory>
using namespace std;
using namespace mfem;
void rb_rhs_bottom(const Vector &x, vector<complex<double>> &fv);
void rb_rhs_top(const Vector &x, vector<complex<double>> &fv);
void rb_rhs_bottom_re(const Vector &x, Vector &v);
void rb_rhs_bottom_im(const Vector &x, Vector &v);
void rb_rhs_top_re(const Vector &x, Vector &v);
void rb_rhs_top_im(const Vector &x, Vector &v);
void rb_rhs_left_re(const Vector &x, Vector &v);
void rb_rhs_left_im(const Vector &x, Vector &v);
int dim;
double k = 10.0;
int main(int argc, char *argv[]){
	int order = 1;
	bool herm_conv = true;
	Mesh *mesh = new Mesh("mycavity_k_10.msh", 1, 1);
   mesh->FinalizeMesh();
	int dim = mesh->Dimension();
	FiniteElementCollection *fec;
   	FiniteElementSpace *fespace;
	fec = new ND_FECollection(order, dim);
    fespace = new FiniteElementSpace(mesh, fec);
	cout<<"the mesh has faces "<<mesh->GetNFaces()<<endl;
	int size = fespace->GetTrueVSize();
	
	Array<int> ess_tdof_list;
	Array<int> ess_bdr;
	Array<int> rbmarker_left, rbmarker_bottom, rbmarker_top, rbmarker;
	rbmarker_left.SetSize(mesh->bdr_attributes.Max());
	rbmarker_left = 0;
	rbmarker_left[4] = 1;
	rbmarker_bottom.SetSize(mesh->bdr_attributes.Max());
        rbmarker_bottom = 0;
        rbmarker_bottom[1] = 1;
	rbmarker_top.SetSize(mesh->bdr_attributes.Max());
        rbmarker_top = 0;
        rbmarker_top[3] = 1;
	rbmarker.SetSize(mesh->bdr_attributes.Max());
	rbmarker = 1;
	rbmarker[0] = 0;

      	ess_bdr.SetSize(mesh->bdr_attributes.Max());
	ess_bdr = 0;
	ess_bdr[0] = 1;
	fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

	VectorFunctionCoefficient f_left_re(dim, rb_rhs_left_re);
	VectorFunctionCoefficient f_left_im(dim, rb_rhs_left_im);
	VectorFunctionCoefficient f_bottom_re(dim, rb_rhs_bottom_re);
	VectorFunctionCoefficient f_bottom_im(dim, rb_rhs_bottom_im);
	VectorFunctionCoefficient f_top_re(dim, rb_rhs_top_re);
	VectorFunctionCoefficient f_top_im(dim, rb_rhs_top_im);
 	// 9. Setup Complex Operator convention
 	ComplexOperator::Convention conv =
 		herm_conv ? ComplexOperator::HERMITIAN : ComplexOperator::BLOCK_SYMMETRIC;
 	ComplexLinearForm *b = new ComplexLinearForm(fespace, conv);
   	b->AddBoundaryIntegrator(new BoundaryTangentialLFIntegrator(f_left_re), new BoundaryTangentialLFIntegrator(f_left_im), rbmarker_left);
	b->AddBoundaryIntegrator(new BoundaryTangentialLFIntegrator(f_bottom_re), new BoundaryTangentialLFIntegrator(f_bottom_im), rbmarker_bottom);
	b->AddBoundaryIntegrator(new BoundaryTangentialLFIntegrator(f_top_re), new BoundaryTangentialLFIntegrator(f_top_im), rbmarker_top);
 	b->Assemble();
 	// 11. Define the solution vector x as a complex finite element grid function
    	//     corresponding to fespace.
    ComplexGridFunction x(fespace);
    x = 0.0;
	// Integrators inside the computational domain
	ConstantCoefficient one(1.0);
	ConstantCoefficient nksqr(-k * k);
	ConstantCoefficient ksi(0.0);
	ConstantCoefficient zero(0.0);
	ConstantCoefficient ik(-k); 
	SesquilinearForm *a = new SesquilinearForm(fespace, conv);
	a->AddDomainIntegrator(new CurlCurlIntegrator(one), new CurlCurlIntegrator(zero));
	a->AddDomainIntegrator(new VectorFEMassIntegrator(nksqr), new VectorFEMassIntegrator(ksi));
 	a->AddBoundaryIntegrator(new BoundaryMassIntegrator(zero), new BoundaryMassIntegrator(ik), rbmarker);
 	a->Assemble(0);
 	OperatorPtr A;
 	Vector B, X;
 	a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);
   	DDMSolver::ASMSolver_Impedence ASM_Imp(4,fespace->GetNDofs());
   	double eta=-k;
//    ASM_Imp.Setup(*mesh,*fespace,0,a,k,eta);
#ifdef MFEM_USE_SUITESPARSE
      	ComplexUMFPackSolver csolver(*A.As<ComplexSparseMatrix>());
	csolver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
        csolver.SetPrintLevel(1);
        csolver.Mult(B, X);
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
	ConstantCoefficient ksqr(k * k);
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
	// 12. Recover the solution as a finite element grid function.
	a->RecoverFEMSolution(X, *b, x);
	{
		ofstream mesh_ofs("cavity.mesh");
		mesh_ofs.precision(8);
		mesh->Print(mesh_ofs);

		ofstream sol_r_ofs("cavity_x_r_k10_xi10_gf.");
		sol_r_ofs.precision(8);
		x.real().Save(sol_r_ofs);
	}

	// 14. Send the solution by socket to a GLVis server.
	char vishost[] = "localhost";
	int  visport   = 19916;
	socketstream sol_sock(vishost, visport);
	sol_sock.precision(8);
	sol_sock << "solution\n" << *mesh << x.real() << "keys amrRljcUUuuu\n" << flush;

	if (fec)
   	{
	delete fespace;
	delete fec;
	}
   	delete mesh;

	return 0;
}

void rb_rhs_bottom_re(const Vector &x, Vector &v)
{
	vector<complex<double>> vval(v.Size());
	rb_rhs_bottom(x, vval);
	for (int i=0; i < v.Size(); ++i)
	{
		v[i] = vval[i].real();
	}
}
void rb_rhs_bottom_im(const Vector &x, Vector &v)
{
        vector<complex<double>> vval(v.Size());
        rb_rhs_bottom(x, vval);
        for (int i=0; i < v.Size(); ++i)
        {
                v[i] = vval[i].imag();
        }

}
void rb_rhs_top_re(const Vector &x, Vector &v)
{   
        vector<complex<double>> vval(v.Size());
        rb_rhs_top(x, vval);
        for (int i=0; i < v.Size(); ++i)
        {
                v[i] = vval[i].real();
        }
}
void rb_rhs_top_im(const Vector &x, Vector &v)
{   
        vector<complex<double>> vval(v.Size());
        rb_rhs_top(x, vval);
        for (int i=0; i < v.Size(); ++i)
        {
                v[i] = vval[i].imag();
        }
    
}   
void rb_rhs_left_re(const  Vector &x, Vector &v)
{
	v[0] = 0.0;
	v[1] = 0.0;
}
void rb_rhs_left_im(const Vector &x, Vector &v)
{
	v[0] = 0.0;
	v[1] = -2.0 * k;


}
void rb_rhs_bottom(const Vector &x, vector<complex<double>> &fv)
{  
	complex<double> zi = complex<double>(0., 1.);
	double x0 = x(0);
	fv[0] = zi * k * exp(zi * k * x0);
	fv[1] = 0.0;	
}
void rb_rhs_top(const Vector &x, vector<complex<double>> &fv)
{
        complex<double> zi = complex<double>(0., 1.);
        double x0 = x(0);
        fv[0] = -zi * k * exp(zi * k * x0);
        fv[1] = 0.0;    

}