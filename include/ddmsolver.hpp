#ifndef MFEM_DDM_SOLVER_HPP
#define MFEM_DDM_SOLVER_HPP
#include "mfem.hpp"
#include <iostream>
namespace mfem
{
namespace DDMSolver
{
class ASMSolver : public Solver
{
public:
	ASMSolver();

	ASMSolver(int n_dom, int nrow_,bool oras=false);

	virtual ~ASMSolver();

	virtual void Mult(const Vector& b, Vector& x) const override;

	virtual void SetOperator(const Operator& op) override;

	virtual void setUnityVector(int ndof);

	virtual void setCoarseSpace(Operator *P,SparseMatrix &A_c);

	int * SplitMesh(const Mesh & mesh, int ncommon, bool if_print=false);
	Array<Array<int> *> GetOverlapindex(Array<Array<int> * > &elem_list,Mesh & mesh,int delta);
	SparseMatrix RART(SparseMatrix & A,Array<int> & ind);
	SparseMatrix* subAssemble(int i,FiniteElementSpace &fespace,BilinearForm* a);
	void subBdrAssemble(int i,Mesh &mesh,FiniteElementSpace &fespace,BilinearForm *a,SparseMatrix &subA);
	virtual void setsubsolver();
	void setsubdofmap();
	void Setup(Mesh &mesh,SparseMatrix A,FiniteElementSpace Space,int delta);
	void Setup(Mesh &mesh,FiniteElementSpace Space,int delta,BilinearForm* a);
	Array<Array<int> * > getsubelemlist();
	Array<SparseMatrix *> getsubMatrixlist();
	
protected:
	int num_domains;
	int nrow;
	bool UseUnity;
	bool Re_Assemble;
	bool UseTwolevel;
	Operator *Prolongation;
	UMFPackSolver *solver_coarse;
	Vector *UnityVector;
	Array<Array<int> * > sub_dof_list;//dof belong to subdomain i is store in sub_dof_list[i]
	Array<SparseMatrix *> sub_A_list;
	Array<Array<int> * > sub_elem_list;
	Array<UMFPackSolver *> sub_solver_list;
	Array<Array<int> *> sub_bd_list;
	Array<std::map<int,int> *> sub_dof_map;
};


class ASMSolver_Impedence : public ASMSolver
{
public:

	ASMSolver_Impedence();

	ASMSolver_Impedence(int n_dom, int nrow_);

	virtual ~ASMSolver_Impedence();

	virtual void Mult(const Vector& b, Vector& x) const override;

	void setsubsolver() override;

	void setCoarseSpace(Operator *P,SparseMatrix &A_c) override;

	void setCoarseSpace(Operator *P,ComplexSparseMatrix &A_c);
	
	virtual void setUnityVector(int ndof) override;

	void Setup(Mesh &mesh,FiniteElementSpace &fespace,int delta,SesquilinearForm* a,double &k,double &eta);

	ComplexSparseMatrix* subAssemble(int i,FiniteElementSpace &fespace,SesquilinearForm* a,double &eta);

	virtual void subBdrAssemble(int i,Mesh &mesh,FiniteElementSpace &fespace,SesquilinearForm* a,SparseMatrix &subA_r,SparseMatrix &subA_i,double &eta);

protected: 
	ComplexUMFPackSolver *solver_coarse;
	Array<ComplexSparseMatrix *> sub_A_list;
	Array<ComplexUMFPackSolver *> sub_solver_list;
};

class ASMSolver_ComplexDir : public ASMSolver
{
	public:

	ASMSolver_ComplexDir();

	ASMSolver_ComplexDir(int n_dom, int nrow_);

	virtual ~ASMSolver_ComplexDir();

	void setsubsolver() override;

	void setCoarseSpace(Operator *P,ComplexSparseMatrix &A_c);

	void setCoarseSpace(Operator *P,SparseMatrix &A_c) override;

	virtual void Mult(const Vector& b, Vector& x) const override;

	protected:
	ComplexUMFPackSolver *solver_coarse;
	Array<ComplexSparseMatrix * > sub_A_complex;
	Array<ComplexUMFPackSolver *> sub_solver_list;
};


class twolevelsolver : public Solver
{
public:
	twolevelsolver();
	twolevelsolver(int nrow,SparseMatrix O_c,SparseMatrix O_f,FiniteElementSpaceHierarchy *FEMH,Solver* csolver,Solver* fsolver);
	virtual ~twolevelsolver();
	virtual void SetOperator(const Operator& op) override;
	void solve(Vector& b, Vector& x,Vector& p);
	virtual void Mult(const Vector& b, Vector& x) const override;
protected:
	SparseMatrix A_c;
	SparseMatrix A_f;
	Solver* coarsesolver;
	Solver* finesolver;
	FiniteElementSpaceHierarchy* FEMHierarchy;
};
}
}
 #endif
