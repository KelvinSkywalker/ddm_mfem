#include <iostream>
#include <cassert>
#include "ddmsolver.hpp"
#include "metis.h"

#include <fstream>
#include <map>
#include <iostream>

namespace mfem
{
namespace DDMSolver
{
void csrtocoo(int *csr_row,int *csr_col,double *csr_val,int *coo_row,int *coo_col,double *coo_val,int height)
{
	for (int i=0,k=0;i<height;i++)
	{
		for(int j=csr_row[i];j<csr_row[i+1];j++)
		{
			coo_row[k]=i;
			coo_col[k]=(csr_col[j]);
			coo_val[k]=(csr_val[j]);
			k++;
		}
	}

}
void cootocsr(int *csr_row,int *csr_col,double *csr_val,int *coo_row,int* coo_col,double *coo_val,int height,int length)
{
	Array<int> perrow(height);
	for(int i=0;i<height;i++)
	{
		perrow[i]=0;
	}
	for(int i=0;i<length;i++)
	{
		perrow[coo_row[i]]+=1;
	}
	Array<int> copyperrow=perrow;
	csr_row[0]=0;
	for(int i=1;i<=height;i++)
	{
		csr_row[i]=csr_row[i-1]+perrow[i-1];
	}
	for(int i=0;i<length;i++)
	{
		int tempi=coo_row[i];
		int tempoffset=csr_row[tempi]+copyperrow[tempi]-perrow[tempi];
		csr_col[tempoffset]=coo_col[i];
		csr_val[tempoffset]=coo_val[i];
		perrow[tempi]-=1;
	}
}

//ASMSolver
ASMSolver::ASMSolver():num_domains(0), nrow(0) {}

ASMSolver::ASMSolver(int n_dom, int nrow_,bool useunity): Solver(nrow_, nrow_),
											num_domains(n_dom),
											nrow(nrow_), 
											sub_dof_list(n_dom), 
											sub_A_list(n_dom),
											sub_elem_list(n_dom),
											sub_solver_list(n_dom),
											sub_bd_list(n_dom),
											sub_dof_map(n_dom) 
{
	UseUnity=useunity;
	Re_Assemble=false;
	UseTwolevel=false;
}

ASMSolver::~ASMSolver() {}

void ASMSolver::SetOperator(const Operator& op)
{
   MFEM_ABORT("SetOperator is not supported in DDM!");
}


// the same as the function Solve() in ddm.py
void ASMSolver::Mult(const Vector& b, Vector& x) const
{
	x=double(0);
	for(int i=0;i<num_domains;i++)
	{
		SparseMatrix subA=*sub_A_list[i];
		Array<int> subdof=*sub_dof_list[i];
		Vector subx(sub_dof_list[i]->Size()),subb;
		b.GetSubVector(subdof,subb);
		// UMFPackSolver umf_solver;
    	// umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
    	// umf_solver.SetOperator(subA);
    	// umf_solver.Mult(subb, subx);
		sub_solver_list[i]->Mult(subb,subx);
		if(UseUnity)
		{
			Vector Di;
			UnityVector->GetSubVector(subdof,Di);
			subx*=Di;
		}
		if(Re_Assemble)
		{
			subx.SetSubVector(*sub_bd_list[i],0.0);
		}
		x.AddElementVector(subdof,subx);
	}
	if(UseTwolevel)
	{
		Vector b_c(Prolongation->Width()),x_c(Prolongation->Width());
		x_c=0.0;
		Prolongation->MultTranspose(b,b_c);
		solver_coarse->Mult(b_c,x_c);
		Prolongation->AddMult(x_c,x);
	}

}
void ASMSolver::setsubsolver()
{
	for(int i=0;i<num_domains;i++)
	{
		sub_solver_list[i]=new UMFPackSolver(1);
		sub_solver_list[i]->SetOperator(*sub_A_list[i]);
		sub_solver_list[i]->Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
	}
}
void ASMSolver::setUnityVector(int ndof)
{
	UseUnity = true;
	Vector *v=new Vector(ndof);
	Vector dofcount(ndof);
	dofcount=0.0;
	*v=1.0;
	for(int i=0;i<num_domains;i++)
	{
		for(int j=0;j<sub_dof_list[i]->Size();j++)
		{
			dofcount[(*sub_dof_list[i])[j]]+=1;
		}
	}
	*v/=dofcount;
	UnityVector=v;
	return;
}
void ASMSolver::setCoarseSpace(Operator *P,SparseMatrix &A_c)
{
	UseTwolevel=true;
	Prolongation=P;
	solver_coarse=new UMFPackSolver(true);
	solver_coarse->SetOperator(A_c);
	solver_coarse->Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
}
void ASMSolver::setsubdofmap()
{
	for(int i=0;i<num_domains;i++)
	{
		std::map<int,int> *subdofmap_i=new std::map<int,int>;
		Array<int> &tempdofarray=*sub_dof_list[i];
		for(int j=0;j<sub_dof_list[i]->Size();j++)
		{
			(*subdofmap_i)[tempdofarray[j]]=j;
		}
		sub_dof_map[i]=subdofmap_i;
	}
	return;
};
int * ASMSolver::SplitMesh(const Mesh & mesh, int ncommon, bool if_print)
{
	int num_elem = mesh.GetNE();
	int num_node = mesh.GetNV();

	// The size of the eptr array is n + 1, where n is the number of elements in the mesh
	int * eptr = new int[num_elem + 1];
	eptr[0] = 0;

	// The size of the eind array is of size equal to the sum of the number of nodes in all the elements of the mesh, which is not the number of nodes in the mesh!
	int len_eind = num_elem * mesh.GetElement(0)->GetNVertices();
	int * eind = new int[len_eind];

	// the eptr and eind are the csr_row_vec and csr_col_vec 
	for(int i=0;i<num_elem;i++)
	{
		const Element * elem = mesh.GetElement(i);
		eptr[i+1] = eptr[i] + elem->GetNVertices();

		const int * indices = elem->GetVertices();
		int offset = 0;
		for(int j=eptr[i];j<eptr[i+1];j++)
		{
			eind[j] = indices[offset];
			offset += 1;
		}
	}

	// use metis default options
	idx_t options[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options);

	// Upon successful completion, this variable stores either the edgecut or the total communication volume of the dual graph’s partitioning.
	int objval = 0;

	// This is a vector of size num_elem that upon successful completion stores the partition vector for the elements of the mesh.
	// for example [0,1,0,1]
	int * epart = new int[num_elem];

	// This is a vector of size num_node that upon successful completion stores the partition vector for the nodes of the mesh.
	int * npart = new int[num_node];

	// call metis
	int result = METIS_PartMeshDual(&num_elem, &num_node,
          eptr, eind, nullptr, nullptr, &ncommon, &num_domains, nullptr, options,
          &objval, epart, npart);
	assert(result == METIS_OK);

	// print the contents of epart
	if (if_print)
	{
		std::cout << "epart = ";
		for(int i=0;i<num_elem;i++)
			std::cout<< epart[i] << " ";
		std::cout << std::endl;
	}

	delete[] eptr;
	delete[] eind;

	// the npart is not needed, then it can be deleted 
	delete[] npart;

	// note that the epart is used to create sub_elem_idx in GetOverlapIndex() function, which should be deleted in the GetOverlapIndex() function! TODO
	
	return epart;
}
Array<Array<int> *> ASMSolver::GetOverlapindex(Array<Array<int> * > &elem_list,Mesh & mesh,int delta)
{
	Array<Array<int> * > bdlist(num_domains);
	Array<Array<int> * > subdomainvertex(num_domains);
	Table *ver2el=mesh.GetVertexToElementTable();
	for(int n=0;n<num_domains;n++)
	{
		bdlist[n]=new Array<int>;
		elem_list[n]->Sort();
		*bdlist[n]=*elem_list[n];
	}
	if(delta==0)
	{
		return elem_list;
	}
	for(int m=0;m<delta;m++)
	{
		for(int n=0;n<num_domains;n++)
		{
			subdomainvertex[n]=new Array<int>; //find the vertex in bdlist
			for(int i=0;i<bdlist[n]->Size();i++)
			{
				Array<int> tempvertex;
				int tempelem=(*bdlist[n])[i];
				mesh.GetElementVertices(tempelem,tempvertex);
				for(int j=0;j<tempvertex.Size();j++)
				{
					subdomainvertex[n]->Union(tempvertex[j]);
				}
			}
			delete bdlist[n];
			bdlist[n]=new Array<int>;
			for(int i=0;i<subdomainvertex[n]->Size();i++)
			{
				int tv=(*subdomainvertex[n])[i];
				Array<int> te;
				ver2el->GetRow(tv,te);
				for(int j=0;j<te.Size();j++)
				{
					int tte=te[j];
					if(elem_list[n]->Find(tte)==-1)
					{
						bdlist[n]->Union(tte);
						elem_list[n]->Append(tte);
					}
				}
			}
			delete subdomainvertex[n];
		}
	}
	subdomainvertex.DeleteAll();
	bdlist.DeleteAll();
	return elem_list;
}
SparseMatrix ASMSolver::RART(SparseMatrix & A,Array<int> & ind)
{
	int *I=A.GetI(),*J=A.GetJ();
	ind.Sort();
	std::map<int,int> dict;
	double *Data=A.GetData();
	int h=ind.Size();
	int l1=0;
	for(int i=0;i<h;i++)
	{
		int temr=ind[i];
		l1+=(I[temr+1]-I[temr]);
		dict[temr]=i;
	}
	int *coo_row1=new int[l1],*coo_col1 =new int[l1];
	double *coo_val1=new double[l1];
	for(int i=0,now=0;i<h;i++)
	{
		int temr=ind[i];
		for(int j=I[temr];j<I[temr+1];j++)
		{
			coo_row1[now]=dict[temr];
			coo_col1[now]=J[j];
			coo_val1[now]=Data[j];
			now++;
		}
	}
	int *csri1=new int[A.Height()+1],*csrj1=new int[l1];
	double *csrv1=new double[l1];
	cootocsr(csri1,csrj1,csrv1,coo_col1,coo_row1,coo_val1,A.Height(),l1);
	int l2=0;
	for(int i=0;i<h;i++)
	{
		int temr=ind[i];
		l2+=(csri1[temr+1]-csri1[temr]);
	}
	int *coo_row2=new int[l2],*coo_col2=new int[l2];
	double *coo_val2=new double[l2];	
	for(int i=0,now=0;i<h;i++)
	{
		int temr=ind[i];
		for(int j=csri1[temr];j<csri1[temr+1];j++)
		{
			coo_row2[now]=dict[temr];
			coo_col2[now]=csrj1[j];
			coo_val2[now]=csrv1[j];
			now++;
		}
	}
	int *csri2=new int[h+1],*csrj2=new int[l2];
	double *csrv2=new double[l2];
	cootocsr(csri2,csrj2,csrv2,coo_col2,coo_row2,coo_val2,h,l2);
	SparseMatrix subA(csri2,csrj2,csrv2,h,h,0,0,0);
	delete [] coo_col1;
	delete [] coo_row1;
	delete [] coo_val1;
	delete [] coo_row2;
	delete [] coo_col2;
	delete [] coo_val2;
	delete [] csri1;
	delete [] csrj1;
	delete [] csrv1;
	return subA;
}
void ASMSolver::Setup(Mesh &mesh,SparseMatrix A,FiniteElementSpace Space,int delta)
{
	int *elemid=SplitMesh(mesh,num_domains);
	int NE=mesh.GetNE();
	for(int i=0;i<num_domains;i++)
	{
		sub_elem_list[i]=new Array<int>;
		sub_dof_list[i]=new Array<int>;
		sub_A_list[i]=new SparseMatrix;
	}
	for(int i=0;i<NE;i++)
	{
		int tdomain=elemid[i];
		sub_elem_list[tdomain]->Append(i);
	}
	delete [] elemid;
	sub_elem_list=GetOverlapindex(sub_elem_list,mesh,delta);
	for(int i=0;i<num_domains;i++)
	{
		// std::cout<<"sub elem list "<<i<<" is "<<std::endl;
		// sub_elem_list[i]->Print();
		for(int j=0;j<sub_elem_list[i]->Size();j++)
		{
			Array<int> tdof;
			int tempe=(*sub_elem_list[i])[j];
			Space.GetElementDofs(tempe,tdof);
			Space.AdjustVDofs(tdof);
			for(int j=0;j<tdof.Size();j++)
			{
				sub_dof_list[i]->Union(tdof[j]);
			}
		}
		*sub_A_list[i]=RART(A,*sub_dof_list[i]);
		sub_A_list[i]->Finalize();
	}
	setsubsolver();
}
void ASMSolver::Setup(Mesh &mesh,FiniteElementSpace Space,int delta,BilinearForm* a)
{
	Re_Assemble=true;
	int *elemid=SplitMesh(mesh,num_domains);
	int NE=mesh.GetNE();
	for(int i=0;i<num_domains;i++)
	{
		sub_elem_list[i]=new Array<int>;
		sub_dof_list[i]=new Array<int>;
	}
	for(int i=0;i<NE;i++)
	{
		int tdomain=elemid[i];
		sub_elem_list[tdomain]->Append(i);
		mesh.SetAttribute(i,num_domains+10);
	}
	mesh.SetAttributes();
	delete [] elemid;
	sub_elem_list=GetOverlapindex(sub_elem_list,mesh,delta);
	for(int i=0;i<num_domains;i++)
	{
		// std::cout<<"sub elem list "<<i<<" is "<<std::endl;
		// sub_elem_list[i]->Print();
		for(int j=0;j<sub_elem_list[i]->Size();j++)
		{
			Array<int> tdof;
			int tempe=(*sub_elem_list[i])[j];
			Space.GetElementDofs(tempe,tdof);
			Space.AdjustVDofs(tdof);
			for(int j=0;j<tdof.Size();j++)
			{
				sub_dof_list[i]->Union(tdof[j]);
			}
		}
		sub_dof_list[i]->Sort();
	}
	setsubdofmap();
	for(int i=0;i<num_domains;i++)
	{
		sub_dof_list[i]->Sort();
		// std::cout<<"sub dof list "<<i<<" is "<<std::endl;
		// sub_dof_list[i]->Print();
		sub_A_list[i]=subAssemble(i,Space,a);
		sub_A_list[i]->Finalize();
		// SparseMatrix rarta=RART(A,*sub_dof_list[i]);
		// DenseMatrix a1,a2;
		// rarta.ToDenseMatrix(a1);
		// a1.Neg();
		// sub_A_list[i]->ToDenseMatrix(a2);
		// a1+=a2;
		// std::cout<<"the difference is"<<std::endl;
		// double *norm=new double;
		// a1.Norm2(norm);
		// std::cout<<a1.MaxMaxNorm()<<std::endl;
	}
	setsubsolver();
}
SparseMatrix* ASMSolver::subAssemble(int i,FiniteElementSpace &fespace,BilinearForm* a)
{
	SparseMatrix *subA=new SparseMatrix(sub_dof_list[i]->Size(),sub_dof_list[i]->Size());
	*subA=0;
	for(int j=0;j<sub_elem_list[i]->Size();j++)
	{
		int telem=(*sub_elem_list[i])[j];
		DenseMatrix elemmatrix;
		a->ComputeElementMatrix(telem,elemmatrix);
		Array<int> glodof,subdof;
		fespace.GetElementDofs(telem,glodof);
		subdof=glodof;
		for(int k=0;k<glodof.Size();k++)
		{
			if(glodof[k]<0)
			{
				int adjustglodof,adjustsubdof;
				adjustglodof=-glodof[k]-1;
				adjustsubdof=(*sub_dof_map[i])[adjustglodof];
				subdof[k]=-adjustsubdof-1;
			}
			else
			{
				subdof[k]=(*sub_dof_map[i])[glodof[k]];
			}
		}
		subA->AddSubMatrix(subdof,subdof,elemmatrix);
		// std::cout<<"subdof is "<<std::endl;
		// subdof.Print();
	}
	subA->Finalize();
	subBdrAssemble(i,*fespace.GetMesh(),fespace,a,*subA);
	return subA;
}
void ASMSolver::subBdrAssemble(int i,Mesh &mesh,FiniteElementSpace &fespace,BilinearForm *a,SparseMatrix &subA)
{
	for(int j=0;j<sub_elem_list[i]->Size();j++)
	{
		mesh.SetAttribute((*sub_elem_list[i])[j],i+1);
	}
	mesh.SetAttributes();
	mesh.Finalize();
	Array<int> t_i(1);
	t_i=i+1;
	SubMesh submesh_i(SubMesh::CreateFromDomain(mesh,t_i));
	// submesh_i.GenerateBoundaryElements();
	Array<int> subelem2face=submesh_i.GetParentFaceIDMap();
	Array<int> *eliminateBdrDof=new Array<int>;
	Array<int> *insideBdrDof=new Array<int>(*eliminateBdrDof);
	// Array<int> *globalBdrDof=new Array<int>;
	for(int j=0;j<submesh_i.GetNBE();j++)
	{
		int *subfaceid = new int, *subfaceorietation= new int;
		submesh_i.GetBdrElementFace(j,subfaceid,subfaceorietation);
		int subbdrattr=submesh_i.GetBdrAttribute(j);
		bool isglobalbdrdof=0;
		if(mesh.bdr_attributes.Find(subbdrattr)!=(-1))
		{
			isglobalbdrdof=true;
		}
		//subfaceid=submesh_i.GetBdrElementEdgeIndex(j);
		int parentfaceid=subelem2face[*subfaceid];
		Array<int> temdof;
		fespace.GetFaceDofs(parentfaceid,temdof);
		fespace.AdjustVDofs(temdof);
		for(int k=0;k<temdof.Size();k++)
		{
			int subdofid=(*sub_dof_map[i])[temdof[k]];
			eliminateBdrDof->Union(subdofid);
			if(!isglobalbdrdof)
			{
				insideBdrDof->Union(subdofid);
			}
		}
		delete subfaceid;
		delete subfaceorietation;
	}
	// for(int k=0;k<globalBdrDof->Size();k++)
	// {
	// 	int el=(*globalBdrDof)[k];
	// 	insideBdrDof->DeleteFirst(el);
	// }
	subA.Finalize();
	for(int k=0;k<eliminateBdrDof->Size();k++)
	{
		int tempdof=(*eliminateBdrDof)[k];
		// int temval=subA.Elem(tempdof,tempdof);
		subA.EliminateRowCol(tempdof,mfem::Operator::DiagonalPolicy::DIAG_ONE);
		// subA.Set(tempdof,tempdof,temval);
	}
	insideBdrDof->Sort();
	// std::cout<<"insidebdrdof of "<<i<<" are "<<insideBdrDof->Size()<<std::endl;
	delete eliminateBdrDof;
	sub_bd_list[i]=insideBdrDof;
	// std::cout<<"bdr dof is "<<std::endl;
	// eliminateBdrDof.Print();
}
Array<Array<int> * > ASMSolver::getsubelemlist()
{
	return sub_elem_list;
}
Array<SparseMatrix *> ASMSolver::getsubMatrixlist()
{
	return sub_A_list;
}

//ASMSolver for Impedence BC

ASMSolver_Impedence::ASMSolver_Impedence():ASMSolver(0,0) {};

ASMSolver_Impedence::ASMSolver_Impedence(int n_dom, int nrow_): 
	ASMSolver(n_dom, nrow_),
	sub_A_list(n_dom),
	sub_solver_list(n_dom){};

ASMSolver_Impedence::~ASMSolver_Impedence(){};

void ASMSolver_Impedence::Mult(const Vector& b, Vector& x) const
{
	x=0.0;
	for(int i=0;i<num_domains;i++)
	{
		ComplexSparseMatrix subA=*sub_A_list[i];
		Array<int> subdof=*sub_dof_list[i];
		Vector subx(2*sub_dof_list[i]->Size()),subb(2*sub_dof_list[i]->Size());
		Array<int> subdof2=subdof;
		subdof2.SetSize(2*subdof.Size());
		int subs=subdof.Size();
		for(int j=0;j<subs;j++)
		{
			subdof2[subs+j]=subdof[j]+x.Size()/2;
		}
		b.GetSubVector(subdof2,subb);
		// std::cout<<"subsolver i has height "<<sub_solver_list[i]->Height()<<std::endl;
		// std::cout<<"subvector i has height "<<subb.Size()<<std::endl;
		sub_solver_list[i]->Mult(subb,subx);
		if(UseUnity)
		{
			Vector Di;
			UnityVector->GetSubVector(subdof2,Di);
			subx*=Di;
		}
		x.AddElementVector(subdof2,subx);
	}
	if(UseTwolevel)
	{
		int ndof_f=Prolongation->Height(),ndof_c=Prolongation->Width();
		Vector b_c(ndof_c*2),x_c(ndof_c*2);
		Vector b0(b);
		Vector b_r(b0,0,ndof_f),b_i(b0,ndof_f,ndof_f);
		Vector bc_r(b_c,0,ndof_c),bc_i(b_c,ndof_c,ndof_c);
		Prolongation->MultTranspose(b_r,bc_r);
		Prolongation->MultTranspose(b_i,bc_i);
		solver_coarse->Mult(b_c,x_c);
		Vector xc_r(x_c,0,ndof_c),xc_i(x_c,ndof_c,ndof_c);
		Vector x_r(x,0,ndof_f),x_i(x,ndof_f,ndof_f);
		Prolongation->AddMult(xc_r,x_r);
		Prolongation->AddMult(xc_i,x_i);
	}
}

void ASMSolver_Impedence::subBdrAssemble(int i,Mesh &mesh,FiniteElementSpace &fespace,SesquilinearForm* a,SparseMatrix &subA_r,SparseMatrix &subA_i,double &eta)
{
	ComplexOperator::Convention conv=a->GetConvention();
	for(int j=0;j<sub_elem_list[i]->Size();j++)
	{
		mesh.SetAttribute((*sub_elem_list[i])[j],i+1);
	}
	mesh.Finalize();
	Array<int> t_i(1);
	t_i=i+1;
	SubMesh submesh_i=SubMesh::CreateFromDomain(mesh,t_i);
	// submesh_i.Finalize();
	// submesh_i.GenerateBoundaryElements();
	int order_i=fespace.GetOrder(0);
	int dim_i=mesh.Dimension();
	FiniteElementCollection *fec_i= new ND_FECollection(order_i, dim_i);;
   	FiniteElementSpace *fespace_i= new FiniteElementSpace(&submesh_i,fec_i);;
	SesquilinearForm *a_i = new SesquilinearForm(fespace_i, conv);
	ConstantCoefficient zero(0.0);
	ConstantCoefficient ieta(-eta);
	a_i->AddBoundaryIntegrator(new VectorFEMassIntegrator(zero), new VectorFEMassIntegrator(ieta));
	Array<int> subelem2face=submesh_i.GetParentFaceIDMap();
	for(int j=0;j<submesh_i.GetNBE();j++)
	{
		int *subfaceid=new int,*subfaceorientation=new int;
		submesh_i.GetBdrElementFace(j,subfaceid,subfaceorientation);
		// int istruebdr;
		// istruebdr=mesh.bdr_attributes.Find(submesh_i.GetBdrAttribute(j));
		int parentfaceid=subelem2face[*subfaceid];
		Array<int> subbdrelemdof,glodof,subdof;
		fespace.GetFaceDofs(parentfaceid,glodof);
		subdof=glodof;
		DofTransformation *doftrans;
		doftrans=fespace_i->GetBdrElementVDofs(j,subbdrelemdof);
		for(int k=0;k<glodof.Size();k++)
		{
			if(glodof[k]<0)
			{
				int adjustglodof,adjustsubdof;
				adjustglodof=-glodof[k]-1;
				adjustsubdof=(*sub_dof_map[i])[adjustglodof];
				subdof[k]=-adjustsubdof-1;
			}
			else
			{
				int adjustsubdof=(*sub_dof_map[i])[glodof[k]];
				subdof[k]=adjustsubdof;
			}
		}
		DenseMatrix subBdrmatrix_r,subBdrmatrix_i;
		BilinearForm &a_image=a_i->imag();
		BilinearForm &a_real=a_i->real();
		a_image.ComputeBdrElementMatrix(j,subBdrmatrix_i);
		a_real.ComputeBdrElementMatrix(j,subBdrmatrix_r);
		if(doftrans)
		{
			doftrans->TransformDual(subBdrmatrix_i);
			doftrans->TransformDual(subBdrmatrix_r);
		}
		subA_i.AddSubMatrix(subdof,subdof,subBdrmatrix_i,0);
		subA_r.AddSubMatrix(subdof,subdof,subBdrmatrix_r,0);
	}
	delete a_i;	
	return;
}

void ASMSolver_Impedence::setUnityVector(int ndof)
{
	UseUnity = true;
	Vector *v=new Vector(ndof);
	Vector dofcount(ndof);
	dofcount=0.0;
	*v=1.0;
	for(int i=0;i<num_domains;i++)
	{
		for(int j=0;j<sub_dof_list[i]->Size();j++)
		{
			dofcount[(*sub_dof_list[i])[j]]+=1;
		}
	}
	*v/=dofcount;
	Vector *v2=new Vector(2*ndof);
	*v2=0.0;
	v2->AddSubVector(*v,0);
	v2->AddSubVector(*v,ndof);
	UnityVector=v2;
	delete v;
	return;
}

void ASMSolver_Impedence::setCoarseSpace(Operator *P,SparseMatrix &A_c)
{
	std::cout<<"the coarse matrix should be complex Sparse Matirx"<<std::endl;
}

void ASMSolver_Impedence::setCoarseSpace(Operator *P,ComplexSparseMatrix &A_c)
{
	UseTwolevel=true;
	Prolongation=P;
	solver_coarse=new ComplexUMFPackSolver(true);
	solver_coarse->SetOperator(A_c);
	solver_coarse->Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
};

ComplexSparseMatrix* ASMSolver_Impedence::subAssemble(int i,FiniteElementSpace &fespace,SesquilinearForm* a,double &eta)
{
	int h=sub_dof_list[i]->Size();
	int skipzeros=0;
	ComplexOperator::Convention conv=a->GetConvention();
	SparseMatrix *subA_r=new SparseMatrix(h,h),*subA_i=new SparseMatrix(h,h);
	*subA_r=0;
	*subA_i=0;
	BilinearForm &a_r=a->real(),&a_i=a->imag();
	for(int j=0;j<sub_elem_list[i]->Size();j++)
	{
		int telem=(*sub_elem_list[i])[j];
		DenseMatrix elemmatrix_r,elemmatrix_i;
		a_r.ComputeElementMatrix(telem,elemmatrix_r);
		a_i.ComputeElementMatrix(telem,elemmatrix_i);
		Array<int> glodof,subdof;
		fespace.GetElementDofs(telem,glodof);
		subdof=glodof;
		for(int k=0;k<glodof.Size();k++)
		{
			if(glodof[k]<0)
			{
				int adjustglodof,adjustsubdof;
				adjustglodof=-glodof[k]-1;
				adjustsubdof=(*sub_dof_map[i])[adjustglodof];
				subdof[k]=-adjustsubdof-1;
			}
			else
			{
				subdof[k]=(*sub_dof_map[i])[glodof[k]];
			}
		}
		subA_r->AddSubMatrix(subdof,subdof,elemmatrix_r,skipzeros);
		subA_i->AddSubMatrix(subdof,subdof,elemmatrix_i,skipzeros);
		// std::cout<<"subdof is "<<std::endl;
		// subdof.Print();
	}
	subBdrAssemble(i,*fespace.GetMesh(),fespace,a,*subA_r,*subA_i,eta);
	//std::cout<<"no zero elem  are "<<subA_i->NumNonZeroElems()<<" + "<<subA_r->NumNonZeroElems()<<std::endl;
	subA_r->Finalize(0);
	subA_i->Finalize(0);
	ComplexSparseMatrix *subA=new ComplexSparseMatrix(subA_r,subA_i,0,0,conv);
	return subA;
};

void ASMSolver_Impedence::setsubsolver()
{
		for(int i=0;i<num_domains;i++)
	{
		sub_solver_list[i]=new ComplexUMFPackSolver(1);
		sub_solver_list[i]->Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
		sub_solver_list[i]->SetOperator(*sub_A_list[i]);
	}
}

void ASMSolver_Impedence::Setup(Mesh &mesh,FiniteElementSpace &fespace,int delta,SesquilinearForm* a,double &k,double &eta)
{
	int *elemid=SplitMesh(mesh,num_domains);
	int NE=mesh.GetNE();
	for(int i=0;i<num_domains;i++)
	{
		sub_elem_list[i]=new Array<int>;
		sub_dof_list[i]=new Array<int>;
	}
	for(int i=0;i<NE;i++)
	{
		int tdomain=elemid[i];
		sub_elem_list[tdomain]->Append(i);
		mesh.SetAttribute(i,num_domains+10);
	}
	delete [] elemid;
	sub_elem_list=GetOverlapindex(sub_elem_list,mesh,delta);
	for(int i=0;i<num_domains;i++)
	{
		// std::cout<<"sub elem list "<<i<<" is "<<std::endl;
		// sub_elem_list[i]->Print();
		for(int j=0;j<sub_elem_list[i]->Size();j++)
		{
			Array<int> tdof;
			int tempe=(*sub_elem_list[i])[j];
			fespace.GetElementDofs(tempe,tdof);
			fespace.AdjustVDofs(tdof);
			for(int j=0;j<tdof.Size();j++)
			{
				sub_dof_list[i]->Union(tdof[j]);
			}
		}
		sub_dof_list[i]->Sort();
		// std::cout<<"sub dof list "<<i<<" is "<<std::endl;
		// sub_dof_list[i]->Print();
	}
	setsubdofmap();
	for(int i=0;i<num_domains;i++)
	{
		sub_A_list[i]=subAssemble(i,fespace,a,eta);
	}
	setsubsolver();
}

//Complex Dirichlet Solver
ASMSolver_ComplexDir::ASMSolver_ComplexDir():ASMSolver(0,0) {};

ASMSolver_ComplexDir::ASMSolver_ComplexDir(int n_dom, int nrow_): 
	ASMSolver(n_dom, nrow_),
	sub_A_complex(n_dom),
	sub_solver_list(n_dom){};

ASMSolver_ComplexDir::~ASMSolver_ComplexDir(){};

void ASMSolver_ComplexDir::setsubsolver()
{
	for(int i=0;i<num_domains;i++)
	{
		SparseMatrix *zerosImagepart=new SparseMatrix(*sub_A_list[i]);
		double *imagedata=zerosImagepart->WriteData();
		for(int j=0;j<zerosImagepart->NumNonZeroElems();j++)
		{
			imagedata[j]=0;
		}
		sub_A_complex[i]= new ComplexSparseMatrix(sub_A_list[i],zerosImagepart,0,0);
		sub_solver_list[i]=new ComplexUMFPackSolver();
		sub_solver_list[i]->Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
		sub_solver_list[i]->SetOperator(*sub_A_complex[i]);
	}
}

void ASMSolver_ComplexDir::setCoarseSpace(Operator *P,SparseMatrix &A_c)
{
	std::cout<<"the coarse matrix should be complex Sparse Matirx"<<std::endl;
}

void ASMSolver_ComplexDir::setCoarseSpace(Operator *P,ComplexSparseMatrix &A_c)
{
	UseTwolevel=true;
	Prolongation=P;
	solver_coarse=new ComplexUMFPackSolver(true);
	solver_coarse->SetOperator(A_c);
	solver_coarse->Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
};

void ASMSolver_ComplexDir::Mult(const Vector& b, Vector& x) const
{
	x=0.0;
	for(int i=0;i<num_domains;i++)
	{
		ComplexSparseMatrix subA=*sub_A_complex[i];
		Array<int> subdof=*sub_dof_list[i];
		Vector subx(2*sub_dof_list[i]->Size()),subb(2*sub_dof_list[i]->Size());
		Array<int> subdof2=subdof;
		subdof2.SetSize(2*subdof.Size());
		int subs=subdof.Size();
		for(int j=0;j<subs;j++)
		{
			subdof2[subs+j]=subdof[j]+x.Size()/2;
		}
		b.GetSubVector(subdof2,subb);
		// std::cout<<"subsolver i has height "<<sub_solver_list[i]->Height()<<std::endl;
		// std::cout<<"subvector i has height "<<subb.Size()<<std::endl;
		sub_solver_list[i]->Mult(subb,subx);
		if(UseUnity)
		{
			Vector Di;
			UnityVector->GetSubVector(subdof2,Di);
			subx*=Di;
		}
		Vector x_r(x,0,x.Size()/2),x_i(x,x.Size()/2,x.Size()/2);
		x_r.SetSubVector(*sub_bd_list[i],0.0);
		x_i.SetSubVector(*sub_bd_list[i],0.0);
		x.AddElementVector(subdof2,subx);
	}
	if(UseTwolevel)
	{
		int ndof_f=Prolongation->Height(),ndof_c=Prolongation->Width();
		Vector b_c(ndof_c*2),x_c(ndof_c*2);
		Vector b0(b);
		Vector b_r(b0,0,ndof_f),b_i(b0,ndof_f,ndof_f);
		Vector bc_r(b_c,0,ndof_c),bc_i(b_c,ndof_c,ndof_c);
		Prolongation->MultTranspose(b_r,bc_r);
		Prolongation->MultTranspose(b_i,bc_i);
		solver_coarse->Mult(b_c,x_c);
		Vector xc_r(x_c,0,ndof_c),xc_i(x_c,ndof_c,ndof_c);
		Vector x_r(x,0,ndof_f),x_i(x,ndof_f,ndof_f);
		Prolongation->AddMult(xc_r,x_r);
		Prolongation->AddMult(xc_i,x_i);
	}
}

//twolevel solver
twolevelsolver::~twolevelsolver(){}
twolevelsolver::twolevelsolver():coarsesolver(),finesolver(),FEMHierarchy(){}
twolevelsolver::twolevelsolver(int nrow,SparseMatrix O_c,SparseMatrix O_f,FiniteElementSpaceHierarchy *FEMH,Solver* csolver,Solver* fsolver):Solver(nrow,nrow)
{
	coarsesolver=csolver;
	finesolver=fsolver;
	FEMHierarchy=FEMH;
	A_c=O_c;
	A_f=O_f;
}
void twolevelsolver::SetOperator(const Operator& op)
{
   MFEM_ABORT("SetOperator is not supported in DDM!");
}
void twolevelsolver::solve(Vector& r, Vector& x,Vector& p)
{
	//PCG in finespace
	Vector Br_f(x.Size()),Ap_f(x.Size());
	double alpha_f,beta_f;
	finesolver->Mult(r,Br_f);
	A_f.Mult(p,Ap_f);
	double Brdr_f=Br_f*r;
	alpha_f=Brdr_f/(p*Ap_f);
	p*=alpha_f;
	x+=p;
	p*=(1/alpha_f);
	Ap_f*=alpha_f;
	r-=Ap_f;
	finesolver->Mult(r,Br_f);
	beta_f=(Br_f*r)/Brdr_f;
	p*=beta_f;
	p+=Br_f;
	//prolong in coarse mesh
	mfem::Operator *P=FEMHierarchy->GetProlongationAtLevel(0);
	int width=P->Width();
	Vector p_c(width),x_c(width),r_c(width);
	Vector Ap_c(width),Br_c(width);
	P->MultTranspose(x,x_c);
	P->MultTranspose(r,r_c);
	P->MultTranspose(p,p_c);
	// //PCG in coarse mesh
	// double alpha_c,beta_c;
	// coarsesolver->Mult(r_c,Br_c);
	// A_c.Mult(p_c,Ap_c);
	// double Brdr_c=Br_c*r_c;
	// alpha_c=Brdr_c/(p_c*Ap_c);
	// p_c*=alpha_c;
	// x_c+=p_c;
	// p_c*=(1/alpha_c);
	// Ap_c*=alpha_c;
	// r_c-=Ap_c;
	// coarsesolver->Mult(r_c,Br_c);
	// beta_c=(Br_c*r_c)/Brdr_c;
	// p_c*=beta_c;
	// p_c+=Br_c;
	//interpolation to fine mesh
	P->Mult(x_c,x);
	P->Mult(r_c,r);
	P->Mult(p_c,p);
	return;
}
void twolevelsolver::Mult(const Vector& b, Vector& x) const
{}
}
}

