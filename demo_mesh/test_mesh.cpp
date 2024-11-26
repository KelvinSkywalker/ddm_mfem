#include "mfem.hpp"
#include <fstream>
#include <iostream>
using namespace std;
using namespace mfem;

void NumInfo(const Mesh & mesh)
{
	cout<< "====================================="<<endl;
	cout << "dim: " << mesh.Dimension() << endl;
	cout << "has boundary: " << mesh.HasBoundaryElements() << endl;

	int NE = mesh.GetNE();
	int NF = mesh.GetNFaces();
	int Ne = mesh.GetNEdges();
	int NV = mesh.GetNV();
	cout << "number of elements: " << NE << endl;
	cout << "number of faces: " << NF << endl;
	cout << "number of edges: " << Ne << endl;
	cout << "number of nodes: " << NV << endl;
	cout << "GetNumFaces(): " << mesh.GetNumFaces() <<endl;

	for(int i=0;i<NE;i++)
	{
		int attr = mesh.GetAttribute(i);
		cout << "the attribute of element "<<i<< " is "<<attr<<endl;
	}
}

void ElemInfo(const Mesh & mesh)
{
	cout<< "====================================="<<endl;
	const Element * elem;
	elem = mesh.GetElement(1);
	int num_edge = elem->GetNEdges();
	cout << "The number of edges in Element 1 is " << num_edge << endl;

	const Element * face;
	face = mesh.GetFace(1);
	cout << "The number of vertices in Face 1 is "<< face->GetNVertices() << endl;

	const double * coord;
	coord = mesh.GetVertex(50);
	cout << "The coordinate of vertex 50 is: "<< *coord << ", "<< *(coord+1) << ", "<< *(coord+2) << endl;
}

void TopoInfo(Mesh & mesh)
{
	cout<< "====================================="<<endl;
	Table * v2e = mesh.GetVertexToElementTable();
	int nrow = v2e->Size();
	int * row = v2e->GetI();
	int * col = v2e->GetJ();
	cout << "the number of vertices is: "<< nrow <<endl;
	for(int i=0;i<nrow;i++)
	{
		int idx_begin = row[i];
		int idx_end = row[i+1];
		cout << "the index of element connected to vertex "<<i<<" is: ";
		for(int j=idx_begin;j<idx_end;j++)
		{
			int idx_elem = col[j];
			cout << idx_elem << " " ;
		}
		cout << endl;
	}

	delete v2e;
}

void BdInfo(const Mesh & mesh)
{
	cout<< "====================================="<<endl;
	int num_bd_elem = mesh.GetNBE();
	cout << "the number of boundary element (triangle): "<<num_bd_elem<<endl;
	
	int bd_attr_size = mesh.bdr_attributes.Size();
	cout << "mesh.bdr_attributes.Size() = " << bd_attr_size << endl;
	cout << "mesh.bdr_attributes.Max() = " << mesh.bdr_attributes.Max() << endl;
	for(int i=0;i<bd_attr_size;i++)
	{
		cout<< "bdr_attr["<<i<<"] = "<< mesh.bdr_attributes[i]<<endl;
	}
}
void MeshPartition(Mesh &mesh,int numpart)
{
	int *elemid=mesh.GeneratePartitioning(numpart);
	cout<<"the elem index list of Partition =";
	int NE=mesh.GetNE();
	for(int i=0;i<NE;i++)
	cout<<"the elem "<<i<<"is in subdomain"<<elemid[i]<<endl;
}

int main()
{
	const char *mesh_file = "./cube.msh";
	/* const char *mesh_file = "./beam-tet.mesh"; */
	Mesh mesh(mesh_file,1,1);

	NumInfo(mesh);
	ElemInfo(mesh);
	MeshPartition(mesh,3);
	TopoInfo(mesh);
	BdInfo(mesh);

	ofstream outfile;
	outfile.open("mesh.vtk");
	mesh.PrintVTK(outfile);
	outfile.close();
	/* mesh.PrintVTU("mesh"); */
	return 0;
}
