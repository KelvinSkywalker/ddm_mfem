import meshio
from fealpy.mesh.tetrahedron_mesh import TetrahedronMesh 


def ReadGmsh(file_name):
    gmesh = meshio.read(file_name,file_format = 'gmsh')
    node = gmesh.points[:,:2]
    cell = gmesh.cells_dict['tetra']
    mesh = TetrahedronMesh(node,cell)
    print('num of element:', mesh.ds.cell.shape[0])
    print('num of face:', mesh.ds.face.shape[0])
    print('num of edge:', mesh.ds.edge.shape[0])
    print('num of node:', mesh.node.shape[0])

def CreateMesh():
    domain = [0,1,0,1,0,1]
    n = 3
    mesh = TetrahedronMesh.from_box(domain,nx=n,ny=n,nz=n)
    print('num of element:', mesh.ds.cell.shape[0])
    print('num of face:', mesh.ds.face.shape[0])
    print('num of edge:', mesh.ds.edge.shape[0])
    print('num of node:', mesh.node.shape[0])

# CreateMesh()
ReadGmsh('./cube.msh')
