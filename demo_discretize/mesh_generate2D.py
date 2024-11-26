import os
import gmsh
import math
def cavity_mesh(k, lc, filename, non_uniform=False):
    pi = math.pi
    leftx = 0
    rightx = 0.36+4*2*pi/k
    bottomy = 0
    uppery = 0.08+4*2*pi/k
    gmsh.initialize()
    model = gmsh.model
    # lc = pi/(10*k)
    # lc = k**(-1.5)
    # 定义外部矩形
    p1 = model.geo.addPoint(0, 0, 0, lc)
    p2 = model.geo.addPoint(rightx, 0, 0, lc)
    p3 = model.geo.addPoint(rightx, uppery, 0, lc)
    p4 = model.geo.addPoint(0, uppery, 0, lc)
 
    l1 = model.geo.addLine(p1, p2, 1)
    l2 = model.geo.addLine(p2, p3, 2)
    l3 = model.geo.addLine(p3, p4, 3)
    l4 = model.geo.addLine(p4, p1, 4)
 
    ll = [l1, l2, l3, l4]
    ll = model.geo.addCurveLoop(ll)
    # 定义内部矩形
    sp1 = model.geo.addPoint(2*2*pi/k, 2*2*pi/k, 0, lc)
    sp2 = model.geo.addPoint(0.36+2*2*pi/k, 2*2*pi/k, 0, lc)
    sp3 = model.geo.addPoint(0.36+2*2*pi/k, 0.08+2*2*pi/k, 0, lc)
    sp4 = model.geo.addPoint(2*2*pi/k, 0.08+2*2*pi/k, 0, lc)

    sl1 = model.geo.addLine(sp2, sp1, 5)
    sl2 = model.geo.addLine(sp3, sp2, 6)
    sl3 = model.geo.addLine(sp4, sp3, 7)
    sl4 = model.geo.addLine(sp1, sp4, 8)
 
    sll = [sl1, sl4, sl3, sl2]
    sll = model.geo.addCurveLoop(sll)
 
    model.geo.addPlaneSurface([ll,sll], 1)
    model.geo.synchronize()
 
    # 设置非一致剖分
    if non_uniform:
        model.mesh.setSize([(0, sp1)], lc/3)  # 在内部顶点处加密
        model.mesh.setSize([(0, sp4)], lc/3)
    gmsh.model.addPhysicalGroup(1, [l1], 2, name="bottom")
    gmsh.model.addPhysicalGroup(1, [l2], 3, name="right")
    gmsh.model.addPhysicalGroup(1, [l3], 4, name="top")
    gmsh.model.addPhysicalGroup(1, [l4], 5, name="left")
    gmsh.model.addPhysicalGroup(1, [sl1, sl2, sl3, sl4], 1, name="inner")

    gmsh.model.addPhysicalGroup(2, [1], name="My surface")

    model.mesh.generate(2)
    gmsh.write(filename)
    gmsh.finalize()
    
def square(n):
    gmsh.initialize()
    model = gmsh.model
    lc = 1/n
    a = 1
    b = 1
    p1 = model.geo.addPoint(0, 0, 0, lc)
    p2 = model.geo.addPoint(a, 0, 0, lc)
    p3 = model.geo.addPoint(a, b, 0, lc)
    p4 = model.geo.addPoint(0, b, 0, lc)
    l1 = model.geo.addLine(p1, p2, 1)
    l2 = model.geo.addLine(p2, p3, 2)
    l3 = model.geo.addLine(p3, p4, 3)
    l4 = model.geo.addLine(p4, p1, 4)
    ll = model.geo.addCurveLoop([l1, l2, l3, l4])
    model.geo.addPlaneSurface(ll, 1)
    model.geo.synchronize()
    model.mesh.generate(2)
    gmsh.write('square' + str(n) +'.msh')
    gmsh.finalize()
 
if __name__ == "__main__":
    k = 10
    filename = "mycavity_"+str(k) + ".msh2"
    cavity_mesh(k, k**(-1.5), filename)
    os.rename(filename, "mycavity_"+str(k) + ".msh")    

    
    
    
