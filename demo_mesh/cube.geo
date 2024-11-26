lc = 0.5;
/* lc = 1.0; */
hx = 1.0;
hy = 1.0;
hz = 1.0;
x0 = 0.0;
y0 = 0.0;
z0 = 0.0;

Point(1) = {x0, y0, z0, lc};
Point(2) = {x0 + hx, y0, z0, lc};
Point(3) = {x0 + hx, y0 + hy, z0, lc};
Point(4) = {x0, y0 + hy, z0, lc};
Point(5) = {x0, y0, z0+hz, lc};
Point(6) = {x0 + hx, y0, z0+hz, lc};
Point(7) = {x0 + hx, y0 + hy, z0+hz, lc};
Point(8) = {x0, y0 + hy, z0+hz, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {5, 6, 7, 8};
Line Loop(3) = {1, 10, -5, -9};
Line Loop(4) = {2, 11, -6, -10};
Line Loop(5) = {3, 12, -7, -11};
Line Loop(6) = {4, 9, -8, -12};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};

Surface Loop(1) = {1,2,3,4,5,6};
Volume(1) = {1};


Physical Surface(11) = {1,2,3};
Physical Surface(12) = {4,5,6};
Physical Volume(101) = {1};

Mesh.MshFileVersion = 2.2;
