# impedance3Dtest
this experiment solve the Maxwell Equations with impedance boundary conditions using RAS_IMP.
$$\nabla\times\nabla\times u - \kappa^2 u = f \text{ in }\Omega$$  
$$(\nabla\times u)\times n - i \kappa (n\times u)\times n = g \text{ on }\partial\Omega $$  
usage:  
make the .cpp code by make impedance3Dtest , and then run the experiment with argv as   
[1] refine_level [2] kappa [3]num_subdomain [4]delta (overlapping size), for example:  
./impedance3Dtest 1 3 4 1  
would run the experiment refining the mesh 1 times with kappa=3,4 sudmains and 1 overlapping size. 
