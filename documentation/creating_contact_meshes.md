1) Read [Smith2018ContactMesh](https://github.com/opensim-jam-org/jam-resources/blob/main/documentation/doxygen/Smith2018ContactMesh_doxygen.pdf) and [Smith2018ArticularContactForce](https://github.com/opensim-jam-org/jam-resources/blob/main/documentation/doxygen/Smith2018ArticularContactForce_doxygen.pdf) doxygen to understand the mesh requirements. This [Smith et al, CMBBE, 2018](https://www.tandfonline.com/doi/full/10.1080/21681163.2016.1172346) paper gives further details on the implementation of the contact model.
2) Convert geometry to a triangular mesh (.stl, .obj or .vtp) file type
3) Transform the mesh to the correct reference frame. This could be the OpenSim parent body reference frame to which the contact mesh will be attached. You can also use a local frame frame and attach the mesh to a PhysicalOffsetFrame located on a body in the OpenSim model. See [DM grand challenge model](https://github.com/opensim-jam-org/jam-resources/blob/main/models/knee_tka/grand_challenge/DM/DM.osim) for an example.  
   * The Transform filter in paraview is helpful for translating, rotating, and scaling meshes.
4) Remesh (if necessary) to generate a clean, smooth, mesh with equilateral triangles of similar areas. You many also need to adjust the number of triangles in the mesh. Meshes with huge numbers of triangles (100,000+) will cause performance and visualization issues. 
   * I find the "Remesh" and "Reduce" tools in MeshMixer to be the best for this step.
   * See Smith CMBBE 2018 for convergence study on triangle size.     
5) Remove extra triangles that are not on the contact surface. This will improve computational performance.
6) Check normal vectors. Remesh contact surface mesh if necessary. 
7) Visualize contact mesh on top of the original mesh geometry to check that the mesh processing 
