# mesh_converter

## Function 

This program is a simple utility to 

1. convert Gambit Neutral mesh file to VTK format files;
2. convert `.msh` file generated by GMSH to VTK format files;
3. convert low-order elements to high-order elements.

Meshes generated are compatible with [svFSI](https://github.com/SimVascular/svFSI).

## Build

Before building `mesh_converter`, make sure your system has `make` and `gfortran` installed. In Linux system, simply type
   ```base
   cd mesh_converter && make
   ```
If successful, an executable `convertMesh.exe` will be present in the folder `./bin`.

## Usage

The general usage of this program is
   ```bash
   <path to bin>/bin/convertMesh.exe <input-file>
   ```
Depending on the format of the original mesh(es), `<input-file>` can be either the name of the original mesh or a text file listing the paths to the original meshes. The output will include a body mesh named `mesh-complete.mesh.vtu` and all the boundary meshes stored in the folder `mesh-surfaces`.

Multiple examples are included here to further demonstrate its usage.

1. Convert Gambit Neutral mesh file
   ```bash
   cd ./example/01-tri3 && ../../bin/convertMesh.exe disk_tri_h0.02.neu
   ```

2. Convert GMSH file to VTK mesh
   ```bash
   cd ./example/11-quad4-gmsh && ../../bin/convertMesh.exe gmsh_flow_past_cylinder.msh
   ```
   Type 2 when asked for the dimension of the body mesh. When generating mesh using GMSH, remember to assign physical groups to both the body mesh and the boundary meshes.

3. Convert low-order elements to high-order elements

   ```bash
   cd ./example/08-tet10 && ../../bin/convertMesh.exe mesh_tet4.txt
   ```

   
