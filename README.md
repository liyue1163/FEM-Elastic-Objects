# FEM-Elastic-Objects
The finite element method for 3D-elastic objects 

## Compilation for Testing

Compiling the code in the above manner will yield working, but very slow executables. To run the code at full speed, you should compile it in release mode. Starting in the **build directory**, do the following:

    cmake .. -DCMAKE_BUILD_TYPE=Release
    
Followed by:

    make 


## Execution

Once built, you can execute the assignment from inside the `build/` using 

    ./a3-finite-elements-3d
 
Or 

    ./a3-finite-elements-3d arma

![Armadillo simulated via Finite Element Elasticity](img/armadillo.gif)

double pressing `S` to switch between the underlying coarse model (default on load-up) and the detailed "skin" model.
