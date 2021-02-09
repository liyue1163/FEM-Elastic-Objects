# FEM-Elastic-Objects
The finite element method for 3D-elastic objects 

## Compilation for Testing

Compiling the code in the above manner will yield working, but very slow executables. To run the code at full speed, you should compile it in release mode. Starting in the **build directory**, do the following:

    cmake .. -DCMAKE_BUILD_TYPE=Release
    
Followed by:

    make 
  
Your code should now run significantly (sometimes as much as ten times) faster. 

If you are using Windows, then running `cmake ..` should have created a Visual Studio solution file
called `a3-finite-elements-3d.sln` that you can open and build from there. Building the project will generate an .exe file.

Why don't you try this right now?

## Execution

Once built, you can execute the assignment from inside the `build/` using 

    ./a3-finite-elements-3d
 
Or 

    ./a3-finite-elements-3d arma

![Armadillo simulated via Finite Element Elasticity](img/armadillo.gif)


