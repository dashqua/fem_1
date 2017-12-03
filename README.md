# fem_1 Stationnary code (Work in progress)

Here I post one of my last work (at the 1st of December, 2017) on a the Finite Element Method.
My goal is to observe the behavior of a couple (wing , air flow) at subsonic velocity.
For simulation purpose, I will consider that the wing is fixed while the flow is moving around it. It is the strict contrary of the natural phenomenon but it as a reasonnable first attemp of a model.
To be able to pursue we have to dive in a little bit of physics consideration, I won't go into details but we have to keep in mind the following:
- The simulation is based on uncompressible Navier-Stokes equation but the solution is reduced to an harmonic function (displacement)
- We consider the airfoil to be in a squarebox. The upper and lower boundaries are considered as physically far away from the airfoil, the left boundary is where the flow is coming from and its trajectory is horizontal and inversely proportionnal to the height. The right boundary forces the solution to be the same as at the left boundary
- Speaking of solution, we computes here the streamline function
- About the wing, I am studying the NACA0012 airfoil which is available [here](http://airfoiltools.com/airfoil/details?airfoil=n0012-il)
- I use a collection of fortran tools, originally provided in my training and based of some fortran libs and GMSH meshtool. All the credits are mentionned in the code

I modified some parts of the code and tidied the directory. It is still not as modular as I want it to be because in one hand the module management (library + some other ones..) are a bit messy for me and in another hand the code has still to be modified to handle different mesh options. Those are axis of improvment, feel free to do so.

## Reqs

GMSH  
MatLab (not Octave)  
gfortran  
libs already included   

### Build mesh

```
gmsh wing_naca.geo
```
then just mesh 2D, no refining or 2nd order polynomials.
File -> Save Mesh
If you are curious you would maybe sneak into the .geo, Important options are commented (fr) and you will note that the NACA is defined with its curve parametric equation. Since the data available on the website of NACA airfoils are sets of points, another axis of improvement would be to replace the equation by the set of points available on the website.  
Modifiation of options are not yet supported.

### Compile

```
make  
```

### Run

```
./pw2mts.exe
```

It generates solutionnaca.dat which is the solution file, obviously.

### Post processing

Fire up MatLab and run trace_naca.m to draw trajectories.








