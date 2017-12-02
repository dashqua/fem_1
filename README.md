# fem_1 Stationnary code

Here I post one of my last work (at the 1st of December, 2017) on a the Finite Element Method.
My goal is to observe the behavior of a couple (wing , air flow) at subsonic velocity.
For simulation purpose, I will consider that the wing is fixed while the flow is moving around it. It is the strict contrary of the natural phenomenon but it as a reasonnable first attemp of a model.
To be able to pursue we have to dive in a little bit of physics consideration, I won't go into details but we have to keep in mind the following:
- The simulation is based on uncompressible Navier-Stokes equation but the solution is reduced to an harmonic function (displacement)
- We consider the airfoil to be in a squarebox. The upper and lower boundaries are considered as physically far away from the airfoil, the left boundary is where the flow is coming from and its trajectory is horizontal and inversely proportionnal to the height. The right boundary forces the solution to be the same as at the left boundary
- Speaking of solution, we computes here the streamline function
- About the wing, I am studying the NACA0012 airfoil which is available [here](http://airfoiltools.com/airfoil/details?airfoil=n0012-il)
- I use a collection of fortran tools, originally provided in my training and based of some fortran libs and GMSH meshtool. All the credits are mentionned in the code

OK let's move to the real stuff:
