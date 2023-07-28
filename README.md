# MongeAmpereDG

1.Installation
dune-mongeampere has the following dependencies
• cimg header
• eigen3 library
• Umfpack
• boost program options
• adolc
• opennurbs

If those are not found in standard paths, further search paths may be added in the files in dune-
mongeampere/cmake/modules.

After all dependecies are installed, you execute the following command in the folder dune-common
to execute cmake and make in all dune modules (compare dune installation )
@. . /dune−common/bin: ./dunecontrol configure
@. . /dune−common/bin: ./dunecontrol all
Optionally you may configure install folder and compiler in an option file - compare the dune
website and the help of dunecontrol on the option --opts=FILE. If the compilation of some dune
modules fails, it is possible to reinit the compilation process of a single target with the parameter
--only=dune-MODULENAME. Note that all included dune modules need to be installed. If the
configure command was executed successfully, the compile command ”make” can be executed in
the local build directores.
1.1 Quickstart
After a successful build you should find the program dune mongeampere image OT in the src folder
of the building directory. It is called with two further configuration files. The control parameters
of the finite element (FE) solver are given by the .ini file which is passed iwth the command line
option -c. The configuration ini file of the problems geometry however is passed by the command
line option -g. In the folder inputData are two of example init files given
Before calling the program it has to be secured that all output folders exist. E.g. in the case of
the example ini files .
@ . . dune−mongeampere/build-cmake: mkdir data
@ . . dune−mongeampere/build-cmake: mkdir data/PS12Splines
@ . . dune−mongeampere/build-cmake: mkdir data/PS12Splines/OT
@ . . dune−mongeampere/build-cmake: mkdir plots
@ . . dune−mongeampere/build-cmake: mkdir plots/PS12Splines
@ . . dune−mongeampere/build-cmake: mkdir plots/PS12Splines/OT

In the standard build folder you can call the program then with
@ . . dune−mongeampere/build-cmake/src: ./dune_mongeampere_image_OT
−c ../../inputData/SolversquareToSquareOT.ini −g ../../inputData/imageOTsetting.ini

If the build folder was chosen manually, the paths need to be adapted accordingly.

2.executable Programs
dune_mongeampere : solves standard MA equations
dune_mongeampere_OT : solves the in Optimal Transport (OT) arising MA equation det(D2 (u)) − f/g(∇u)= 0 with the sensible transport boundary conditions. Input distributions f and g are fixed in the source code.
dune_mongeampere_image_OT : solves the in Optimal Transport (OT) arising MA equation det(D2 (u)) − f/g(∇u)= 0 with the sensible transport boundary conditions. Input of the distributions f and g are image files
dune_mongeampere_reflector : calculates a reflective surface considering a point light source
dune_mongeampere_reflector_parallel : calculates a reflective surface considering parallel light
dune_mongeampere_refractor : calculates a refractive surface considering a point light source
dune_mongeampere_refractor parallel : calculates a refractive surface considering parallel light
read_and_write_OT : possibility to read the FE solution given its coefficients, in order to generate plots
export_refractor_points : generates points on a cartesian grid given a FE solution
fit_surface : fits a NURBS surface given a cartesian point grid
convert_dat_to_bmp : convertes a .data file to a .bmp file

3. Parameters for executables starting with dune_mongeampere_...

All necessary input files are listed by callling the command
@ . . dune−mongeampere/build−cmake/src:./dune_mongeampere_... −h

Usually one configuration file for the geometry, passed by -g and one file for the FE solver, passed by -c, is required.

3.1 Configuration

Note that the some configurations have to be chosen before compiling: In the source code you can easily change
1. choice of the finite element space
2. using automatic differentiation
3. parallel light or point source
4. optimisation in case that the source domain Ω is rectangular
For more information see appendix A. All other parameters are given during the run time by the
ini files.

23.1.1
Parameter in the FE solver ini file (passed with -c)
[solver]
initValueFromFile: indicates if the initialguess is read from file
initValue: path to initialguess (given in ASCII coefficient values)
startlevel: start refinement of initial grid
nonlinearSteps: steps of refinement
maxSteps: maximum number of steps in a nonlinear step
Dirichlet: indicate Dirichlet boundary conditions (deprecated)
JacSimultaneously: indicate combined evaluation of function and Jacobian
lambda: penalisation value for determinant
sigma: penalisation value for jumps in DG method (deprecated)
sigmaGrad: penalisation value for gradient jumps in ansatz (deprecated)
sigmaBoundary: penalisation value for weakly forced boundary conditions (deprecated)
[dogled]
iradius: initial trust region radius
stopcriteria0: inf_norm(f') ≤ stopcriteria(1)
stopcriteria1: 2_norm(dx) <= stopcriteria(2) ∗ (stopcriteria(2) + 2_norm(x) )
stopcriteria0: inf_norm(f) <= stopcriteria(3)
silendmode: suppress output of dogleg solver
exportJacobianIfSingular: activate matlab export to std::cout in case of singular Jacobian
check_Jacobian: activate a check of the Jacobian with Finite Differences (computationally expensive)
exportJacobianifFalse: activate matlab export to std::cout in case of wrong Jacobian
[output]
directory: output folder for data as coefficients
plotdirectory: output folder for output data as vtks, rayTracer files, 3dm files ...
prefix: output prefix
refinement: specify additional grid refinement for plots
cartesianGridN: specify the elements for a cartesian export grid
write_vtk: write solution to vtk-ASCII files

3.1.2 Parameter of the geometry ini file (passed with -g)
[geometry optic]
xMin: lower left x value of optical surface
xMax: upper right x value of optical surface
yMin: lower left y value of optical surface
yMax: upper right y value of optical surface
gridfile: path to initial grid file (.msh format)
plotgridfile: path to an additional plotgrid file (.msh format)
boundaryN: number of direction in approximation of source boundary
initialOpticDistance: distance between source and opt.surface for the initial guess
[geometry target]
xMin: lower left x value of targetet
xMax: upper right x value of target
yMin: lower left y value of target
yMax: upper right x value of target
gridfile: path to grid of the target (.msh format)
target_is_xy_plane: indicates if target is parallel to xy plane (otherwise parallel to xz is assumed)
z: z value of the target
boundaryN: number of direction in approximation of target boundary
[light in]
imageName: path to image of source distribution
[lightout out]
targetImageName: path to image of target distribution
[solver]
epsDivide: controls blurring between steps
epsEnd: controls blurring in the final step
minPixelVale: mininmal pixel value the target is brightened up to
[povray]
cameraAngle: camera angle for ray tracer config file
jitter: jitter for ray tracer config file
nPhotons: number of photons for ray tracer config file
lightSourceRadius: light source radius for ray tracer config file
lightSourceFalloff: fall of of light ray for ray tracer config file
lightSourceTightness: light source tightness for ray tracer config file
lightSourceIntensity: light source intensity for ray tracer config file
writeAperture: indicate whether light source is additionally cropped in ray tracer file 

4. output files
Since the FE solver solves the problem on different grid with different target distributions, output is generated of every intermediate step and indicated by a counter. Imagine we are in the second intermediate step
If not specified otherwise in the ini file the solution u is written to a vtk-file called output folder/EXAMPLE2numericalSolu
in xml format. For more information about those plots see B.

4.1 Output in the OT case
Additionally the transport of the source grid Ω, i.e. nablau(Ω), is written in another .vtk-Datei
called output folder/EXAMPLE2outputGrid.vtu.

4.2 Output of the optical surfaces
The file output folder/EXAMPLEnumericalSolutionreflector2.pov can be used to render the target under the current optical surface with the ray tracer povray.
The file output folder/EXAMPLEnumericalSolutionreflector2.3dm exports a triangular mesh of the current optical surface.
The file output folder/EXAMPLEnumericalSolutionreflectorPoints2.3dm exports a point cloud of the optical surface.

4.3 Conversion in 3dm NURBS files
Note that output folder/EXAMPLEnumericalSolutionreflector2.3dm is only an linear approxima-
tion on the actual optical surface. To export a better approximation EXAMPLEnumericalSolu-
tionrefractorPoints4.3dm in case of a rectangular grid do the following.
@. . /src: ./exportrefractorpoints −g ../inputData/Optic/lensSetting.ini −o ../inputData/Optic/ExportOptions.ini
@. ./src: ./fitsurface ../plots/Lens/EXAMPLEnumericalSolutionrefractorPoints4.txt ../plots/Lens/EXAMPLEnumericalSolutionrefractorPoints4.3dm

This generates a Cartesian point cloud of the example and fits a cubic NURBS surface.

Apppendix A : Parameter during compile time

A.1 Choice oft the finite elements
The desired finite element is enabled per a preprocessor variable. Up-to-date is probably only the variable PS12, which enables C 1 -continuous S-Splines based on the Powell-Sabin-12-Split.

A.2 Automatic differentiation
Calculations without automatic differentiation by adolc is not possible at the moment ... In the implementation with known derivative (i.e. only in the OT case) there is some bug ...

A.3 Kind of light source
To change the point light source to a source with parallel light, the preprocessor variable PARALLEL LIGHT has to be defined. This however, is not made manually, but the cmake file compiles both versions separately.

A.4 Rectangular grid
To allow some simplifications and optimisation in case of a rectangular source grid you should define in the file solver/solver config.h die preprocessor variable RECTANGULAR GRID.

Appendix B:  Information in the vtk file

output folder/EXAMPLE2numericalSolution.vtu contains not only the values of u, but also the additional attributes
gradx : ∇u
HessianIJ : entries of the Hessian D2_uIJ
Residual : the residual (e.g. det(D2 u) − f/g(∇u) in the OT case)
EV1/2 : eigen values of D^2 u
Curvature: the Gaussian curvature of the surface (only for in case of optics)

In the OT case output folder/EXAMPLE2outputGrid.vtu includes the attribut est integral:

In every triangle ∆ of the grid of Ω the integral 1/2 int_∆ f dx was estimated by taking the average of f at the corners and multiplying by the area |∆|. 
With the OT problem property it holds that the integral 1/2 int_tau(∆) g dy is the same. Hence we store the value 1/2 int_∆ f dx/ |tau(∆)| which
should resemble g for small ∆.
