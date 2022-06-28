###############################################################################
#
# Use HOHQMesh to build the model curves and generate the mesh
#
###############################################################################

# Wobbly box around a cylinder
#
# Create an outer boundary with two vertical sides and two sinusoidal
# "wobbly" sides. Inner boundary is a circle with a refinement line
# placed behind the cylinder to capture the wake region.
#
# Keywords: Outer boundary, inner boundary, paramteric equations,
#           circular arc, manual refinement region

using HOHQMesh, GLMakie

# Create a new HOHQMesh model project. The project name
# "box_around_circle" will be the name of the mesh file
# saved in the directory "out".
cylinder_flow = newProject("box_around_circle", "out")

# Reset polynomial order of the mesh model curves and output format.
# The "ABAQUS" mesh file format is needed for the adaptive mesh
# capability of Trixi.jl.
setPolynomialOrder!(cylinder_flow, 3)
setMeshFileFormat!(cylinder_flow, "ABAQUS")

# A background grid is required for the mesh generation. In this example we lay a
# background grid of Cartesian boxes with size 0.4.
addBackgroundGrid!(cylinder_flow, [0.4, 0.4, 0.0])

# Add outer boundary curves in counter-clockwise order.
# Note, the curve names are those that will be present in the mesh file.
left = newEndPointsLineCurve("Left", [0.0, 1.0, 0.0], [0.0, -1.0, 0.0])

xEqn = "x(t) = 4*t"
yEqn = "y(t) = -1.0 - 0.04*sin(8*pi*t)"
zEqn = "z(t) = 0.0"
bottom = newParametricEquationCurve("Bottom", xEqn, yEqn, zEqn)

right = newEndPointsLineCurve("Right", [4.0, -1.0, 0.0], [4.0, 1.0, 0.0])

xEqn = "x(t) = -4*(t-1)"
yEqn = "y(t) = 1.0 - 0.04*sin(8*pi*(t-1))"
zEqn = "z(t) = 0.0"
top = newParametricEquationCurve("Top", xEqn, yEqn, zEqn)

addCurveToOuterBoundary!(cylinder_flow, bottom)
addCurveToOuterBoundary!(cylinder_flow, right)
addCurveToOuterBoundary!(cylinder_flow, top)
addCurveToOuterBoundary!(cylinder_flow, left)

# Add inner boundary curve
cylinder = newCircularArcCurve("Circle",        # curve name
                               [0.6, 0.0, 0.0], # circle center
                               0.25,            # circle radius
                               0.0,             # start angle
                               360.0,           # end angle
                               "degrees")       # angle units

addCurveToInnerBoundary!(cylinder_flow, cylinder, "inner1")

# Add a refinement line for the wake region.
refine_line = newRefinementLine("wake_region", "smooth", [0.3,0.0,0.0], [4.0,0.0,0.0], 0.1, 0.25)

addRefinementRegion!(cylinder_flow, refine_line)

# Visualize the model, refinement region and background grid
# prior to meshing.
plotProject!(cylinder_flow, MODEL+REFINEMENTS+GRID)

@info "Press enter to generate the mesh and update the plot."
readline()

# Generate the mesh. Saves the mesh file to the directory "out".
generate_mesh(cylinder_flow)

###############################################################################
#
# Use Trixi to build the spatial discretization for the compressible Euler
# equations on the unstructured quadrilatreal mesh generated above.
# Use OridinaryDiffEq to perform time integration.
#
###############################################################################

# Channel flow around a cylinder at Mach 2
#
# Boundary conditions are supersonic Mach 2 inflow at the left portion of the domain
# and supersonic outflow at the right portion of the domain. The top and bottom of the
# channel as well as the cylinder are treated as Euler slip wall boundaries.
# This flow results in strong shock refletions / interactions as well as Kelvin-Helmholtz
# instabilities at later times as two Mach stems form above and below the cylinder.
#
# Keywords: supersonic flow, shock capturing, AMR, unstructured curved mesh,
#           positivity preservation, compressible Euler, 2D

@info "Run Mach 2 flow simulation with Trixi. Press enter to continue..."
readline()

using OrdinaryDiffEq
using Trixi

###############################################################################
# Semidiscretization of the compressible Euler equations

# Setup the governing equations, initial condition and boundary conditions
equations = CompressibleEulerEquations2D(1.4)

@inline function initial_condition_mach2_flow(x, t, equations::CompressibleEulerEquations2D)
  # set the freestream flow parameters
  rho = 1.4
  v1 = 2.0
  v2 = 0.0
  p = 1.0

  prim = SVector(rho, v1, v2, p)
  return prim2cons(prim, equations)
end

# Supersonic inflow boundary condition.
# Calculate the boundary flux entirely from the external solution state, i.e., set
# external solution state values for everything entering the domain.
@inline function boundary_condition_supersonic_inflow(u_inner, normal_direction::AbstractVector, x, t,
                                                      surface_flux_function, equations::CompressibleEulerEquations2D)
  u_boundary = initial_condition_mach2_flow(x, t, equations)
  flux = Trixi.flux(u_boundary, normal_direction, equations)

  return flux
end


# Supersonic outflow boundary condition.
# Calculate the boundary flux entirely from the internal solution state. Analogous to supersonic inflow
# except all the solution state values are set from the internal solution as everything leaves the domain
@inline function boundary_condition_supersonic_outflow(u_inner, normal_direction::AbstractVector, x, t,
                                                       surface_flux_function, equations::CompressibleEulerEquations2D)
  flux = Trixi.flux(u_inner, normal_direction, equations)

  return flux
end


# Create the unstructured quad mesh from the mesh file that
# was created by the interactive HOHQMesh.jl.
mesh_file = joinpath(@__DIR__, "out", "box_around_circle.inp")
mesh = P4estMesh{2}(mesh_file)

# Boundary names are those we assigned in HOHQMesh.jl
boundary_conditions = Dict( :Bottom  => boundary_condition_slip_wall,
                            :Circle  => boundary_condition_slip_wall,
                            :Top     => boundary_condition_slip_wall,
                            :Right   => boundary_condition_supersonic_outflow,
                            :Left    => boundary_condition_supersonic_inflow )


# Set the numerical fluxes for the volume and the surface contributions.
volume_flux = flux_ranocha_turbo
surface_flux = flux_lax_friedrichs

# Create the spatial discretization that uses SBP flux differencing ansatz.
# Shock capturing inidicator and combined with special volume integral blends
# low-order and high-order methods near discontinous (shocked!) flow regions.
polydeg = 4
basis = LobattoLegendreBasis(polydeg)
shock_indicator = IndicatorHennemannGassner(equations, basis,
                                            alpha_max=0.5,
                                            alpha_min=0.001,
                                            alpha_smooth=true,
                                            variable=density_pressure)
volume_integral = VolumeIntegralShockCapturingHG(shock_indicator;
                                                 volume_flux_dg=volume_flux,
                                                 volume_flux_fv=surface_flux)
solver = DGSEM(polydeg=polydeg, surface_flux=surface_flux, volume_integral=volume_integral)


# Combine all the spatial discretization components into a high-level descriptions.
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_mach2_flow, solver,
                                    boundary_conditions=boundary_conditions)

###############################################################################
# Setup an ODE problem
tspan = (0.0, 0.0) # For testing just initialize everything
# tspan = (0.0, 2.25)
ode = semidiscretize(semi, tspan)

# Callbacks
summary_callback = SummaryCallback()

# Prints solution errors to the screen at check-in intervals.
analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval=analysis_interval)

# Prints information to the screen to show that a simulation is still running.
alive_callback = AliveCallback(analysis_interval=analysis_interval)

# Saves the solution and mesh data to HDF5 files for postprocessing.
save_solution = SaveSolutionCallback(interval=1000,
                                     save_initial_solution=true,
                                     save_final_solution=true,
                                     solution_variables=cons2prim)

# Indicator and control of how P4est and Trixi handle adaptive mesh refinement.
amr_indicator = IndicatorLöhner(semi, variable=Trixi.density)

amr_controller = ControllerThreeLevel(semi, amr_indicator,
                                      base_level=0,
                                      med_level=3, med_threshold=0.05,
                                      max_level=5, max_threshold=0.1)

amr_callback = AMRCallback(semi, amr_controller,
                           interval=1,
                           adapt_initial_condition=true,
                           adapt_initial_condition_only_refine=true)

# Combine all the callbacks into a description.
callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        save_solution,
                        amr_callback)

# Positivity limiter for density and pressure necessary with strong shocks.
stage_limiter! = PositivityPreservingLimiterZhangShu(thresholds=(5.0e-7, 1.0e-6),
                                                     variables=(pressure, Trixi.density))

###############################################################################
# Run the simulation using Strong Stability Preserving Runge-Kutta (SSPRK).
sol = solve(ode, SSPRK43(stage_limiter!),
            save_everystep=false, callback=callbacks);

# Print the timer and simulation summary
summary_callback()

###############################################################################
#
# Use Trixi2Vtk to postprocess the Trixi `.h5` solution files and save them
# in the directory `out`. Can be visualized in external software like ParaView.
#
###############################################################################

@info "Postprocess output files with Trixi2Vtk. Press enter to continue..."
readline()

using Trixi2Vtk

trixi2vtk("out/solution*", output_directory="plot_files")

###############################################################################
#
# Use Makie (already loaded during the mesh generation step) to create a first
# visualization of the Trixi solution.
#
###############################################################################

@info "Visualize the solution with Makie."

pd = PlotData2D(sol);
plot(pd["rho"])
