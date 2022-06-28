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
tspan = (0.0, 2.25)
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
