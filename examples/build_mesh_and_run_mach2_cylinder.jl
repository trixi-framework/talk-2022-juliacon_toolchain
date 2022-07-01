###############################################################################
#
# Use HOHQMesh to build the model curves and generate the mesh
#
###############################################################################

include("interactive_cylinder_with_sine_walls.jl")

###############################################################################
#
# Use Trixi to build the spatial discretization for the compressible Euler
# equations on the unstructured quadrilatreal mesh generated above.
# Use OridinaryDiffEq to perform time integration.
#
# As written in `elixir_euler_mach2_cylinder.jl`, this only initializes the
# simulation, i.e., `tspan = (0.0, 0.0)`. To change the final time, edit
# line 94 of this file, e.g., `tspan = (0.0, 0.5)`.
#
###############################################################################

@info "Run Mach 2 flow simulation with Trixi. Press enter to continue..."
readline()

include("elixir_euler_mach2_cylinder.jl")

###############################################################################
#
# Use Trixi2Vtk to postprocess the Trixi `.h5` solution files and save them
# in the directory `plot_files`. Can be visualized in external software
# like ParaView.
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
