# Box around a cylinder with sinusoidal walls
#
# Create an outer boundary with two vertical sides and two sinusoidal
# "wobbly" sides. Inner boundary is a circle with a refinement line
# placed behind the cylinder to capture the wake region.
#
# Keywords: Outer boundary, inner boundary, paramteric equations,
#           circular arc, manual refinement region

using HOHQMesh, GLMakie

# Create a new HOHQMesh model project. The project name
# "cylinder_sine_walls" will be the name of the mesh file
# saved in the directory "out".
cylinder_flow = newProject("cylinder_sine_walls", "out")

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

# Outer boundary curve chain is created to have counter-clockwise
# orientation, as required by HOHQMesh generator
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
