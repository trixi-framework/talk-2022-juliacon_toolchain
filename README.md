# JuliaCon 2022: From Mesh Generation to Adaptive Simulation: A Journey in Julia

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)

<p align="center">
  <a href="https://www.youtube.com/watch?v=hoViWRAhCBE" target="_blank" rel="noopener noreferrer"><img
    src="https://user-images.githubusercontent.com/25242486/176432903-668ce8bf-4119-4d15-a46e-a1df90944e14.png"
    width="500px" /></a>
</p>

This is the companion repository for the [JuliaCon 2022](https://juliacon.org/2022) talk

**From Mesh Generation to Adaptive Simulation: A Journey in Julia**<br>
*Andrew R. Winters*<br>

(see abstract [below](#abstract)). Here you can find the presentation slides
in [talk.pdf](talk.pdf) as well as Julia scripts in the [examples](examples/)
directory to locally create results of mesh generation with
[HOHQMesh.jl](https://github.com/trixi-framework/HOHQMesh.jl)
and [Trixi.jl](https://github.com/trixi-framework/Trixi.jl) simulations
presented in the talk.

To reduce the size of the file for the [talk.pdf](talk.pdf), the video is not
embedded. Instead, the video `cylinder_flow_sine_walls.mp4` shown in the
presentation is located in [media](media/)
or it is available below
```@raw html
  <style>.embed-container { position: relative; padding-bottom: 56.25%; height: 0; overflow: hidden; max-width: 100%; } .embed-container iframe, .embed-container object, .embed-container embed { position: absolute; top: 0; left: 0; width: 100%; height: 100%; }</style><div class='embed-container'><iframe src='https://www.youtube.com/embed/Q3Pi41gbOkI' frameborder='0' allowfullscreen></iframe></div>
```

## Abstract

We present a Julia toolchain for the adaptive simulation of hyperbolic PDEs
such as flow equations on complex domains. It begins with using HOHQMesh.jl to
create a curved, unstructured mesh. This mesh is then used in Trixi.jl, a
numerical simulation framework for conservation laws. We visualize the
results using Julia’s plotting packages. We highlight select features
in Trixi.jl, like adaptive mesh refinement (AMR) or shock capturing,
useful for practical applications with complex transient behavior.


### More detailed description

Applications of interest in computational fluid mechanics typically occur
on domains with curved boundaries. Further, the solution of a non-linear
physical model can develop complex phenomena such as discontinuities,
singularities, and turbulence.

Attacking such complex flow problems may seem daunting. In this talk,
however, we present a toolchain with components entirely available in
the Julia ecosystem to do just that. In broad strokes the workflow is:

1. Use [HOHQMesh.jl](https://github.com/trixi-framework/HOHQMesh.jl)
   to interactively prototype and visualize a domain with curved boundaries.
2. HOHQMesh generates an all quadrilateral mesh amenable for high-order numerical
   methods.
3. The mesh file is passed to [Trixi.jl](https://github.com/trixi-framework/Trixi.jl),
   a numerical simulation framework for conservation laws.
4. Solution-adaptive refinement of the mesh within Trixi is handled by
   [P4est.jl](https://github.com/trixi-framework/P4est.jl).
5. After the simulation, interactive visualization can be done using
   [Makie.jl](https://makie.juliaplots.org/stable/).
6. Solution data can also be exported with
   [Trixi2Vtk.jl](https://github.com/trixi-framework/Trixi2Vtk.jl)
   for visualization in
   external software like [ParaView](https://www.paraview.org/).

The strength and simplicity of this workflow is through the combination
of several packages either originally written in Julia, like Trixi.jl,
or wrappers, like P4est.jl or HOHQMesh.jl, that provide Julia users access
to powerful, well-developed numerical libraries and tools written in other
programming languages.

## Getting started


### Installing Julia
To obtain Julia, go to https://julialang.org/downloads/ and download the latest
stable release (v1.7.3 as of 2022-06-28). Then, follow the
[platform-specific instructions](https://julialang.org/downloads/platform/)
to install Julia on your machine. Note that there is no need to compile anything
if you are using Linux, MacOS, or Windows.

After the installation, open a terminal and start the Julia *REPL*
(i.e., the interactive prompt) with
```shell
julia
```

### Installing the required Julia packages
To run the scripts in the [examples](examples/) directory and allow for
fully reproducible results, we have used Julia's package manager
to pin all packages to a fixed release. This makes it straightforward to
reproduce the Julia environment in which all the results presented were created.

If you have not done it yet, clone the repository where this code is stored:
```shell
git clone https://github.com/trixi-framework/talk-2022-juliacon_toolchain.git
```
Then, navigate to your repository folder and install the required packages:
```shell
cd talk-2022-juliacon_toolchain
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```
This will download and build all required packages, including the ODE package
[OrdinaryDiffEq](https://github.com/SciML/OrdinaryDiffEq.jl), the visualization
package [GLMakie](https://github.com/JuliaPlots/Makie.jl/tree/master/GLMakie),
the mesh generator [HOHQMesh.jl](https://github.com/trixi-framework/HOHQMesh.jl),
and [Trixi](https://github.com/trixi-framework/Trixi.jl).
The `--project=.` argument tells Julia to use the `Project.toml`
and `Manifest.toml` files from this repository to figure out which packages to install.

Once the initialization and installation is complete you must start Julia with the
`--project` flag set to your local clone of this repository
```shell
julia --project=@.
```


## Mesh created with tools from HOHQMesh.jl
To reproduce the figures and create the mesh file output execute from the REPL
```julia
include(joinpath("examples", "interactive_cylinder_with_sine_walls.jl"))
```
This will create the directory `out` where the mesh file is saved.

## Simulation with Trixi.jl
The elixir file described in the presentation to setup and run a simulation
of Mach 2 flow over a cylinder is `elixir_euler_mach2_cylinder.jl`.
To run the simulation up to a final time of 0.5, execute from the REPL
```julia
using Trixi
trixi_include(joinpath("examples", "elixir_euler_mach2_cylinder.jl"), tspan=(0.0,0.5))
```
where the final time is adjusted within the `trixi_include` call.
This simulation to the final time 0.5 takes approximately 20 minutes on a single thread.
To visualize the solution `sol` at the final time execute
```julia
using GLMakie
pd = PlotData2D(sol)
plot(pd["rho"])
```

## Combined script
The script `build_mesh_and_run_mach2_cylinder.jl` executes the entire toolchain
described in the talk. That is, the script generates the mesh, runs the simulation,
visualizes the final result in Makie, converts the output files to VTK format
using Trixi2Vtk, and saves them to the `plot_files` directory. Execute this script with
```julia
include(joinpath("examples", "build_mesh_and_run_mach2_cylinder.jl"))
```

As written, the script `elixir_euler_mach2_cylinder.jl` sets `tspan = (0.0, 0.0)` on line 94.
This can be adjusted
to take a different final time, e.g., the final time for the video in [media](media/)
is 2.25.

To reproduce the ParaView visualization, first open ParaView (after
[downloading and installing](https://www.paraview.org/download/) it if necessary).
Then load the ParaView state by clicking
through `File -> Load State` and open `supersonic_cylinder_state.pvsm.pvsm`.
Next, from the prompt "Load State Data File Options" select "Choose File Names",
navigate to the `plot_files` directory and select the appropriate
`solution.pvd` and `solution_celldata.pvd` files.


## Authors
This repository was initiated by
[Andrew Winters](https://liu.se/en/employee/andwi94).


## License
The contents of this repository are licensed under the MIT license
(see [LICENSE.md](LICENSE.md)).
