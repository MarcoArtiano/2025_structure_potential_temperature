# Numerical experiments
The test cases have been run with Julia v1.11.6. It is strongly recommended to start julia with multiple threads, as some simiulations take hours to be completed. To do that, start julia with the following command:
```bash
julia --project=. --threads=num
```
where `num` is the number of threads that you want to use for that session. Alternatively, you can also set `--threads=auto` to use a reasonable default number of threads for your system.

The equations are defined inside the folder `PotentialTemperature/`,  as long as other utilities function.

Before running the simulations, please activate the project and instantiate, i.e., start Julia as described above and run the following code in the Julia REPL:
```julia
julia> using Pkg
julia> Pkg.activate(".")
julia> Pkg.instantiate()
```

## Test case 1: Conservation of total energy and physical entropy

Run the following command to reproduce the results about 1D density wave test case
```julia
julia> include("test_cases/conservation/density_wave.jl")
```
Run the following command to reproduce the results about 3D Taylor-Green vortex test case
```julia
julia> include("test_cases/conservation/tgv.jl")
```
Run the following command to create the
```julia
julia> include("test_cases/conservation/post_process.jl")
```
## Test case 2: Well-balancedness
Run the following commands to reproduce the results about 2D well balancedness on curvilinear mesh
```julia
julia> include("test_cases/balance/isothermal_rhotheta.jl")

julia> include("test_cases/balance/isothermal.jl")

julia> include("test_cases/balance/potential_temperature_rhotheta.jl")

julia> include("test_cases/balance/potential_temperature.jl")
```
and
```julia
julia> include("test_cases/balance/plot.jl")
```
to create the figures shown in the paper.

## Test case 3: Inertia gravity waves
Run the following command to reproduce the results about inertia gravity waves
```julia
julia> include("test_cases/gravitywaves/inertia_gravity_waves.jl")
```
julia> include("test_cases/gravitywaves/convergence_analysis.jl")
```
and
```julia
julia> include("test_cases/gravitywaves/plots.jl")
```
to create the figures shown in the paper.

## Test case 4: Linear hydrostatic mountain
Run the following command to reproduce the results about linear hydrostatic mountain
```julia
julia> include("test_cases/linearhydrostatic/run_cases.jl")
```julia
julia> include("test_cases/linearhydrostatic/run_cases.jl")
```
and
```julia
julia> include("test_cases/linearhydrostatic/plots.jl")

julia> include("test_cases/linearhydrostatic/contours_hydrostatic.jl")
```
to create the figures shown in the paper.

## Test case 5: Linear nonhydrostatic mountain
Run the following command to reproduce the results about linear nonhydrostatic mountain
```julia
julia> include("test_cases/linearnonhydrostatic/run_cases.jl")
```
and
```julia
julia> include("test_cases/linearnonhydrostatic/plots.jl")

julia> include("test_cases/linearnonhydrostatic/contours_nonhydrostatic.jl")
```
to create the figures shown in the paper.

## Test case 6: Schär mountain
Run the following command to reproduce the results about Schär mountain
```julia
julia> include("test_cases/mountain/contours_mountain.jl")
```

## Test case 7: Baroclinic instability
Run the following command to reproduce the results about baroclinic instability
```julia
julia> include("test_cases/baroclinic/baroclinic.jl")
```
