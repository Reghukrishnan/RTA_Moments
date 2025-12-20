using Plots
using ProgressMeter


include("Auxiliary/AuxAnalytical.jl")
include("Auxiliary/AuxMB.jl")
include("Auxiliary/AuxRK4.jl")

T₀  = 1             #  GeV     
m    = 1.01      # m corresponds to m = 0.1  Gev and T = 1GeV
α₀   = 0 

tₛ  =  0.1           # τ fm
tₑ  = 10
τ₀  = 0.1           # τ/τᵣ

N = 1000    # Number of time steps
tspan = exp.(range(log(tₛ), log(tₑ), length=N))   # time step array

n = 1
lmax = 50  # How many l moments to take




l = 0
g0 = [ ρeq(T₀,α₀,n,l+k) for k in 0:lmax] #Inital l moments
ρnl0 = Vector{Float64}(undef, N)         # Array to which time evolution is stored

Analytical_BFS!(ρnl0,tspan,g0,n,l,lmax)


l = 1      # Which l moment are we computing
g1 = [ ρeq(T₀,α₀,n,l+k) for k in 0:lmax] #Inital l moments
ρnl1 = Vector{Float64}(undef, N) 

Analytical_BFS!(ρnl1,tspan,g1,n,l,lmax)

plot(tspan,ρnl1./ρnl0, 
        xaxis   =:log ,
        #yaxis   =:log ,
        xlabel  =   "τ",  
        ylabel  =   "ρnl",
        dpi     =   300)




