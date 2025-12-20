
using CairoMakie


include("Auxiliary/AuxHydro.jl")
include("Auxiliary/AuxRK4.jl")

# Initial proper time
tₛ = 0.1
tₑ = 10

# Initial conditions: (ε, n, Π, π)
u0 = [
    0.3,   # ε(τ0)
    0.1,    # n(τ0)
    0.0,    # Π(τ0)
    0.0     # π(τ0)
]

# Transport and EOS parameters
cₛ² = 1/3
η   = 0.2
ζ   = 0.0
τπ  = 0.5
τΠ  = 0.5

p = (cₛ², η, ζ, τπ, τΠ)


tspan = exp.(range(log(tₛ), log(tₑ), length=N))   # time step array

u = RK4(u0,tspan,H_IS!,p)

ε = u[:,1]
n = u[:,2]

T = @. ((π^2 )* (ε / 3))^(1/4)

τ = ((T).*(tspan)./((5)*η))  # \tau/\tau_R scaled time variable

plot(τ,T, 
        xaxis   =:log ,
        #yaxis   =:log ,
        xlabel  =   "τ",  
        ylabel  =   "T",
        dpi     =   300)

