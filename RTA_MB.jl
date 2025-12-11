using Plots
using ProgressMeter


include("AuxRK4.jl")
include("AuxMB.jl")


#-----------SCALED INITIAL CONDITIONS---------------------

T₀  = 1             #  GeV     
m    = 0.001         # m corresponds to m = 0.1  Gev and T = 1GeV
α₀   = 0            # μ/T 

r    = 0            # MB =0, FD = 1, BE = -1


tₛ  =  0.1           # τ fm
tₑ  = 100
τ₀  = 0.1           # τ/τᵣ

η₀      =   (tₛ/τ₀)*(T₀/5)  # viscosity to entropy ratio
#---------------------------------------------------
τᵣ⁰     =   (η₀/T₀)*(5)
ωᵣ⁰     =   (1/5)*(T₀/η₀)  # τᵣ⁰ = 5η₀/T₀. --> ωᵣ⁰ = 1/τᵣ⁰ = T₀/5η₀ --> ωᵣ = ωᵣ⁰(T/T₀)



#----------------
nₘₐₓ = 5 
nₐᵣ = [i for i in -1:nₘₐₓ]
nₙ = findall(x->x ==0,nₐᵣ)[1]   # Finding the location of 'n' for number desnity
nₑ = findall(x->x ==1,nₐᵣ)[1]   # Finding the lcoation of 'n' for energy density

N = size(nₐᵣ)[1]
L = 100

#------------------ Initial desnity-------------------------------------------------
ρ₀ = Matrix{Float64}(undef, N, L+1)

#Init_χ_Eq_q!(ρ₀,nₐᵣ,L,m,T₀,α₀,r)   # Initialises the moments with equilibrium initial conditions.
Init_ρ_Eq_b!(ρ₀,nₐᵣ,L,m,T₀,α₀)


ρ₀[nₙ,L+1] = α₀
ρ₀[nₑ,L+1] = T₀
#println((1/3)*ρ₀[nₑ,1]/ρ₀[nₙ,1])

#------------------Printing the parameter values------------------------------

if r ==  0 
    println("Maxwell Boltzmann moments")
elseif r ==  1
    println("Fermi Dirac moments")
elseif r == -1
    println("Bose Einstein moments")
end

println("m     : ", m)
println("T₀    : ",T₀)
println("t₀    : ",tₛ)
println("τ₀    : ",τ₀)
println("η₀/s₀ : ",η₀)


Nₚ = 10000                        # Number of time steps
ProgressBar = Progress(Nₚ)
tspan   = trange2((tₛ,tₑ),Nₚ,"exp")  # Exponentially scaled time steps


#------------------------------------------------
println("\nSolving Scaled RTA.")
p = (N,L,ωᵣ⁰,nₐᵣ,nₙ,nₑ,1)

ρ = RK4(ρ₀,tspan,RTAB,p)     # Rk4 solver for differential equations. Returns the moments as a function of time χ[t,n,l]


#---------------------------------------------------
T = ρ[:,nₑ,L+1]
α = @. log((ρ[:,nₙ,1]*((2*π^2)))/((T^3)*Γ(3)))
τ = ((T).*(tspan)./((5)*η₀))  # \tau/\tau_R scaled time variable

#neq = @.   exp(α)*(T^3)/((π^2))
#eeq = @. (3*exp(α)*(T^4))/(π^2)
plot!(τ,T, 
        xaxis   =:log ,
        #yaxis   =:log ,
        xlabel  =   "τ",  
        ylabel  =   "T",
        dpi     =   300)



