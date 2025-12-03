using Plots
using ProgressMeter
include("Auxiliary.jl")


#-----------SCALED INITIAL CONDITIONS---------------------

T₀  = 1             #  GeV     
m    = 0.1          # m corresponds to m = 0.1  Gev and T = 1GeV
α₀   = 0            # μ/T 

r    = 0            # MB =0, FD = 1, BE = -1


tₛ  =  0.1           # fm
tₑ  = 100
τ₀  = 0.1

η₀      =   (tₛ/τ₀)*(T₀/5)
#---------------------------------------------------
τᵣ⁰     =   (η₀/T₀)*(5)
ωᵣ⁰     =   (1/5)*(T₀/η₀)



#----------------
nₘₐₓ = 5 
nₐᵣ = [i for i in -1:nₘₐₓ]
nₙ = findall(x->x ==0,nₐᵣ)[1]   # Finding the location of 'n' for number desnity
nₑ = findall(x->x ==1,nₐᵣ)[1]   # Finding the lcoation of 'n' for energy density

N = size(nₐᵣ)[1]
L = 20

#------------------ Initial desnity-------------------------------------------------
ρ₀ = Matrix{Float64}(undef, N, L+1)

Init_χ_Eq_q!(ρ₀,nₐᵣ,L,m,T₀,α₀,r)   # Initialises the moments with equilibrium initial conditions.

ρ₀[nₙ,L+1] = α₀
ρ₀[nₑ,L+1] = T₀

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
println("η₀/s₀ : ",η₀)


Nₚ = 5000                        # Number of time steps
ProgressBar = Progress(Nₚ)
tspan   = trange((tₛ,tₑ),Nₚ,"exp")  # Exponentially scaled time steps
T       = zeros(Float64, Nₚ)
α       = zeros(Float64, Nₚ)

#------------------------------------------------
println("\nSolving Scaled RTA.")
p = (N,L,ωᵣ⁰,nₐᵣ,nₙ,nₑ,1)

ρ = RK4(ρ₀,tspan,RTA,p)     # Rk4 solver for differential equations. Returns the moments as a function of time χ[t,n,l]


#---------------------------------------------------
T = ρ[:,nₑ,L+1]
τ = ((T).*(tspan)./((5)*η₀))  # \tau/\tau_R scaled time variable


plot(τ,T, 
        xaxis   =:log ,
        xlabel  =   "τ",  
        ylabel  =   "T",
        dpi     =   300)