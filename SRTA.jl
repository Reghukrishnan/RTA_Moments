using Plots
using ProgressMeter
include("Auxiliary.jl")


#-----------SCALED INITIAL CONDITIONS---------------------

T₀  = 1             #  GeV     
m    = 0.1          # m corresponds to m = 0.1  Gev and T = 1GeV
α₀   = 0            # μ/T 


tₛ  =  0.1           # fm
tₑ  = 100
τ₀  = 0.1

η₀      =   (tₛ/τ₀)*(T₀/5)
#---------------------------------------------------
τᵣ⁰     =   (η₀/T₀)*(5)
ωᵣ⁰     =   (1/5)*(T₀/η₀)



#----------------
nₘₐₓ = 5 
nₐᵣ = [i for i in 0:nₘₐₓ]
nₙ = findall(x->x ==0,nₐᵣ)[1]   # Finding the location of 'n' for number desnity
nₑ = findall(x->x ==1,nₐᵣ)[1]   # Finding the lcoation of 'n' for energy density

N = size(nₐᵣ)[1]
L = 20



#------------------ Initial desnity-------------------------------------------------
χ₀ = Matrix{Float64}(undef, N, L+1)

Init_χ_Eq!(χ₀,nₐᵣ,L,m,T₀,α₀)   # Initialises the moments with equilibrium initial conditions. All values set to 1.

χ₀[nₙ,L+1] = α₀
χ₀[nₑ,L+1] = T₀


#------------------Printing the parameter values------------------------------

println("m     : ", m)
println("T₀    : ",T₀)
println("t₀    : ",tₛ)
println("η₀/s₀ : ",η₀)


Nₚ = 5000                        # Number of time steps
ProgressBar = Progress(Nₚ)
tspan = trange((tₛ,tₑ),Nₚ,"exp")  # Exponentially scaled time steps

#------------------------------------------------
println("\nSolving Scaled RTA.")
p = (N,L,ωᵣ⁰,nₐᵣ,nₙ,nₑ)

χ = RK4(χ₀,tspan,SRTA,p)     # Rk4 solver for differential equations. Returns the moments as a function of time χ[t,n,l]


#---------------------------------------------------
T = χ[:,nₑ,L+1]
τ = ((T).*(tspan)./((5)*η₀))  # \tau/\tau_R scaled time variable


plot(τ,χ[:,nₑ,1], 
        xaxis   =:log ,
        xlabel  =   "τ",  
        ylabel  =   "χ",
        ylims   =   (0.9975,1.0025),
        #yticks  =   [0,0.25,0.5,0.75,1],
        label   =   " - χ₁,₁",
        dpi     =   300)
plot(τ,χ[:,nₙ,1], 
        xaxis   =:log ,
        xlabel  =   "τ",  
        ylabel  =   "χ",
        ylims   =   (0.9975,1.0025),
        #yticks  =   [0,0.25,0.5,0.75,1],
        label   =   " - χ₀,₁",
        dpi     =   300)

plot(τ,T, 
        xaxis   =:log ,
        xlabel  =   "τ",  
        ylabel  =   "T",
        #ylims   =   (0,1),
        #yticks  =   [0,0.25,0.5,0.75,1],
        label   =   " T",
        dpi     =   300)