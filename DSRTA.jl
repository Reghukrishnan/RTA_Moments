using Plots
using ProgressMeter
include("Auxiliary.jl")



# In this code we use scaled variables Χnl for the moments.
# u = τ/τ₀ 
# 


fm_to_Gev = 1.973
Gev_to_fm = 5.068
#-----------INITIAL CONDITIONS---------------------------------

m       =  0.0          # GeV           -  m corresponds to m = 0.1  Gev and T = 1GeV
α₀      =  0.0          # dimentionless -  μ/T 

T₀      =  1.0          # GeV           -  Temperature
η₀      =  0.02         # dimentionless -  η/s specific viscosity 
 

#-----------SCALED INITIAL CONDITIONS--------------------------
ζ₀      =  m/T₀          # dimentionless proxy for mass
uᵣ⁰     =  10            # dimentionless relaxation strength - 5*η₀/τ₀T₀ = τ_R⁰/τ_0
w₀      =  1/uᵣ⁰         # dimentionless - τ₀/τ_R 

#--------------------------------------------------------------

h₀      =  0            # dimentionless -  h = ln(T/T₀)
u₀      =  1.0          # dimentionless -  u = τ/τ₀
ξ₀      =  0            # dimentionless - proxy for time, decay strength ∫dτ/τᵣ = ∫du/uᵣ
#---------------------------------------------------------------



#-----------Initialising the moment array-------------------------
nₘₐₓ    = 5 
nₐᵣ     = [i for i in 0:nₘₐₓ]
nₙ      = findall(x->x ==0,nₐᵣ)[1]   # Finding the location of 'n' for number desnity
nₑ      = findall(x->x ==1,nₐᵣ)[1]   # Finding the lcoation of 'n' for energy density

N       = size(nₐᵣ)[1]
L       = 20



#------------------ Initialising the scaled moments-------------------------------------------------
χ₀      = Matrix{Float64}(undef, N, L+1)

# Init_χ_Eq!(χ₀,nₐᵣ,L,m,T₀,α₀)   # Initialises the moments with equilibrium initial conditions. All values set to 1.


p = (m,T,α,Πₚ,πₚ)

Init_χ_AIso!(χ₀,nₐᵣ,L,p)  

χ₀[nₙ,L+1]      = α₀           # dimentionless initial chemical potential
χ₀[nₑ,L+1]      = h₀           # dimentionless initial proxy for temperature
χ₀[nₑ+1,L+1]    = u₀           # dimentionless initial proxy for time u = τ/τ₀ -> u₀ = 1.0


#------------------Printing the parameter values------------------------------

println("m      : ", m)
println("T₀     : ", T₀)
println("ζ₀     : ", ζ₀)

println("uᵣ⁰    : ", uᵣ⁰, " (relaxation strength)")

println("------------------------")
println("h₀     : ", h₀, " ( log(T/T₀) ) ")
println("u₀     : ", u₀)



Nₚ = 5000                       # Number of time steps
ProgressBar = Progress(Nₚ)

ξ = trange((ξ₀,10),Nₚ,"exp",1)  # Exponentially scaled time steps

#------------------------------------------------
println("\n Solving Scaled RTA.")
p = (nₐᵣ,nₙ,nₑ,ζ₀,uᵣ⁰)

χ = RK4(χ₀,ξ,DSRTA,p)           # Rk4 solver for differential equations. Returns the moments as a function of time χ[t,n,l]

println("\nSolving Completed.")
#---------------------------------------------------
T  = exp.(χ[:,nₑ,L+1])          # exp(h)/τ₀ -> h = log(τ₀T)
α  = χ[:,nₙ,L+1] 


uᵣ = uᵣ⁰./T                     # uᵣ = τᵣ/τ₀
u  = χ[:,nₑ+1,L+1]              # u  = τ/τ₀    
       
w  = u./uᵣ                      # τ = τ/τᵣ = u/uᵣ



r = (w.-w₀)./ξ

plot( ξ,χ[:,nₑ,2],
        xticks  =   [0,1,2,3,4,5,6,7,8,9,10], 
        #yticks  =   [0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0],
        label   = "$(ζ₀),$(uᵣ⁰) ",
        xlabel  = "ξ",
        ylabel  = "(w.-w₀)/ξ")


