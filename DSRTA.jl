using Plots
using ProgressMeter
include("Auxiliary.jl")



# In this code we use scaled variables Χnl for the moments.
# u = τ/τ₀ 
# 


fm_to_Gev = 1.973
Gev_to_fm = 5.068
#-----------INITIAL CONDITIONS---------------------------------

m       =  0.4          # GeV           -  m corresponds to m = 0.1  Gev and T = 1GeV
α₀      =  0.0          # dimentionless -  μ/T 

T₀      =  1.0          # GeV           -  Temperature
η₀      =  0.02         # dimentionless -  η/s specific viscosity 




#-----------SCALED INITIAL CONDITIONS--------------------------
ζ₀      =  m/T₀          # dimentionless proxy for mass
uᵣ⁰     =  10            # dimentionless relaxation strength - 5*η₀/τ₀T₀ = τ_R⁰/τ_0
w₀      =  1/uᵣ⁰         # dimentionless - τ₀/τ_R 

println("\n----------------------------------\n")

println("m      : ", m)
println("T₀     : ", T₀)
println("ζ₀     : ", ζ₀)


ϵ₀ = Gnl(4,0,ζ₀)
P  = Gnl(4,1,ζ₀)/3

println("-1 ≤ Πₚ ≤ $(-1 + (ϵ₀/(3P)))")

#----------INITIAL ANISOTROPIES----------------------------------
                        # PhyΠπP(ζ)     # Returns random compatible πₚ and Πₚ given a ζ 
Πₚ      = 0.0           # PhyπPΠ(ζ,πₚ)   # Returns a random compatible Πₚ given πₚ and ζ           
πₚ      = 0.0           # PhyΠPπ(ζ,Πₚ)   # Returns a random compatible πₚ given Πₚ and ζ

      
πₚ ,Πₚ = PhyΠπP(ζ₀)


#----------------------------------------------------------------


h₀      =  0.0          # dimentionless -  h = ln(T/T₀)
u₀      =  1.0          # dimentionless -  u = τ/τ₀
γ₀      =  0.0          # dimentionless - proxy for time, decay strength ∫dτ/τᵣ = ∫du/uᵣ
#---------------------------------------------------------------


#------------------Printing the parameter values------------------------------

println("uᵣ⁰    : ", uᵣ⁰, " (relaxation strength)")

println("\n----------------------------------\n")

println("h₀     : ", h₀, " ( log(T/T₀) ) ")
println("u₀     : ", u₀)

println("\n----------------------------------\n")
#-----------Initialising the moment array-------------------------
nₘᵢₙ    = -1
nₘₐₓ    = 5 
nₐᵣ     = [i for i in nₘᵢₙ:nₘₐₓ]
nₙ      = findall(x->x ==0,nₐᵣ)[1]      # Finding the location of 'n' for number desnity
nₑ      = findall(x->x ==1,nₐᵣ)[1]      # Finding the lcoation of 'n' for energy density

N       = size(nₐᵣ)[1]
L       = 20



#------------------ Initialising the scaled moments-------------------------------------------------
χ₀      = Matrix{Float64}(undef, N, L+1)

p = (ζ₀,πₚ,Πₚ)

InitχAIso!(χ₀,nₐᵣ,L,p)  

χ₀[nₙ,L+1]      = α₀           # dimentionless initial chemical potential
χ₀[nₑ,L+1]      = h₀           # dimentionless initial proxy for temperature
χ₀[nₑ+1,L+1]    = u₀           # dimentionless initial proxy for time u = τ/τ₀ -> u₀ = 1.0



Nₚ = 8000                       # Number of time steps
ProgressBar = Progress(Nₚ)

γ = trange((γ₀,10),Nₚ,"exp",10)  # Exponentially scaled time steps

#------------------------------------------------
println("\n Solving Scaled RTA.")
p = (nₐᵣ,nₙ,nₑ,ζ₀,uᵣ⁰)

χ = RK4(χ₀,γ,DSRTA,p)           # Rk4 solver for differential equations. Returns the moments as a function of time χ[t,n,l]


println("\nSolving Completed.")
#---------------------------------------------------
T  = exp.(χ[:,nₑ,L+1])          # exp(h) -> h = log(T/T₀ )
α  = χ[:,nₙ,L+1] 


uᵣ = uᵣ⁰./T                     # uᵣ = τᵣ/τ₀
u  = χ[:,nₑ+1,L+1]              # u  = τ/τ₀    
       
w  = u./uᵣ                      # τ = τ/τᵣ = u/uᵣ



r = (w.-w₀)./γ

plot!( w,T,
        xaxis   =:log,
        #xticks  =   [0,1,2,3,4,5,6,7,8,9,10], 
        #yticks  =   [0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0],
        label   = "T ",
        xlabel  = "w",
        ylabel  = "T")


