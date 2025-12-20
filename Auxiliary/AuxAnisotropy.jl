
using SpecialFunctions
using Integrals
using NLsolve


#Integrant 
using FastGaussQuadrature


Γ = gamma

uᵢ,uwᵢ  = gausslaguerre(100)
vᵢ,vwᵢ  = gausslegendre(50)


#--------------------Modified Maximum Entropy Distribution-------------------------------------------------------
function mnl(u,v,p)
     
    n = p[1] 
    l = p[2]
    ζ = p[3]
    β = p[4] 
    ξ = p[5]
    ϑ = p[6]
                            # v ∈ [-1,1], θ \in [0,π]
    y = ζ + u               # u ∈ [0,∞], y \in [ζ,∞]
    p = sqrt(y^2 - ζ^2)
    p² = (y^2 - ζ^2)

    
    s = (β*y + (p²/y)*( ϑ*(1 - (v^2)) + ξ*(v^2) )) # + δ*Pₗ⁴ 
    

    return  ( v^(2l) ) * ( y^(n-2l-2) ) * ( ( y^2 - ζ^2 )^(l + 0.5) ) * exp(u-s)
end



function Mnl(n,l,ζ,p)
    β,ξ,ϑ = p
    p = (n,l,ζ,β,ξ,ϑ)

    Int = 0.0
    for i in eachindex(uᵢ), j in eachindex(vᵢ)
        Int+= uwᵢ[i]*vwᵢ[j]*mnl(uᵢ[i],vᵢ[j],p)
    end
    return ((2l+1)*0.5 )*Int
end


#---------------------------- Romatshke Strickland Distribution-----------------------------------------
function rnl(u,v,p)
     
    n = p[1] 
    l = p[2]
    ζ = p[3]
    β = p[4]
    ξ = p[5]
    ϑ = p[6]
    

                        # v ∈ [-1,1], θ \in [0,π]
    y   = ζ + u             # u ∈ [0,∞], y \in [ζ,∞]
    p²  = (y^2 - ζ^2)
    return  ( v^(2l) ) * ( y^(n-2l-2) ) * ( ( y^2 - ζ^2 )^(l + 0.5) ) * exp(  u - ( sqrt( (β^2)*y^2 + p²*((ξ^2 )- 1)*(v^2) + ((ϑ^2) - 1)*ζ^2 )  )  ) 
end

function Rnl(n,l,ζ,p)
    β,ξ,ϑ = p 
    p = (n,l,ζ,β,ξ,ϑ)

    Int = 0.0
    for i in eachindex(uᵢ), j in eachindex(vᵢ)
        Int+= uwᵢ[i]*vwᵢ[j]*rnl(uᵢ[i],vᵢ[j],p)
    end
    return ((2l+1)*0.5)*Int
end
















#--------------------Initialises moments for an isotropic Maxwell Juttner equilibrium distibution --------------------------------
function InitχEq!(χ₀,nₐᵣ,L,m,T,α)

    for (n,nv) in enumerate(nₐᵣ)    
        for l in 0:L-1  
            #println(nv," ",l)      
            χ₀[n,l+1] = 1#+ ρeq(T₀,μ₀,nv,0)*(l/(2*l + 1)^(2))
        end
    end
end

#--------------------Initialises moments for an an-isotropic distibution --------------------------------
function InitχAIso!(χ₀,nₐᵣ,L,p,Idst = Mnl) # p = (ζ,α,Πₚ,πₚ)
    ζ   = p[1]
    πₚ  = p[2]      # Π/P ratio
    Πₚ  = p[3]      # π/P ratio

    β,ξ,ϑ = PAnIso(ζ,πₚ,Πₚ,Idst)   # Finding parameters for corresponding values of  Πₚ and πₚ

    println("Parameters Found \n")
    ΕnΠπP(ζ,[β,ξ,ϑ],Idst)         # Prints the thermodynamic parameters as a check

    p = (β,ξ,ϑ )
    println("\nInitialising scaled Moments\n")
    for (n,nv) in enumerate(nₐᵣ)    
        for l in 0:L-1  
            #println(nv," ",l)      
            χ₀[n,l+1] = ( Idst(n,l,ζ,p)/Gnl(n,l,ζ) )     # ρeq(T₀,μ₀,nv,0)*(l/(2*l + 1)^(2))
        end
    end

end



#------------------------------------------------------------------------------------------------
#                    PHYSICAL ANISOTROPY GENERATOR
#------------------------------------------------------------------------------------------------

function PhyπPζ(ζ)           # Returns a random compatible πₚ given a ζ
    ϵ       = Gnl(4,0,ζ)
    P₃      = Gnl(4,1,ζ) # 3P
    λₘₐₓ    = ϵ/P₃

    println("πₚ has to be between $(-2λₘₐₓ) and $(λₘₐₓ)")
    πₚ = -2λₘₐₓ + 3λₘₐₓ*rand()

    return πₚ
end

function PhyΠPζ(ζ)           # Returns a random compatible Πₚ given a ζ
    ϵ       = Gnl(4,0,ζ)
    P₃      = Gnl(4,1,ζ) # 3P
    λₘₐₓ    = ϵ/P₃

    println("Πₚ     = -1 + λ" )
    print("\n")
    if ζ == 0
        println("Mass is zero. System has no bulk pressure.")
        λ = 1  
    else 
        ϵ       = Gnl(4,0,ζ)
        P₃      = Gnl(4,1,ζ) # 3P
        λₘₐₓ    = ϵ/P₃

        println("λₘₐₓ   = ϵ/3P :",λₘₐₓ )
        print("\n")
        println("0 <= λ <= λₘₐₓ" )
        println("-1 <= Πₚ <= $(λₘₐₓ)" )
        λ = rand()*λₘₐₓ  
    end

    println("λ      : ",λ )
    Πₚ = -1 + λ
    println("Πₚ     : ",Πₚ )

    return Πₚ
end

function PhyΠπP(ζ)              # Returns random compatible πₚ and Πₚ given a ζ
    println("Πₚ     = -1 + λ" )
    print("\n")
    if ζ == 0
        println("Mass is zero. System has no bulk pressure.")
       λ = 1  
    else 
        ϵ   = Gnl(4,0,ζ)
        P₃   = Gnl(4,1,ζ) # 3P
        λₘₐₓ = ϵ/P₃
        println("λₘₐₓ   = ϵ/3P :",λₘₐₓ )
        print("\n")
        println("0 <= λ <= λₘₐₓ" )
        println("-1 <= Πₚ <= $(λₘₐₓ)" )
        println("$(-2λₘₐₓ) <= πₚ <= $(λₘₐₓ)")
        λ = rand()*λₘₐₓ  
    end

    println("λ      : ",λ )
    Πₚ = -1 + λ
    println("Πₚ     : ",Πₚ )

    print("\n")

    println("-2λ <= πₚ <= λ")
    πₚ = -2λ + 3λ*rand()
    println(" πₚ    : ",πₚ )
    return πₚ,Πₚ
end


function PhyΠPπ(ζ,Πₚ)           # Returns a random compatible πₚ given Πₚ and ζ

    if ζ == 0 && Πₚ != 0
        println("Πₚ cannot be non-zero for ζ = $(ζ)")
        println("Setting Πₚ = ",0)
        Πₚ = 0
    end

    println("Πₚ     = -1 + λ" )
    print("\n")

    λ = 1 + Πₚ
    println("λ      : ",λ )

    println("-2λ <= πₚ <= λ")
    πₚ = -2λ + 3λ*rand()
    println(" πₚ    : ",πₚ )
    return πₚ
end

function PhyπPΠ(ζ,πₚ)           # Returns a random compatible Πₚ  given πₚ and ζ

    # Check consistency 

    ϵ   = Gnl(4,0,ζ)
    P₃   = Gnl(4,1,ζ) # 3P
    λₘₐₓ = ϵ/P₃

    if !(-2λₘₐₓ <=  πₚ <= λₘₐₓ)
         throw(ArgumentError("πₚ has to be between $(-2λₘₐₓ) and $(λₘₐₓ)"))
    end

    if ζ == 0 
        println("Πₚ :",0 )
        return 0
    end
    

    println("Πₚ     = -1 + λ" )
    print("\n")
    println("0.5πₚ -1   <= Πₚ ")
    println("   πₚ -1   <= Πₚ ")
    print("\n")
    println("0.5πₚ  <= λ ")
    println("   πₚ  <= λ ")
    println("   0  <= λ ")

    λₘᵢₙ = max(πₚ,0)
    println("λₘᵢₙ = max(πₚ,0) :",λₘᵢₙ )

   
    println("λₘₐₓ   = ϵ/3P :",λₘₐₓ )


    λ = λₘᵢₙ + (λₘₐₓ - λₘᵢₙ)*rand()
    println("λ      : ",λ )
    Πₚ = -1 + λ

    println("Πₚ :",Πₚ )
    return πₚ
end


function ΠπTab(m)       # Table of shear and bulk pressures in Sunil Jaiswal's paper.  https://doi.org/10.1103/PhysRevC.105.024911

    πₚ = [-1.00 ,-1.00,-1.00, 0.99,-1.80, 0.00, 0.00] 
    Πₚ = [ 0.00 ,-0.25,-0.37, 0.00, 0.00,-0.25,-0.85]

    return [πₚ[m],Πₚ[m]]
end


#------------------------------------------------------------------------------------------------
#                    ANISOTROPY PARAMETER INITIALISATION
#------------------------------------------------------------------------------------------------
function PAnIso(ζ,πₚ,Πₚ,Idst=Mnl)       # Returns the parameter set of β,ξ,ϑ for given values of ζ,πₚ,Πₚ for a model initial distribution.
                                        # Idst is the function for initial distribution. In this case the default is Mnl ~ The Minimum entropy distribution.                                   

    if ζ == 0 && Πₚ != 0
        throw(ArgumentError("Πₚ cannot be non-zero for ζ = $(ζ)"))
    end

    ϵ   = Gnl(4,0,ζ)    # Scaled by T⁴/2π²
    P   = (Gnl(4,1,ζ)/3)


    PL  = (1 + Πₚ -     πₚ )*P
    PT  = (1 + Πₚ + 0.5*πₚ )*P

    function residuals!(res,x,p)
        β,ξ,ϑ = x
        p = (β,ξ,ϑ)
        
        
        ϵᵍ     =  Idst(4,0,ζ,p)
        PLᵍ    = (Idst(4,1,ζ,p)/3)
        PTᵍ    = ( ϵᵍ - PLᵍ - ((ζ^2)*Idst(2,0,ζ,p)) )*0.5

        res[1] =     (ϵᵍ/ϵ)   - 1.0
        res[2] =     (PLᵍ/PL) - 1.0
        res[3] =     (PTᵍ/PT) - 1.0 
    end
    
    
    
    prob = NonlinearProblem(residuals!,[1.0,0.0,0.0]; maxiters=2000)
    sol = solve(prob,RobustMultiNewton())   # RobustMultiNewton()

    println(" β     : ", sol.u[1] )
    println(" ξ     : ", sol.u[2] )
    println(" ϑ     : ", sol.u[3] )

    return sol.u
end



#------------------------------------------------------------------------------------------------
#                      COMPUTIG THE ANISOTROPIES
#------------------------------------------------------------------------------------------------

function ΠP(ζ=0.0001,p = [1.0,0.0,0.0],Idst=Mnl)    # Returns the scaled bulk pressure Π/P  for a given distribution with the p[arameter set β,ξ,ϑ

    return (  ( Idst(4,0,ζ,p ) -(ζ^2)*Idst(2,0,ζ,p ) )/Gnl(4,1,ζ)  ) -1

end


function πP(ζ=0.0001,p = [1.0,0.0,0.0],Idst=Mnl)    # Returns the scaled shear pressure π/P  for a given distribution with the p[arameter set β,ξ,ϑ

    return 1 + ΠP(ζ,p,Idst)  - ( Idst(4,1,ζ,p )/Gnl(4,1,ζ))

end 
#---------------------------------------------------------------------
function ΕnΠπP(ζ=0.0001,p=[1.0,0.0,0.0],Idst=Mnl)   # Returns ϵ,n,πₚ,Πₚ  for a given distribution with the p[arameter set β,ξ,ϑ
                                                    #  p = (β,ξ,ϑ) of Modified maximum entropy distribution.

    ϵ   = (Idst(4,0,ζ,p)/Gnl(4,0,ζ))
    n   = (Idst(3,0,ζ,p)/Gnl(3,0,ζ))
    
    πₚ  =   πP(ζ,p,Idst)
    Πₚ  =   ΠP(ζ,p,Idst)

    println(" ϵ     : ", ϵ )
    println(" n     : ", n )
    println(" πₚ    : ", πₚ )
    println(" Πₚ    : ", Πₚ )

    return ϵ,n,πₚ,Πₚ 
end