using SpecialFunctions
using Integrals
using NLsolve




#Integrant 
using FastGaussQuadrature


Γ = gamma

uᵢ,uwᵢ  = gausslaguerre(100)
vᵢ,vwᵢ  = gausslegendre(50)
#Integrant 
function gnl(u,p)
    ζ = p[1]
    n = p[2]
    l = p[3]
    y = ζ .+ u

    return @. ( y^(n-2.0*l-2.0))*(y^2.0 - ζ^2.0)^(l + 0.5)
    y = ζ .+ u

    return @. ( y^(n-2.0*l-2.0))*(y^2.0 - ζ^2.0)^(l + 0.5)
end

function Gnl(n,l,ζ)
    #println("---------",ζ)
function Gnl(n,l,ζ)
    #println("---------",ζ)
    if ζ == 0
        #println("Called - 1")
        return Γ(n)        
        #println("Called - 1")
        return Γ(n)        
    else
        gx  = gnl(uᵢ,(ζ,n,l))
        return exp(-ζ) * sum( (uwᵢ).*gx )
        gx  = gnl(uᵢ,(ζ,n,l))
        return exp(-ζ) * sum( (uwᵢ).*gx )
    end
end




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




#--------------------Initialises moments for an isotropic Maxwell Juttner equilibrium distibution --------------------------------
function InitχEq!(χ₀,nₐᵣ,L,m,T,α,Πₚ,πₚ)

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




#---------------------------------------------------------------------
function trange(tspan,N,interpolation,scale= 1.0)
    lin = range(0,stop=1,length=N)
    tₛ = tspan[1]
    tₑ = tspan[2]
    Δt = tₑ - tₛ 

    if interpolation == "exp" # Exponential spacing. Dense near first point.
        
        t = tₛ .+ ((exp.(scale*lin).-1) .*(Δt/(exp(scale*1.0) -1.0)))
        return t

    elseif interpolation == "lin"  # linear spacing. Equally spaced at all points
        tₛ       = tspan[1]
        tₑ      = tspan[2]
        step    = (tₑ -  tₛ)/(N)
        return [tₛ + (i)*step for i in 0:N-1 ]
    elseif interpolation == "pol"
        γ = (1)/n
        tₛ       = (tspan[1])^γ
        tₑ      = (tspan[2])^γ
        step    = (tₑ -  tₛ)/(N)
        return [(tₛ + (i)*step)^n for i in 0:N-1 ]
    elseif interpolation == "cos"  # cosine spacing. Densly spaced at end points
        strength = 3
        t = 0.5*(1 .- cos.(π*lin)).^strength
        t .= tₛ .+ t.*Δt
    else
        println("Invalid Interpolation input.") 
    end    

    
end


#-------------------------------------------------------------------------------------------------------

function SRTA(dχ::Matrix,χ,t::Float64,p)
    N   = p[1]
    L   = p[2]
    ωᵣ⁰   = p[3]
    nₐᵣ = p[4]
    nₙ  = p[5]
    nₑ  = p[6]
    γ  = p[7]

    α = χ[nₙ,L+1]
    T = χ[nₑ,L+1]
    #println(χ[nₑ,1])
    
    
    #println(m, "\n\n\n\n\n")
    if m == 0
        ζ = 0
    else
        ζ = m/T
    end
    #println("------------",T)
    
    
    
    ωᵣ = ωᵣ⁰*((T/T₀))
    G30 = Gnl(3.0,0,ζ)
    G40 = Gnl(4.0,0,ζ)
    G41 = Gnl(4.0,1.0,ζ)
    G50 = Gnl(5.0,0,ζ)
    
    Gd = G40^2 - G30*G50
    

    for (n,nᵥ) in enumerate(nₐᵣ)        
        #--------------------------------------------------------------------------
        for l = 1:L+1           
            #--------------------------------------------------------------------------

            # Here for lmax = L-1 we put the truncation condition.
            if l < L           
                h = (1/3)*( (( Gnl(ζ,nᵥ + 4.0,l-1)*G30 -  G40*Gnl(ζ,nᵥ + 3.0,l-1) )/(Gd))*(G41/Gnl(ζ,nᵥ + 3.0,l-1))*(χ[nₑ,2]* χ[n,l])  +  3*(nᵥ-2(l-1))*((2*(l-1) +1.0)/(2*(l-1) +3.0))*( Gnl(ζ,nᵥ + 3.0,l)/Gnl(ζ,nᵥ + 3.0,l-1) )* χ[n,l+1] )
                
                f = (2*(l-1))*χ[n,l] 
                r = γ*ωᵣ*( χ[n,l] - 1 )

                dχ[n,l] = -( h/t)- ( f/t)  - ( r )  
            elseif l == L
                h = (1/3)*( (( Gnl(ζ,nᵥ + 4,l-1)*G30 -  G40*Gnl(ζ,nᵥ + 3,l-1) )/(Gd))*(G41/Gnl(ζ,nᵥ + 3,l-1))*(χ[nₑ,2]* χ[n,l])  +  3* (nᵥ-2(l-1)) *( (2*(l-1) +1)/(2*(l-1) +3) )*( Gnl(ζ,nᵥ + 3.0,l)/Gnl(ζ,nᵥ + 3.0,l-1) )* 1 )
                
                f = (2*(l-1))*χ[n,l] 
                r = γ*ωᵣ*( χ[n,l] - 1 )

                dχ[n,l] = -( h/t)- ( f/t)  - ( r )  
            else
                dχ[n,L+1] = 0  
            end         
        end
                
    end
    dχ[nₑ,L+1] = ( G41*G30/Gd)*χ[nₑ,2]*χ[nₑ,L+1]/(3t)

    dχ[nₙ,L+1] = -( (G41*G40)/Gd)*(χ[nₑ,2]/(3t))- (1/t)
    return dχ
end

#-------------------------------------------------------------------------------------------------------

function DSRTA(dχ::Matrix,χ,t::Float64,p)

    # p = (N,L,η₀,κ₀,nₐᵣ,nₙ,nₑ)
    nₐᵣ = p[1]
    nₙ  = p[2]
    nₑ  = p[3]
    ζ₀  = p[4]
    uᵣ⁰ = p[5]
  

    α = χ[nₙ,L+1]       # μ/T
    h = χ[nₑ,L+1]       # ln(τ₀T)
    u = χ[nₑ+1,L+1]     # τ/τ₀
    

    T   = exp(h) 
    uᵣ  = uᵣ⁰/T         # τᵣ/τ_₀ = η₀/(τ_₀T)
    #print(T)
    if ζ₀ == 0
        ζ   = 0
    else
        ζ   = ζ₀/T      # m/T, (m/T₀)*(T₀/T) = 
    end
    

    w = uᵣ/u
    

    G30 = Gnl(3.0,0,ζ)
    G40 = Gnl(4.0,0,ζ)
    G41 = Gnl(4.0,1.0,ζ)
    G50 = Gnl(5.0,0,ζ)
    
    
    Gd = G40^2 - G30*G50
    

    for (n,nᵥ) in enumerate(nₐᵣ)        
        #--------------------------------------------------------------------------
        for l = 1:L+1           
            #--------------------------------------------------------------------------

            # Here for lmax = L-1 we put the truncation condition.
            if l < L           
                g = (1/3)*( (( Gnl(ζ,nᵥ + 4.0,l-1)*G30 -  G40*Gnl(ζ,nᵥ + 3.0,l-1) )/(Gd))*(G41/Gnl(ζ,nᵥ + 3.0,l-1))*(χ[nₑ,2]* χ[n,l])  +  3*(nᵥ-2(l-1))*((2*(l-1) +1.0)/(2*(l-1) +3.0))*( Gnl(ζ,nᵥ + 3.0,l)/Gnl(ζ,nᵥ + 3.0,l-1) )* χ[n,l+1] )
                
                f = (2*(l-1))*χ[n,l] 
                r = ( χ[n,l] - 1 )

                dχ[n,l] = -( g +  f)*(w)  - ( r )  
            elseif l == L
                g = (1/3)*( (( Gnl(ζ,nᵥ + 4,l-1)*G30 -  G40*Gnl(ζ,nᵥ + 3,l-1) )/(Gd))*(G41/Gnl(ζ,nᵥ + 3,l-1))*(χ[nₑ,2]* χ[n,l])  +  3* (nᵥ-2(l-1)) *( (2*(l-1) +1)/(2*(l-1) +3) )*( Gnl(ζ,nᵥ + 3.0,l)/Gnl(ζ,nᵥ + 3.0,l-1) )* 1 )
                
                f = (2*(l-1))*χ[n,l] 
                r = ( χ[n,l] - 1 )

                dχ[n,l] = -( g +  f)*(w)  - ( r )   
            else
                dχ[n,L+1] = 0  
            end         
        end
                
    end
    dχ[nₑ,L+1] = ( G41*G30/Gd)*χ[nₑ,2]*(w/3)  # h ~ τ₀T
    
    dχ[nₙ,L+1] = -( (G41*G40)/Gd)*(χ[nₑ,2]*(w/3))- (w)

    #println("dχ[nₑ,L+1]     :",dχ[nₑ,L+1] )
    dχ[nₑ+1,L+1] = uᵣ

    return dχ
end


#---------------------------------------







function g_q(y,p)
    ζ = p[1]
    α = p[2]
    n = p[3]
    l = p[4]
    λ = p[5]
    r = p[5]
    return ( y^(n-2.0*l-2))*(y^2.0 - 1)^(l + 0.5)/( (exp(ζ*y - α) + r )^λ )
end

function Gnl_q(ζ,α,n,l,λ=1,r = 1)       # r = 0 MB, r = 1 FD, r = -1 BE
    
    if ζ == 0
        return Γ(n)
    else
        domain  = (1.0,ζ*10 + 1000.0)
        
        prob    = IntegralProblem(g, domain,[ζ,α,n,l,λ,r] )
        Gₙₗ      = solve(prob, QuadGKJL())
        return (ζ^(n+1))*Gₙₗ[1]
    end
end

function Init_χ_Eq_q!(ρ₀,nₐᵣ,L,m,T,α,r)
    ζ = m/T
    for (n,nv) in enumerate(nₐᵣ)    
        for l in 1:L  
            #println(nv," ",l)      
            ρ₀[n,l] = ρeq_q(T,α,n,l,r)  #+ ρeq(T₀,μ₀,nv,0)*(l/(2*l + 1)^(2))
        end
    end
end

function ρeq_q(T,α,n,l,r=0,λ=1)
    ζ = m/T
    return ((T^(n+3))/(2*(l-1)+1)*(2π)^2)*Gnl_q(ζ,α,n,l-1,r )
end




function Tα!(F,x)
    T = x[1]
    α = x[2]
    ζ = m/T                                 # m is defined globally
    F[1] = (T^4)*Gnl_q(ζ,α,1,0,r)/(2π)      # r is defined globally
    F[2] = (T^3)*Gnl_q(ζ,α,0,0,r)/(2π)
end





#---------------
function RTA(dρ::Matrix,ρ,t::Float64,p)
    N   = p[1]
    L   = p[2]
    ωᵣ⁰   = p[3]
    nₐᵣ = p[4]
    nₙ  = p[5]
    nₑ  = p[6]      
    γ  = p[7]       #γ = 0 implies free streaming system

    ϵd  =  ρ[nₑ,1]
    nd  =  ρ[nₙ,1]
    
    
    T = ρ[nₑ,L+1]
    α = ρ[nₙ,L+1]
    
    # USe a root finding algorith to find the Temperature and α
    # Initial Guess, we use conformal MB results
    #sol = nsolve(Tα!,[Tₚ,αₚ])
    #T, α = sol.zero
    
    ζ = m/T
    
    ωᵣ = ωᵣ⁰*((T/T₀))
    

    for (n,nᵥ) in enumerate(nₐᵣ)        
        #--------------------------------------------------------------------------
        for l = 1:L+1           
            #--------------------------------------------------------------------------

            # Here for lmax = L-1 we put the truncation condition.
            if l < L           

                free = (2(l-1)+1)*ρ[n,l] + (n-2l)*ρ[n,l+1] # Free streaming part

                relx = γ*ωᵣ*( ρ[n,l] - ρeq_q(T,α,n,l-1,r,1) )        # Relaxation part

                dρ[n,l] = - ( free/t)  - ( relx )  

            elseif l == L

                free = (2(l-1)+1)*ρ[n,l] + (n-2l)*ρeq_q(T,α,n,l-1,r,1) # L+1 moment is at equilibrium (closure condition)

                relx = γ*ωᵣ*( ρ[n,l] - ρeq_q(T,α,n,l-1,r,1) )

                dρ[n,l] = - ( free/t)  - ( relx )  
            else
                dρ[n,L+1] = 0  
            end         
        end
                
    end

    v = 1/T
    dTϵ = 3v*ϵd + (v^2)*ρeq_q(T,α,2,0,r,1) + r*(v^2 )*ρeq_q(T,α,2,0,r,2) 
    dTn = 3v*nd + (v^2)*ρeq_q(T,α,1,0,r,1) + r*(v^2 )*ρeq_q(T,α,1,0,r,2)

    dαϵ = -ϵd + r*ρeq_q(T,α,1,0,r,2) 
    dαn = -nd + r*ρeq_q(T,α,1,0,r,2)

    D = dTϵ * dαn - dαϵ * dTn   # Determinant of the matrix

    dρ[nₑ,L+1] =  (dαn * dρ[nₑ,1] - dαϵ * dρ[nₙ,1] )/D

    dρ[nₙ,L+1] =  (-dTn * dρ[nₑ,1] + dTϵ * dρ[nₙ,1] )/D

    return dρ
end


#-------------------------------------------------------------------------------------------------------
function RK4(y₀::Matrix,ts::Vector,fun,p = nothing)
    print("Array RK4 Solver Called. \n")
    
    ysize = (size(ts)...,size(y₀)...)  # Get the size of the time series amtrix
    y  = Array{Float64}(undef,ysize )  # Create the time series matrix
    dy = similar(y₀)                   # Create matrix to store dy
    view(y,1,:,:) .= y₀                # Initialise t = tᵢ value of y to y₀
    
    # Checks if the function takes any extra parameters.
    if isnothing(p) 
        for (i,t) in enumerate(ts[begin:end-1])

            dt = ts[i+1] - ts[i]
            
            yᵢ = view(y,i,:,:)              # Creates a view of the Matrix of variables

            # Computes the RK4 parameters
            k₁ = fun(dy, yᵢ,t)
            k₂ = fun(dy, yᵢ .+ 0.5*dt*k₁, t + 0.5*dt)
            k₃ = fun(dy, yᵢ .+ 0.5*dt*k₂, t + 0.5*dt)
            k₄ = fun(dy, yᵢ .+     dt*k₃, t +     dt)

            view(y,i+1,:,:) .= yᵢ .+ (dt*( k₁ .+ 2k₂ .+ 2k₃ .+ k₄ )/6)
            next!(ProgressBar)  
        end

    else
        for (i,t) in enumerate(ts[begin:end-1])

            dt = ts[i+1] - ts[i]
            
            yᵢ = view(y,i,:,:)             # Creates a view of the Matrix of variables

            # Computes the RK4 parameters
            k₁ = fun(dy, yᵢ,t,p)
            k₂ = fun(dy, yᵢ .+ 0.5*dt*k₁, t + 0.5*dt,p)
            k₃ = fun(dy, yᵢ .+ 0.5*dt*k₂, t + 0.5*dt,p)
            k₄ = fun(dy, yᵢ .+     dt*k₃, t +     dt,p)

            view(y,i+1,:,:) .= yᵢ .+ (dt*( k₁ .+ 2k₂ .+ 2k₃ .+ k₄ )/6)
            next!(ProgressBar)  
        end
    end
    return y
    
end