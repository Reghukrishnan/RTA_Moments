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
end


function Gnl(n,l,ζ)
    #println("---------",ζ)
    if ζ == 0
        #println("Called - 1")
        return Γ(n)              
    else
        gx  = gnl(uᵢ,(ζ,n,l))
        return exp(-ζ) * sum( (uwᵢ).*gx )
    end
end

#----------------------------------------------------------------
function ρeq(T,α,n,l)
    ζ = m/T
    return (exp(α))*((T^(n+3))/((2l+1)*(2*(π^2)))) * Gnl(n+3,l,ζ)
end


function Init_ρ_Eq_b!(ρ₀,nₐᵣ,L,m,T,α)
    ζ = m/T
    for (n,nv) in enumerate(nₐᵣ)    
        for l in 1:L  
            #println(nv," ",l)      
            ρ₀[n,l] = ρeq(T,α,nv,l-1)  #+ ρeq(T₀,μ₀,nv,0)*(l/(2*l + 1)^(2))
        end
    end
end





function RTAB(dρ::Matrix,ρ,t::Float64,p)
    N   = p[1]
    L   = p[2]
    ωᵣ⁰  = p[3]
    nₐᵣ = p[4]
    nₙ  = p[5]
    nₑ  = p[6]      
    γ  = p[7]       #γ = 0 implies free streaming system

    if m==5
        T = (1/3)*(   ρ[nₑ,1]/ρ[nₙ,1] )
        α = log((ρ[nₙ,1]*(π^2))/(T^3))
    else
        T = ρ[nₑ,L+1]
        α = ρ[nₙ,L+1]
    end
    
    #T = ρ[nₑ,L+1]
    #α = ρ[nₙ,L+1]

    ζ = m/T
    
    # Use a root finding algorith to find the Temperature and α
    # Initial Guess, we use conformal MB results
    #sol = nsolve(Tα!,[Tₚ,αₚ])
    #T, α = sol.zero

    ωᵣ = (ωᵣ⁰)*(T/T₀)     # τᵣ⁰ = 5η₀/T₀. --> ωᵣ⁰ = 1/τᵣ⁰ = T₀/5η₀ --> ωᵣ = ωᵣ⁰(T/T₀)

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

                free = (  (2*(l-1)+1) *ρ[n,l] ) + ( (nᵥ-2*(l-1)) *ρ[n,l+1] )# Free streaming part

                relx = γ*ωᵣ*( ρ[n,l] - ρeq(T,α,nᵥ,l-1) )        # Relaxation part
                

                dρ[n,l] = - ( free/t)  - ( relx )  

            elseif l == L

                free = ((2*(l-1)+1)*ρ[n,l] ) + ((nᵥ-2*(l-1))*ρeq(T,α,nᵥ,l)) # L+1 moment is at equilibrium (closure condition)

                relx = γ*ωᵣ*( ρ[n,l] - ρeq(T,α,nᵥ,l-1) )

                dρ[n,l] = - ( free/t )  - ( relx )  
            else
                dρ[n,L+1] = 0  
            end         
        end
                
    end

    χ11 = ρ[nₑ,2]/ρeq(T,α,1,1)


    dρ[nₑ,L+1] = ( G41*G30/Gd)*(χ11*ρ[nₑ,L+1])/(3t)

    dρ[nₙ,L+1] = -( (G41*G40)/Gd)*(χ11/(3t))- (1/t)

    return dρ
end










