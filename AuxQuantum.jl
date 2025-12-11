using SpecialFunctions
using Integrals
using NLsolve


#Integrant 
using FastGaussQuadrature


Γ = gamma

uᵢ,uwᵢ  = gausslaguerre(100)
vᵢ,vwᵢ  = gausslegendre(50)




#---------------------------------------

function gnl(u,p)
    ζ = p[1]
    α = p[2]
    n = p[3]
    l = p[4]
    λ = p[5]
    r = p[6]
    y = ζ .+ u
    
    #println("λ : ", λ)
    #println("r : ", r)

    return @. ( (y^(n-2.0*l-2.0)) * (((y^2.0) - (ζ^2.0))^(l + 0.5))    )*( exp(-(λ-1)y + λ*α)/( ( 1 + r*exp(-y + α) )^λ) )  #*( exp(u)/(( exp(ζ*y - α) + r )^λ) )
end


function Gnl(n,l,ζ,α=0,r = 1,λ=1)
    #println("---------",ζ)
    if ζ == 0 && r==0 && λ ==1
        #println("Called - 1")
        return Γ(n)              
    else
        gx  = gnl(uᵢ,[ζ,α,n,l,λ,r])
        return exp(-ζ) *sum( (uwᵢ).*gx )
    end
end




function Init_ρ_Eq!(ρ₀,nₐᵣ,L,m,T,α,r)
    ζ = m/T
    for (n,nv) in enumerate(nₐᵣ)    
        for l in 1:L  
            #println(nv," ",l)      
            ρ₀[n,l] = ρeq(T,α,nv,l-1,r,1)  #+ ρeq(T₀,μ₀,nv,0)*(l/(2*l + 1)^(2))
        end
    end
end

function ρeq(T,α,n,l,r=0,λ=1)
    ζ = m/T
    return ((T^(n+3))/((2*l+1)*(2*π^2))) * Gnl(n+3,l,ζ,α,r,λ )
end



function Tα!(F,x)
    T = x[1]
    α = x[2]
    ζ = m/T                                 # m is defined globally
    F[1] = (T^4)*Gnl(1,0,ζ,α,r)/(2*π^2)      # r is defined globally
    F[2] = (T^3)*Gnl(0,0,ζ,α,r)/(2*π^2)
end





#---------------
function RTA(dρ::Matrix,ρ,t::Float64,p)
    N   = p[1]
    L   = p[2]
    ωᵣ⁰  = p[3]
    nₐᵣ = p[4]
    nₙ  = p[5]
    nₑ  = p[6]      
    γ  = p[7]       #γ = 0 implies free streaming system

    
    #T = (1/3)*(   ρ[nₑ,1]/ρ[nₙ,1] )
    #α = log((ρ[nₙ,1]*(π^2))/(T^3))

    T = ρ[nₑ,L+1]
    α = ρ[nₙ,L+1]
    
    # USe a root finding algorith to find the Temperature and α
    # Initial Guess, we use conformal MB results

    #Tₚ = (1/3)*(   ρ[nₑ,1]/ρ[nₙ,1] )
    #αₚ = log((ρ[nₙ,1]*(π^2))/(Tₚ^3))

    #sol = nlsolve(Tα!,[Tₚ,αₚ])
    #T, α = sol.zero
    
    ζ = m/T
    
    ωᵣ = ωᵣ⁰*((T/T₀))
    

    for (n,nᵥ) in enumerate(nₐᵣ)        
        #--------------------------------------------------------------------------
        for l = 1:L+1           
            #--------------------------------------------------------------------------

            # Here for lmax = L-1 we put the truncation condition.
            if l < L           

                free = (2*(l-1)+1)*ρ[n,l] + (nᵥ-2*(l-1))*ρ[n,l+1] # Free streaming part

                relx = γ*ωᵣ*( ρ[n,l] - ρeq(T,α,nᵥ,l-1,r,1) )        # Relaxation part

                dρ[n,l] = - ( free/t)  - ( relx )  

            elseif l == L

                free = (2(l-1)+1)*ρ[n,l] + (nᵥ-2(l-1))*ρeq(T,α,nᵥ,l,r,1) # L+1 moment is at equilibrium (closure condition)

                relx = γ*ωᵣ*( ρ[n,l] - ρeq(T,α,nᵥ,l-1,r,1) )

                dρ[n,l] = - ( free/t )  - ( relx )  
            else
                dρ[n,L+1] = 0  
            end         
        end
                
    end

    if r!= 0
        v = 1/T
        ϵd = ρ[nₑ,1]
        nd = ρ[nₙ,1]

        dTϵ = (3*v*ϵd) + ((v^2)*ρeq(T,α,2,0,r,1)) + (r*(v^2 )*ρeq(T,α,2,0,r,2)) 
        dTn = (3*v*nd) + ((v^2)*ρeq(T,α,1,0,r,1)) + (r*(v^2 )*ρeq(T,α,1,0,r,2))

        dαϵ = ϵd - r*ρeq(T,α,1,0,r,2) 
        dαn = nd - r*ρeq(T,α,0,0,r,2)

        D = dTϵ * dαn - dαϵ * dTn   # Determinant of the matrix

        dρ[nₑ,L+1] =  ((dαn * dρ[nₑ,1]) - (dαϵ * dρ[nₙ,1]) )/D

        dρ[nₙ,L+1] =  ((-dTn * dρ[nₑ,1]) + (dTϵ * dρ[nₙ,1]) )/D

    else
        
        χ11 = ρ[nₑ,2]/ρeq(T,α,1,1)

        dρ[nₑ,L+1] = ( G41*G30/Gd)*(χ11*ρ[nₑ,L+1])/(3t)

        dρ[nₙ,L+1] = -( (G41*G40)/Gd)*(χ11/(3t))- (1/t)
    end

    return dρ
end