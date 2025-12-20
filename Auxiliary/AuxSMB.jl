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
                h = (1/3)*( (( Gnl(nᵥ + 4.0,l-1,ζ )*G30 -  G40*Gnl(nᵥ + 3.0,l-1,ζ) )/(Gd))*(G41/Gnl(nᵥ + 3.0,l-1,ζ))*(χ[nₑ,2]* χ[n,l])  +  3*(nᵥ-2(l-1))*((2*(l-1) +1.0)/(2*(l-1) +3.0))*( Gnl(nᵥ + 3.0,l,ζ)/Gnl(nᵥ + 3.0,l-1,ζ) )* χ[n,l+1] )
                
                f = (2*(l-1))*χ[n,l] 
                r = γ*ωᵣ*( χ[n,l] - 1 )

                dχ[n,l] = -( h/t)- ( f/t)  - ( r )  
            elseif l == L
                h = (1/3)*( (( Gnl(nᵥ + 4,l-1,ζ)*G30 -  G40*Gnl(nᵥ + 3,l-1,ζ) )/(Gd))*(G41/Gnl(nᵥ + 3,l-1,ζ))*(χ[nₑ,2]* χ[n,l])  +  3* (nᵥ-2(l-1)) *( (2*(l-1) +1)/(2*(l-1) +3) )*( Gnl(nᵥ + 3.0,l,ζ)/Gnl(nᵥ + 3.0,l-1,ζ) )* 1 )
                
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