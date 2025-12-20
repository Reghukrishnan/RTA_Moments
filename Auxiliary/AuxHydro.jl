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



function H_IS!(du, u, τ, p)
    # Unpack variables
    ε, n, Π, π = u

    # Unpack parameters
    cₛ², η, ζ, τπ, τΠ = p

    # Equation of state
    P = cₛ² * ε

    θ = 1 / τ   # Bjorken expansion rate

    # Energy density
    du[1] = -(ε + P + Π - π) / τ

    # Number density
    du[2] = -n / τ

    # Bulk pressure
    du[3] = -(Π + ζ * θ) / τΠ

    # Shear stress
    du[4] = -(π - (4η/3) * θ) / τπ
end



#--------------------------------------------------
# DNMR Bjorken hydro (conformal, with number density)
# u = [ε, n, π]
#--------------------------------------------------
function DNMR_BC!(du, u, τ, p)

    ε, n, π = u
    η = p[1]

    # EOS
    P = ε / 3
    θ = 1 / τ

    # Transport coefficients (RTA)
    τπ = 5η / (ε + P)
    λ1 = (5/7) * η * τπ

    # Energy density
    du[1] = -(ε + P - π) * θ

    # Number density
    du[2] = -n * θ

    # Shear stress (DNMR)
    du[3] = -π / τπ
            + (4η / (3τπ)) * θ
            - (4/3) * π * θ
            - (((λ1 / (2η^2)) * π^2) / τπ)

    return nothing
end






#--------------------------------------------------
# DNMR Bjorken hydro (non-conformal)
# u = [ε, n, Π, π]
#--------------------------------------------------

function EoS(ε::Real, n::Real)

end



function DNMR_BNC!(du, u, τ, p)

    T, α, Π, πₛ = u
    m, η, ζ = p

    z = m/T
    # EOS: user-supplied P(ε,n)
    ε = exp(α)*(T^4)*Gnl(1,0,z)/(2π^2)
    P = exp(α)*(T^4)*Gnl(1,1,z)/(6π^2)

    #P = EoS(ε, n)   
    #θ = 1 / τ

    # Transport coefficients
    βΠ = ((1/3) - cₛ²)(ε + P) - (2/9)(ε - 3P)
    τπ = 5η / (ε + P)
    τΠ = ζ / βΠ

    λ1 = (5/7) * η * τπ

    δππ = 4/3
    δΠΠ = 1.0
    λπΠ = 6/5
    λΠπ = 8/5

    # Temperature ~~Energy density
    du[1] = ( G41*G30/Gd)*(χ11*u[1])/(3t)    #du[1] = -(ε + P + Π - π) * θ

    # Checmical potential ~~ Number density
    du[2]= -( (G41*G40)/Gd)*(χ11/(3t))- (1/t)      #du[2] = -n * θ

    # Bulk pressure
    du[3] = -Π / τΠ
            - ζ * θ / τΠ
            - δΠΠ * Π * θ
            + λΠπ * πₛ * θ

    # Shear stress
    du[4] = -πₛ / τπ
            + (4η / (3τπ)) * θ
            - δππ * πₛ * θ
            + λπΠ * Π * θ
            - (λ1 / (2η^2)) * πₛ^2 / τπ

    return nothing
end