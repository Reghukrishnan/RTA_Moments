using SpecialFunctions 
using AbstractAlgebra

Γ = gamma
Pr = rising_factorial

"""
    K_nlk(n,l,k,τ,τ0)

Kernel appearing in the exact free-streaming solution.
"""
function K_nlk(n::Real, l::Int, k::Int, τ::Real, τ0::Real)
    prefactor = (τ0/τ)^(2l + 1)
    x = 1 - (τ0/τ)^2
    return Pr(l - n/2, k) / Γ(k+1) * prefactor * x^k
end



"""
    Bjorken_FS(n,l,τ,τ0,g0; kmax=50, tol=1e-12)

Compute g_{n,l}(τ) from initial data g_{n,l+k}(τ0).
"""
function Bjorken_FS(
    n::Real,
    l::Int,
    τ::Real,
    τ0::Real,
    g0::Function;
    kmax::Int = 50,
    tol::Real = 1e-12
)
    s = 0.0
    for k in 0:kmax
        term = K_nlk(n,l,k,τ,τ0) * g0(l + k)
        s += term
        abs(term) < tol && break
    end
    return s
end



#-----Time Evolution of the n,l moment for free streaming.

function Analytical_BFS!(    ρ::Vector,
                            ts::Vector,
                            g0::Vector,
                            n::Int,
                            l::Int,
                            lmax::Int=50)
    τ0 = ts[1]

    for (i,τ) in enumerate(ts)

        s = 0.0
        for k in 0:lmax
            term = K_nlk(n,l,k,τ,τ0) * g0[k + 1] # g0 already starts from l
            s += term
        end

        ρ[i] = s

    end
    
end






