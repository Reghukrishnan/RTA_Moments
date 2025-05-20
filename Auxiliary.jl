using SpecialFunctions
using Integrals
using NLsolve




#Integrant 
function g(y,p)
    ζ = p[1]
    n = p[2]
    l = p[3]
    return ( y^(n-2.0*l-2))*(y^2.0 - 1)^(l + 0.5)*exp(-ζ*y)
end

function Gnl(ζ,n,l)
    
    if ζ == 0
        return Γ(n)
    else
        domain  = (1.0,ζ*10 + 1000.0)
        
        prob    = IntegralProblem(g, domain,[ζ,n,l] )
        Gₙₗ      = solve(prob, QuadGKJL())
        return (ζ^(n+1))*Gₙₗ[1]
    end
end


function dGnl(ζ,n,l)
    if ζ == 0
        return 0
    else
        return   ((n+1)*Gnl(ζ,n,l) - Gnl(ζ,n+1,l))/ζ
    end
end

function Init_χ_Eq!(χ₀,nₐᵣ,L,m,T,α)
    for (n,nv) in enumerate(nₐᵣ)    
        for l in 0:L-1  
            #println(nv," ",l)      
            χ₀[n,l+1] = 1#+ ρeq(T₀,μ₀,nv,0)*(l/(2*l + 1)^(2))
        end
    end
end




#---------------------------------------------------------------------
function trange(tspan,N,interpolation,n = nothing)
    if interpolation == "exp"
        tₛ       = log(tspan[1])
        tₑ      = log(tspan[2])
        step    = (tₑ -  tₛ)/(N)
        
        return [exp(tₛ + (i)*step) for i in 0:N-1 ]

    elseif interpolation == "lin"
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
    elseif interpolation == "sig"
        #h       = (tspan[1] + tspan[2])/2
        a       = (8)
        A       = (tspan[1] + tspan[2])
        γ       = log((A - tspan[1])/tspan[1])/a        
        tₛ       = (0)
        tₑ      = (2)*a
        step    = (tₑ -  tₛ)/(N)
        return [A/((1) + exp(-γ*( (i)*step - a)) ) for i in 0:N-1 ]
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
  

    α = χ[nₙ,L+1]
    T = χ[nₑ,L+1]
    #println(χ[nₑ,1])
    
    #println(m, "\n\n\n\n\n")
    ζ = m/T
    
    
    
    ωᵣ = ωᵣ⁰*((T/T₀))
    G30 = Gnl(ζ,3.0,0)
    G40 = Gnl(ζ,4.0,0)
    G50 = Gnl(ζ,5.0,0)
    
    Gd = G40^2 - G30*G50
    

    for (n,nᵥ) in enumerate(nₐᵣ)        
        #--------------------------------------------------------------------------
        for l = 1:L+1               
            #--------------------------------------------------------------------------

            # Here for lmax = L-1 we put the truncation condition.
            if l < L           
                h = (1/3)*( (( Gnl(ζ,nᵥ + 4.0,l-1)*G30 -  G40*Gnl(ζ,nᵥ + 3.0,l-1) )/(Gd))*(G40/Gnl(ζ,nᵥ + 3.0,l-1))*(χ[nₑ,2]* χ[n,l])  +  3*(nᵥ-2(l-1))*((2*(l-1) +1.0)/(2*(l-1) +3.0))* χ[n,l+1] )
                
                f = (2*(l-1))*χ[n,l] 
                r = ωᵣ*( χ[n,l] - 1 )

                dχ[n,l] = -( h/t)- ( f/t)  - ( r )  
            elseif l == L
                h = (1/3)*( (( Gnl(ζ,nᵥ + 4,l-1)*G30 -  G40*Gnl(ζ,nᵥ + 3,l-1) )/(Gd))*(G40/Gnl(ζ,nᵥ + 3,l-1))*(χ[nₑ,2]* χ[n,l])  +  3*(nᵥ-2(l-1))*((2*(l-1) +1)/(2*(l-1) +3))* 1 )
                
                f = (2*(l-1))*χ[n,l] 
                r = ωᵣ*( χ[n,l] - 1 )

                dχ[n,l] = -( h/t)- ( f/t)  - ( r )  
            else
                dχ[n,L+1] = 0  
            end         
        end
                
    end
    dχ[nₑ,L+1] = ( G30*G40/Gd)*χ[nₑ,2]*χ[nₑ,L+1]/(3t)

    dχ[nₙ,L+1] = -( (G40^2)/Gd)*(χ[nₑ,2]/(3t))- (1/t)
    return dχ
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