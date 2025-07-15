using SpecialFunctions
using Integrals
using NLsolve
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

function rnl(u,v,p)
     
    n = p[1] 
    l = p[2]
    ζ = p[3]
    Λ = p[4] 
    ξ = p[5]
    ϑ = p[6]

    θ = (π/2 )*(1 + v )    # v ∈ [-1,1], θ \in [0,π]
    y = ζ + u             # u ∈ [0,∞], y \in [ζ,∞]

    return sin(θ) * ( cos(θ)^(2l) ) * ( y^(n-2l-2) ) * ( ( y^2 - ζ^2 )^(l+0.5) ) * exp( u - (sqrt( (y^2) - (y^2 -ζ^2) * (ξ^2) * (cos(θ)^2) + (1+ϑ)*ζ^2)/Λ) )
end

function Rnl(n,l,ζ=0.0001,Λ=1.0,ξ=0.0,ϑ=0.0)
    
    p = (n,l,ζ,Λ,ξ,ϑ)

    Int = 0.0
    for i in eachindex(uᵢ), j in eachindex(vᵢ)
        Int+= uwᵢ[i]*vwᵢ[j]*rnl(uᵢ[i],vᵢ[j],p)
    end
    return ((2l+1)*π*0.25 )*Int
end


function dGnl(n,l,ζ)
    if ζ == 0
        return 0
    else
        return   (n*Gnl(ζ,n,l) - Gnl(ζ,n+1,l))/ζ
    end
end

function Init_χ_Eq!(χ₀,nₐᵣ,L,m,T,α,Πₚ,πₚ)

    for (n,nv) in enumerate(nₐᵣ)    
        for l in 0:L-1  
            #println(nv," ",l)      
            χ₀[n,l+1] = 1#+ ρeq(T₀,μ₀,nv,0)*(l/(2*l + 1)^(2))
        end
    end
end

function Init_χ_AIso!(χ₀,nₐᵣ,L,p) # p = (ζ,α,Πₚ,πₚ)
    ζ   = p[0]
    α   = p[1]
    Πₚ  = p[2]
    πₚ  = p[3]
    Εₛ   = 1    # E/Eeq
    Pₛ   = 1    # P/Peq

    PL  = Pₛ  + Πₚ - πₚ
    PT  = Pₛ  + Πₚ - (πₚ/2)

    # Finding parameters for corresponding values of  Πₚ and πₚ

    


    for (n,nv) in enumerate(nₐᵣ)    
        for l in 0:L-1  
            #println(nv," ",l)      
            χ₀[n,l+1] = Rnl(n,l,ζ,Λ,ξ,ϑ)/Gnl(n,l,ζ)       # ρeq(T₀,μ₀,nv,0)*(l/(2*l + 1)^(2))
        end
    end
end




#---------------------------------------------------------------------
function trange(tspan,N,interpolation,scale= nothing)
    lin = range(0,stop=1,length=N)
    tₛ = tspan[1]
    tₑ = tspan[2]
    Δt = tₑ - tₛ 

    if interpolation == "exp" # Exponential spacing. Dense near first point.
        if scale == nothing
            scale = 1.0
        end
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
    G30 = Gnl(ζ,3.0,0)
    G40 = Gnl(ζ,4.0,0)
    G41 = Gnl(ζ,4.0,1.0)
    G50 = Gnl(ζ,5.0,0)
    
    Gd = G40^2 - G30*G50
    

    for (n,nᵥ) in enumerate(nₐᵣ)        
        #--------------------------------------------------------------------------
        for l = 1:L+1           
            #--------------------------------------------------------------------------

            # Here for lmax = L-1 we put the truncation condition.
            if l < L           
                h = (1/3)*( (( Gnl(ζ,nᵥ + 4.0,l-1)*G30 -  G40*Gnl(ζ,nᵥ + 3.0,l-1) )/(Gd))*(G41/Gnl(ζ,nᵥ + 3.0,l-1))*(χ[nₑ,2]* χ[n,l])  +  3*(nᵥ-2(l-1))*((2*(l-1) +1.0)/(2*(l-1) +3.0))*( Gnl(ζ,nᵥ + 3.0,l)/Gnl(ζ,nᵥ + 3.0,l-1) )* χ[n,l+1] )
                
                f = (2*(l-1))*χ[n,l] 
                r = ωᵣ*( χ[n,l] - 1 )

                dχ[n,l] = -( h/t)- ( f/t)  - ( r )  
            elseif l == L
                h = (1/3)*( (( Gnl(ζ,nᵥ + 4,l-1)*G30 -  G40*Gnl(ζ,nᵥ + 3,l-1) )/(Gd))*(G41/Gnl(ζ,nᵥ + 3,l-1))*(χ[nₑ,2]* χ[n,l])  +  3* (nᵥ-2(l-1)) *( (2*(l-1) +1)/(2*(l-1) +3) )*( Gnl(ζ,nᵥ + 3.0,l)/Gnl(ζ,nᵥ + 3.0,l-1) )* 1 )
                
                f = (2*(l-1))*χ[n,l] 
                r = ωᵣ*( χ[n,l] - 1 )

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

    if ζ₀ == 0
        ζ   = 0
    else
        ζ   = ζ₀/T         # m/T, (m/T₀)*(T₀/T) = 
    end
    

    w = uᵣ/u
    

    G30 = Gnl(ζ,3.0,0)
    G40 = Gnl(ζ,4.0,0)
    G41 = Gnl(ζ,4.0,1.0)
    G50 = Gnl(ζ,5.0,0)
    
    Gd = G40^2 - G30*G50
    #println("Gd   :",Gd )

    for (n,nᵥ) in enumerate(nₐᵣ)        
        #--------------------------------------------------------------------------
        for l = 1:L+1           
            #--------------------------------------------------------------------------

            # Here for lmax = L-1 we put the truncation condition.
            if l < L           
                h = (1/3)*( (( Gnl(ζ,nᵥ + 4.0,l-1)*G30 -  G40*Gnl(ζ,nᵥ + 3.0,l-1) )/(Gd))*(G41/Gnl(ζ,nᵥ + 3.0,l-1))*(χ[nₑ,2]* χ[n,l])  +  3*(nᵥ-2(l-1))*((2*(l-1) +1.0)/(2*(l-1) +3.0))*( Gnl(ζ,nᵥ + 3.0,l)/Gnl(ζ,nᵥ + 3.0,l-1) )* χ[n,l+1] )
                
                f = (2*(l-1))*χ[n,l] 
                r = ( χ[n,l] - 1 )

                dχ[n,l] = -( h +  f)*(w)  - ( r )  
            elseif l == L
                h = (1/3)*( (( Gnl(ζ,nᵥ + 4,l-1)*G30 -  G40*Gnl(ζ,nᵥ + 3,l-1) )/(Gd))*(G41/Gnl(ζ,nᵥ + 3,l-1))*(χ[nₑ,2]* χ[n,l])  +  3* (nᵥ-2(l-1)) *( (2*(l-1) +1)/(2*(l-1) +3) )*( Gnl(ζ,nᵥ + 3.0,l)/Gnl(ζ,nᵥ + 3.0,l-1) )* 1 )
                
                f = (2*(l-1))*χ[n,l] 
                r = ( χ[n,l] - 1 )

                dχ[n,l] = -( h +  f)*(w)  - ( r )   
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