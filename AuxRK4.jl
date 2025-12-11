
#---------------------------------------------------------------------
function trange2(tspan,N,interpolation,n = nothing)
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