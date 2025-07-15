plot(ξ ,T)
plot(τₛ,ξ ,
        #xaxis = :log,
        xticks  =   [0,1,2,3,4,5,6,7,8,9,10], 
        yticks  =   [0,1,2,3,4,5,6,7,8,9,10],
        xlabel  = "τ",
        ylabel  = "ξ")
plot(ξ[begin+1:1000] ,T[begin+1:1000],
                xaxis   =:log)
plot(ξ[begin+1:end] ,T[begin+1:end],
        xaxis   =:log)

plot(w[begin:1000],ξ[begin:1000])
plot(w,uᵣ, 
        xaxis   =:log,
        yaxis   =:log)

        plot!( ξ ,w,
        xticks  =   [0,1,2,3,4,5,6,7,8,9,10], 
        yticks  =   [0,1,2,3,4,5,6,7,8,9,10])


        plot!(ξ ,T,
        xticks  =   [0,1,2,3,4,5,6,7,8,9,10])

        plot!( ξ ,w.-w₀,
        xticks  =   [0,1,2,3,4,5,6,7,8,9,10], 
        yticks  =   [0,1,2,3,4,5,6,7,8,9,10],
        xlabel  = "ξ",
        ylabel  = "w")


        plot!( ξ[begin+1:end],r[begin+1:end],
        xticks  =   [0,1,2,3,4,5,6,7,8,9,10], 
        #yticks  =   [0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0],
        label   = "$(ζ₀),$(uᵣ⁰) ",
        xlabel  = "ξ",
        ylabel  = "(w.-w₀)/ξ")




function DSRTA(dχ::Matrix,χ,t::Float64,p)

    # p = (N,L,η₀,κ₀,nₐᵣ,nₙ,nₑ)
    N   = p[1]
    L   = p[2]
    η₀   = p[3]
    κ₀  = p[4]
    nₐᵣ = p[5]
    nₙ  = p[6]
    nₑ  = p[7]
  

    α = χ[nₙ,L+1]       # μ/T
    h = χ[nₑ,L+1]       # ln(τ₀T)
    u = χ[nₑ+1,L+1]     # τ/τ₀
    #println(χ[nₑ,1])
    
    #println(m, "\n\n\n\n\n")
    uᵣ  = 5*η₀*exp(-h)  # τᵣ/τ_₀ = η₀/(τ_₀T)

    if κ₀ == 0
        ζ = 0
    else
        ζ   = κ₀*uᵣ         # m/T, κ₀ = mτ₀/5η₀, uᵣ = 5*η₀/τ₀T
    end
    

    #println("h   :",h)
    #println("uᵣ   :",uᵣ)
    #println("ζ   :",ζ)

    w = uᵣ/u
    #println("w   :",w)

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
