plot(τ,ρ[:,nₙ,1]./neq, 
        xaxis   =:log ,
        #yaxis   =:log ,
        xlabel  =   "τ",  
        ylabel  =   "T",
        dpi     =   300)
plot(τ,ρ[:,nₑ,1]./eeq, 
        xaxis   =:log ,
        #yaxis   =:log ,
        xlabel  =   "τ",  
        ylabel  =   "T",
        dpi     =   300)