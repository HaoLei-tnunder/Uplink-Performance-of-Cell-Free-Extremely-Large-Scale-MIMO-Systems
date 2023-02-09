function [RecLength,TraLength] = generateSetup_XLMIMO(M,K,Nr_X,Nr_Y,RecSpacing,Ns_X,Ns_Y,TraSpacing)

cellRange =  1000; 
APheigth = 12;  
UEheigth = 2; 

    for m = 1 : M
    
        RecLength(m).L_x=Nr_X*RecSpacing;
        RecLength(m).L_y=Nr_Y*RecSpacing;
        RecLength(m).distance=round([rand(1)*cellRange, APheigth ,rand(1)*cellRange]);
    
    end
    
    for k = 1 : K
    
        TraLength(k).L_x=Ns_X*TraSpacing;
        TraLength(k).L_y=Ns_Y*TraSpacing;
        TraLength(k).distance=round([rand(1)*cellRange,UEheigth,rand(1)*cellRange]);
    
    end


end





  








    
