function result=IntegerPartTest(phi_lower,phi_upper,scalar_term,flag)
% This function computes the intergration part of \sqrt(1-c^2/sin^2 \phi)
 
scalar_term=abs(scalar_term);
if flag==1 % sin term integration (c^2/sin^2 \phi)  
    lower=max(phi_lower,asin(scalar_term));
    upper=max(phi_upper,asin(scalar_term));
    if lower==upper
        result=0;
    else
        if scalar_term==0
            result= phi_upper-phi_lower;
        else
            res_part=@(phase_term) scalar_term*atan(cos(phase_term)/sqrt(sin(phase_term)^2/(scalar_term^2)-1)) ...
                     -1/scalar_term*asin(cos(phase_term)/sqrt(1-scalar_term^2));
            result=res_part(upper)-res_part(lower);
        end   
    end
else  % cos term integration (c^2/cos^2 \phi)    
    lower=min(phi_lower,acos(scalar_term));
    upper=min(phi_upper,acos(scalar_term));
    if lower==upper
        result=0;
    else
        if scalar_term==0
            result= phi_upper-phi_lower;
        else
            res_part=@(phase_term) -scalar_term*atan(sin(phase_term)/sqrt(cos(phase_term)^2/(scalar_term^2)-1)) ...
                     +asin(sin(phase_term)/sqrt(1-scalar_term^2));
            result=res_part(upper)-res_part(lower);
        end
    end
end

end
