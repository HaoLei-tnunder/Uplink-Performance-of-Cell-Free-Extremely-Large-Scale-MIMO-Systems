function ScalarVar=IntVariance(RecX,RecY,lambda,RecLength)
% This function computes variance in four orthants
SampCord=[RecX,RecY];
a=lambda*SampCord(1)/RecLength.L_x;
b=lambda*(SampCord(1)+1)/RecLength.L_x;
c=lambda*SampCord(2)/RecLength.L_y;
d=lambda*(SampCord(2)+1)/RecLength.L_y;

OrthFlagX=SampCord(1)<0; % the negative axis is 1
OrthFlagY=SampCord(2)<0;
ComXY=abs(SampCord(1))>abs(SampCord(2)); % ell_x>ell_y is 1
orthant=strcat(string(OrthFlagX*1),string(OrthFlagY*1));
point=string(ComXY*1);

switch(orthant)
    case '00' % the first orthant
        phi_1=atan(c/b); phi_2=min(atan(c/a),atan(d/b));
        phi_3=max(atan(c/a),atan(d/b));phi_4=atan(d/a);
        switch(point) 
            case '1' % ell_x>ell_y is 1
                % flag==1,sin; flag==0,cos
                temp1_11=IntegerPart(phi_1,phi_2,c,1);
                temp1_22=IntegerPart(phi_1,phi_2,b,0);
                temp1=temp1_11-temp1_22;
                temp2_11=IntegerPart(phi_2,phi_3,a,0);
                temp2_22=IntegerPart(phi_2,phi_3,b,0);
                temp2=temp2_11-temp2_22;
                temp3_11=IntegerPart(phi_3,phi_4,a,0);
                temp3_22=IntegerPart(phi_3,phi_4,d,1);
                temp3=temp3_11-temp3_22;
                ScalarVar=temp1+temp2+temp3;
            case '0' % ell_x<ell_y is 0
                temp1=IntegerPart(phi_1,phi_2,c,1) ...
                      -IntegerPart(phi_1,phi_2,b,0);
                temp2=IntegerPart(phi_2,phi_3,c,1) ...
                      -IntegerPart(phi_2,phi_3,d,1);
                temp3=IntegerPart(phi_3,phi_4,a,0) ...
                      -IntegerPart(phi_3,phi_4,d,1);
                ScalarVar=temp1+temp2+temp3;
        end
        
        case '10' % the second orthant
        phi_1=pi-atan(d/abs(b)); phi_2=pi-max(atan(c/abs(b)),atan(d/abs(a)));
        phi_3=pi-min(atan(c/abs(b)),atan(d/abs(a)));phi_4=pi-atan(c/abs(a));
        switch(point) 
            case '1' % ell_x>ell_y is 1
                % flag==1,sin; flag==0,cos
                temp1=IntegerPart(pi-phi_1,pi-phi_2,b,0) ...
                    -IntegerPart(pi-phi_1,pi-phi_2,d,1);
                temp2=IntegerPart(pi-phi_2,pi-phi_3,b,0) ...
                    -IntegerPart(pi-phi_2,pi-phi_3,a,0);
                temp3=IntegerPart(pi-phi_3,pi-phi_4,c,1) ...
                      -IntegerPart(pi-phi_3,pi-phi_4,a,0);
                ScalarVar=temp1+temp2+temp3;
            case '0' % ell_x<=ell_y is 0
                temp1=IntegerPart(pi-phi_1,pi-phi_2,b,0) ...
                    -IntegerPart(pi-phi_1,pi-phi_2,d,1);
                temp2=IntegerPart(pi-phi_2,pi-phi_3,c,1) ...
                    -IntegerPart(pi-phi_2,pi-phi_3,d,1);
                temp3=IntegerPart(pi-phi_3,pi-phi_4,c,1) ...
                      -IntegerPart(pi-phi_3,pi-phi_4,a,0);
                ScalarVar=temp1+temp2+temp3;
        end
        
        case '11' % the third orthant
        phi_1=pi+atan(abs(d)/abs(a)); phi_2=pi+min(atan(abs(d)/abs(b)),atan(abs(c)/abs(a)));
        phi_3=pi+max(atan(abs(d)/abs(b)),atan(abs(c)/abs(a)));phi_4=pi+atan(abs(c)/abs(b));
        switch(point) 
            case '1' % ell_x>ell_y is 1
                % flag==1,sin; flag==0,cos
                temp1=IntegerPart(phi_1-pi,phi_2-pi,d,1) ...
                    -IntegerPart(phi_1-pi,phi_2-pi,a,0);
                temp2=IntegerPart(phi_2-pi,phi_3-pi,b,0) ...
                    -IntegerPart(phi_2-pi,phi_3-pi,a,0);
                temp3=IntegerPart(phi_3-pi,phi_4-pi,b,0) ...
                    -IntegerPart(phi_3-pi,phi_4-pi,c,1);
                ScalarVar=temp1+temp2+temp3;
            case '0' % ell_x<=ell_y is 0
                temp1=IntegerPart(phi_1-pi,phi_2-pi,d,1) ...
                    -IntegerPart(phi_1-pi,phi_2-pi,a,0);
                temp2=IntegerPart(phi_2-pi,phi_3-pi,d,1) ...
                    -IntegerPart(phi_2-pi,phi_3-pi,c,1);
                temp3=IntegerPart(phi_3-pi,phi_4-pi,b,0) ...
                    -IntegerPart(phi_3-pi,phi_4-pi,c,1);
                ScalarVar=temp1+temp2+temp3;
        end
        
        
        case '01' % the fourth orthant
        phi_1=2*pi-atan(abs(c)/abs(a)); phi_2=2*pi-max(atan(abs(d)/abs(a)),atan(abs(c)/abs(b)));
        phi_3=2*pi-min(atan(abs(d)/abs(a)),atan(abs(c)/abs(b)));phi_4=2*pi-atan(abs(d)/abs(b));
        switch(point) 
            case '1' % ell_x>ell_y is 1
                % flag==1,sin; flag==0,cos
                temp1=IntegerPart(2*pi-phi_1,2*pi-phi_2,a,0) ...
                    -IntegerPart(2*pi-phi_1,2*pi-phi_2,c,1);
                temp2=IntegerPart(2*pi-phi_2,2*pi-phi_3,a,0) ...
                    -IntegerPart(2*pi-phi_2,2*pi-phi_3,b,0);
                temp3=IntegerPart(2*pi-phi_3,2*pi-phi_4,d,1) ...
                    -IntegerPart(2*pi-phi_3,2*pi-phi_4,b,0);
                ScalarVar=temp1+temp2+temp3;
            case '0' % ell_x<ell_y is 0
                temp1=IntegerPart(2*pi-phi_1,2*pi-phi_2,a,0) ...
                    -IntegerPart(2*pi-phi_1,2*pi-phi_2,c,1);
                temp2=IntegerPart(2*pi-phi_2,2*pi-phi_3,d,1) ...
                    -IntegerPart(2*pi-phi_2,2*pi-phi_3,c,1);
                temp3=IntegerPart(2*pi-phi_3,2*pi-phi_4,d,1) ...
                    -IntegerPart(2*pi-phi_3,2*pi-phi_4,b,0);
                ScalarVar=temp1+temp2+temp3;
        end
 

end