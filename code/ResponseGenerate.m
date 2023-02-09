function [RecResVector,RecNumNr]=ResponseGenerate(spacing,lambda,RecLength,WD_RecX_vec,WD_RecY_vec,flag)

Nr_X=round(RecLength.L_x/spacing); % The elements in x axis
Nr_Y=round(RecLength.L_y/spacing); % The elements in y axis

RecPointTemp_X=(RecLength.distance(1):spacing:spacing*(Nr_X-1)+RecLength.distance(1));
RecPointTemp_Y=(spacing*(Nr_Y-1)+RecLength.distance(2):-spacing:RecLength.distance(2));
[RecPoint_X,RecPoint_Y]=meshgrid(RecPointTemp_X,RecPointTemp_Y); 

RecNumNr=Nr_X*Nr_Y;

RecPointVec_X=RecPoint_X(:);
RecPointVec_Y=RecPoint_Y(:);

 
gamma_z=2*pi/lambda*sqrt(1- (lambda/RecLength.L_x*WD_RecX_vec).^2 - ...
        (lambda/RecLength.L_y*WD_RecY_vec).^2);

switch(flag)
    case 'Transmit'
        temp_x=WD_RecX_vec*RecPointVec_X'*2*pi/RecLength.L_x;
        temp_y=WD_RecY_vec*RecPointVec_Y'*2*pi/RecLength.L_y;
        RecResVector=exp(-1i*(temp_x+temp_y+gamma_z*(RecLength.distance(3)*ones(1,RecNumNr)) ));
    case 'Receive'
        temp_x=RecPointVec_X*WD_RecX_vec'*2*pi/RecLength.L_x;
        temp_y=RecPointVec_Y*WD_RecY_vec'*2*pi/RecLength.L_y;
        RecResVector=exp(1i*(temp_x+temp_y+RecLength.distance(3)*ones(RecNumNr,1)*gamma_z' ));
end

end
