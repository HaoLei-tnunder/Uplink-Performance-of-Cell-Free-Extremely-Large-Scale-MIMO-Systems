function [RecSamplePoint,WD_x_temp_vec,WD_y_temp_vec]=WD_eclipse_sample(lambda,RecLength)
% eclipse_sample: Create a sampled eclipse in wavenumber domain
% Input: lambda is wavelength; RecLength is the structure that includes the
%        length in x and length in y axis;
% Output: RecCor: the cell that includes the coordinates of sampled square;
%         SamplePoint: the logical matrix of sampled eclipse (0 or 1);
%         nr_act: the number of sampled eclipse points (may smaller than theo.)

x_range=(-RecLength.L_x/lambda:1:RecLength.L_x/lambda-1); % points in kx axis
y_range=(RecLength.L_y/lambda-1:-1:-RecLength.L_y/lambda);% points in ky axis
x_range=round(x_range);  
y_range=round(y_range);
% 四舍五入

[x_temp,y_temp]=meshgrid(x_range,y_range); 

WD_x_temp_vec=x_temp(:);
WD_y_temp_vec=y_temp(:);
RecCordNormVec=(WD_x_temp_vec*lambda/RecLength.L_x).^2+(WD_y_temp_vec*lambda/RecLength.L_y).^2;
%离散后的求和取值范围
RecSamplePoint=RecCordNormVec<=1;

% nr_act=sum(sum(RecSamplePoint));
% nr_theo=pi*RecLength.L_x*RecLength.L_y/(lambda^2);

end