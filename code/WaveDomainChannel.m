function [RecVarianceVec,WD_RecX_vec,WD_RecY_vec,nr_act]=WaveDomainChannel(RecSamplePoint,WD_RecX_vec,WD_RecY_vec,lambda,RecLength)
% This function computes the variances of received vector
% Input: RecCor: the cell that includes the coordinates of sampled square;
%         SamplePoint: the logical matrix of sampled eclipse (0 or 1);
%         nr_act: the number of sampled eclipse points (may smaller than theo.)
% Output:

% Compute the variance of n_r variances
RecVarianceVec=arrayfun(  @(RecX_vec,RecY_vec)(IntVariance(RecX_vec,RecY_vec,lambda,RecLength))  ,  WD_RecX_vec ,   WD_RecY_vec,'UniformOutput',true);
%将离散后的求和点代入IntVariance计算方差
RecVarianceVec=abs(RecVarianceVec);   
RecZeroLoc=find(RecSamplePoint==0);
RecVarianceVec(RecZeroLoc)=[];
WD_RecX_vec(RecZeroLoc)=[];
WD_RecY_vec(RecZeroLoc)=[];
RecMiniLoc=find(RecVarianceVec<9e-4);
RecVarianceVec(RecMiniLoc)=[];
WD_RecX_vec(RecMiniLoc)=[];
WD_RecY_vec(RecMiniLoc)=[];
nr_act=length(RecVarianceVec);

% % Compute the variance of n_s variances
% TraVarianceVec=arrayfun(@(TraX_vec,TraY_vec)(IntVariance(TraX_vec,TraY_vec,lambda,TraLength)),TraX_vec,TraY_vec,'UniformOutput',true);
% TraVarianceVec=abs(TraVarianceVec);
% TraZeroLoc=find(TraVarianceVec==0);
% TraVarianceVec(TraZeroLoc)=[];
% TraX_vec(TraZeroLoc)=[];
% TraY_vec(TraZeroLoc)=[];
% TraMiniLoc=find(TraVarianceVec<9e-4);
% TraVarianceVec(TraMiniLoc)=[];
% TraX_vec(TraMiniLoc)=[];
% TraY_vec(TraMiniLoc)=[];
% ns_act=length(TraVarianceVec);



end









