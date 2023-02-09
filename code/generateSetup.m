function [RecVarianceVec,TraVarianceVec,RecResVector,TraResVector,ns,nr] = generateSetup(M,K,Nr_X,Nr_Y,RecSpacing,Ns_X,Ns_Y,TraSpacing,lambda)

%     double = 0;

RecNumNr=Nr_X*Nr_Y;
TraNumNs=Ns_X*Ns_Y;

[RecLength,TraLength] = generateSetup_XLMIMO(M,K,Nr_X,Nr_Y,RecSpacing,Ns_X,Ns_Y,TraSpacing);

WD_RecX_vec=[];
WD_RecY_vec=[];
WD_TraX_vec=[];
WD_TraY_vec=[];
RecVarianceVec=[];
TraVarianceVec=[];
TraResVector=[];
RecResVector=[];

% Obtain the sampled eclipse and their coordinates given a surface
for m = 1 : M

    [RecSamplePoint,WD_RecX_vec_item,WD_RecY_vec_item]=WD_eclipse_sample(lambda,RecLength(m));

    [RecVarianceVec_item,WD_RecX_vec_item,WD_RecY_vec_item,nr_act]=WaveDomainChannel(RecSamplePoint,WD_RecX_vec_item,WD_RecY_vec_item,lambda,RecLength(m));

    %     ChannelSigma_item=RecVarianceVec_item*TraVarianceVec_item';

    RecVarianceVec=[RecVarianceVec,RecVarianceVec_item];
    WD_RecX_vec=[WD_RecX_vec,WD_RecX_vec_item];
    WD_RecY_vec=[WD_RecY_vec,WD_RecY_vec_item];

end


for k = 1 : K

    [TraSamplePoint,WD_TraX_vec_item,WD_TraY_vec_item]=WD_eclipse_sample(lambda,TraLength(k));

    [TraVarianceVec_item,WD_TraX_vec_item,WD_TraY_vec_item,ns_act]=WaveDomainChannel(TraSamplePoint,WD_TraX_vec_item,WD_TraY_vec_item,lambda,TraLength(k));

    % average
    %         if double~=1
    %             tt1=mean(TraVarianceVec_item);
    %             TraVarianceVec_item=tt1*ones(ns_act,1);
    %         end

    TraVarianceVec=[TraVarianceVec,TraVarianceVec_item];
    WD_TraX_vec=[WD_TraX_vec,WD_TraX_vec_item];
    WD_TraY_vec=[WD_TraY_vec,WD_TraY_vec_item];

end

TraResVector=zeros(ns_act,TraNumNs,K);
RecResVector=zeros(RecNumNr,nr_act,M);


for m = 1 : M
    [RecResVector(:,:,m)]=ResponseGenerate(RecSpacing,lambda,RecLength(m),WD_RecX_vec(:,m),WD_RecY_vec(:,m),'Receive');
end

for k = 1 : K
    %响应矢量
    [TraResVector(:,:,k)]=ResponseGenerate(TraSpacing,lambda,TraLength(k),WD_TraX_vec(:,k),WD_TraY_vec(:,k),'Transmit'); %响应矢量

end

ns = ns_act;
nr = nr_act;



end


