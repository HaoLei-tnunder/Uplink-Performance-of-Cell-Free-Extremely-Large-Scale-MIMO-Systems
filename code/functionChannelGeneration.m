function [Channel,bate] = functionChannelGeneration(  RecVarianceVec,TraVarianceVec,RecResVector,TraResVector,M,K,ns,nr,RecNumNr,TraNumNs,nbrOfRealizations , betaVal  )

% R = zeros(RecNumNr*TraNumNs,RecNumNr*TraNumNs);
% RR = zeros(RecNumNr*TraNumNs,RecNumNr*TraNumNs);


Channel=zeros(RecNumNr,TraNumNs,M,K,nbrOfRealizations);
bate = zeros(M,K);


for n = 1 : nbrOfRealizations

    ChannelRandom=zeros(nr,ns,M,K);
    ChannelHa=zeros(nr,ns,M,K);

    for m = 1 : M
        for k = 1 : K

            ChannelRandom(:,:,m,k)=sqrt(1/2)*(randn(nr,ns)+1i*randn(nr,ns));       %w
            ChannelHa(:,:,m,k)= diag(RecVarianceVec(:,m))*ChannelRandom(:,:,m,k)*diag(TraVarianceVec(:,k));
%             ChannelHa(:,:,m,k)=sqrt(betaVal(m,k))*diag(RecVarianceVec(:,m))*ChannelRandom(:,:,m,k)*diag(TraVarianceVec(:,k));            
            Channel(:,:,m,k,n) =RecResVector(:,:,m)* ChannelHa(:,:,m,k)* TraResVector(:,:,k);

%             H = Channel(:,:,m,k,n);
%             bate(m,k) = 1/(TraNumNs*RecNumNr)*trace( H(:)*H(:)'  );
            
        end
    end

end







% for m = 1 : M
%     for k = 1 : K
% m=1;k=1;
%         R_lambda_r=diag(reshape(RecVarianceVec(:,m).*RecVarianceVec(:,m),[],1));
%         R_Cor_r=RecResVector(:,:,m) *R_lambda_r*RecResVector(:,:,m)';
%         R_lambda_s=diag(reshape(TraVarianceVec(:,k).*TraVarianceVec(:,k),[],1));
%         R_Cor_s=TraResVector(:,:,k)'*R_lambda_s*TraResVector(:,:,k);
%         R =kron(R_Cor_s,R_Cor_r)      ;%R_H
%         RR = (R )^(1/2) ;
% 
%         RRR = kron(R_lambda_s,R_lambda_r)      ;%R_Ha
%         RRRR = (R)^(1/2)     ;        


%     end
% end


end


