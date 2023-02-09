function [ SE_MR_Level_2] = functionComputeMonteCarlo_SE_UL_Level_2(H,var_noise,M,K,TraNumNs,p,nbrOfRealizations)

SE_MR_Level_2 = zeros(K,1);
numerator_E = zeros(TraNumNs,TraNumNs,K);
interfer_E = zeros(TraNumNs,TraNumNs,K);
noize_E=zeros(TraNumNs,TraNumNs,K);

Gamma_E=zeros(TraNumNs,TraNumNs,K);
Z_E = zeros(TraNumNs,TraNumNs,K);

for n = 1:nbrOfRealizations

    for m = 1 : M
        for k = 1 : K
            for l = 1 : K
                for mm = 1 : M

                    Gamma_E(:,:,k) = Gamma_E(:,:,k) + p *  H(:,:,m,k,n)' * H(:,:,m,l,n) * H(:,:,mm,l,n)' * H(:,:,mm,k,n)/nbrOfRealizations;

                end
            end

            Z_E(:,:,k) = Z_E(:,:,k) + H(:,:,m,k,n)'*H(:,:,m,k,n)/nbrOfRealizations;

        end
    end
    
end


    for k = 1 : K

        numerator_E(:,:,k) =  sqrt(p)* Z_E(:,:,k)  ;

        noize_E(:,:,k) =    var_noise * Z_E(:,:,k) ;

    end


for k = 1 : K

    interfer_E(:,:,k) =   Gamma_E(:,:,k)  - numerator_E(:,:,k)*numerator_E(:,:,k)'+ noize_E(:,:,k) ;

end

for k = 1 : K

   SE_MR_Level_2(k) = real( log2( det( eye(TraNumNs) + numerator_E(:,:,k)'* ( pinv  ( interfer_E(:,:,k)) ) * numerator_E(:,:,k) )  ) );
   
end


end


