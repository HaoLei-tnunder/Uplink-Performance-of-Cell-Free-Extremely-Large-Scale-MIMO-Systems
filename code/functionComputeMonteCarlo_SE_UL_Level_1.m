function [ SE_MR_Level_1   ] = functionComputeMonteCarlo_SE_UL_Level_1(H,var_noise,M,K,TraNumNs,p,nbrOfRealizations)


SE_MR = zeros(M,K);

for n = 1:nbrOfRealizations

    Gamma_E=zeros(TraNumNs,TraNumNs,M,K);
    D_E = zeros(TraNumNs,TraNumNs,M,K);

    for m = 1 : M
        for k = 1 : K
            for l = 1 : K

                if l ~= k

                    Gamma_E(:,:,m,k) = Gamma_E(:,:,m,k) +p * H(:,:,m,k,n)' * H(:,:,m,l,n) * H(:,:,m,l,n)'* H(:,:,m,k,n) ;

                end
            end

            Gamma_E(:,:,m,k) = Gamma_E(:,:,m,k) + var_noise *H(:,:,m,k,n)' *H(:,:,m,k,n);
        end

    end

    for m = 1 : M
        for k = 1 : K

            D_E(:,:,m,k) = sqrt(p)* H(:,:,m,k,n)'*H(:,:,m,k,n);

        end
    end


    for m = 1 : M
        for k = 1 : K

            SE_MR(m,k) = SE_MR(m,k) + real(log2( det( eye(TraNumNs) + D_E(:,:,m,k)'*( pinv  ( Gamma_E(:,:,m,k)) ) * D_E(:,:,m,k) )))/nbrOfRealizations ;

        end
    end

end

SE_MR_Level_1 = max(   SE_MR,[],1  );

end


