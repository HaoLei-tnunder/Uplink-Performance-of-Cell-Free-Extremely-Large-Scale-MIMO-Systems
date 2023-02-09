function [SE_MR_th_Level_2  ,numerator, interfer,Gamma,Z  ] = functionComputeTheoretical_SE_UL_Level_2(R,RR,var_noise,M,K,RecNumNr,TraNumNs,p)


SE_MR_th_Level_2 = zeros(K,1);
numerator = zeros(TraNumNs,TraNumNs,K);
interfer = zeros(TraNumNs,TraNumNs,K);
noize=zeros(TraNumNs,TraNumNs,K);

Gamma=zeros(TraNumNs,TraNumNs,K);
Z = zeros(TraNumNs,TraNumNs,M,K);

% A = zeros(TraNumNs,TraNumNs,M,K,M,K);
B = zeros(TraNumNs,TraNumNs,M,K,M);
C = zeros(TraNumNs,TraNumNs,M,K,K);
D = zeros(TraNumNs,TraNumNs,M,K);
AB=zeros(TraNumNs,TraNumNs,K);


for m = 1 : M
    for k = 1 : K
        for i = 1 : TraNumNs
            for ii = 1 : TraNumNs

                Z(i,ii,m,k) = trace(R( (ii-1)*RecNumNr+1:ii*RecNumNr,(i-1)*RecNumNr+1:i*RecNumNr));

            end
        end

    end
end


for m = 1 : M
    for k = 1 : K
        for l = 1 : K
            for mm = 1 : M

                if (m~=mm) && (l~=k)


                elseif (m~=mm) && (l==k)

                    B(:,:,m,k,mm) = B(:,:,m,k,mm) +  Z(:,:,m,k) *  Z(:,:,mm,k) ;

                elseif (m==mm) && (l~=k)

                    for i = 1 : TraNumNs
                        for ii = 1 : TraNumNs
                                for x = 1 :TraNumNs

                                    C(i,ii,m,k,l) = C(i,ii,m,k,l) +  trace(  R( (x-1)*RecNumNr+1:x*RecNumNr,(x-1)*RecNumNr+1:x*RecNumNr)  ...
                                        * R( (ii-1)*RecNumNr+1:ii*RecNumNr,(i-1)*RecNumNr+1:i*RecNumNr)     );

                                end

                        end
                    end

                elseif (m==mm) && (l==k)

                    for i = 1 : TraNumNs
                        for x = 1:TraNumNs
                            for xx = 1 :TraNumNs
                                for nnn = 1:TraNumNs
                                    for nnnn = 1 : TraNumNs

                                        if any(nnn == nnnn)

                                            D(nnn,nnnn,m,k) = D(nnn,nnnn,m,k) + norm( RR( (x-1)*RecNumNr+1:x*RecNumNr ,(nnn-1)*RecNumNr+1:nnn*RecNumNr)*RR( (ii-1)*RecNumNr+1:ii*RecNumNr , (xx-1)*RecNumNr+1:xx*RecNumNr),'fro')^2 ...
                                                +     trace(  RR( (x-1)*RecNumNr+1:x*RecNumNr     ,       (nnn-1)*RecNumNr+1:nnn*RecNumNr)  ...
                                                *  RR( (i-1)*RecNumNr+1:i*RecNumNr       ,        (x-1)*RecNumNr+1:x*RecNumNr))  ...
                                                *     trace(  RR( (nnnn-1)*RecNumNr+1:nnnn*RecNumNr     ,(xx-1)*RecNumNr+1:xx*RecNumNr) ...
                                                *  RR( (xx-1)*RecNumNr+1:xx*RecNumNr      ,      (i-1)*RecNumNr+1:i*RecNumNr) )             ;
                                        else

                                            D(nnn,nnnn,m,k) = D(nnn,nnnn,m,k) +  trace(  RR( (x-1)*RecNumNr+1:x*RecNumNr  ,  (nnn-1)*RecNumNr+1:nnn*RecNumNr   )  ...
                                                *  RR( (i-1)*RecNumNr+1:i*RecNumNr       ,          (xx-1)*RecNumNr+1:xx*RecNumNr)  ...
                                                *  RR( (xx-1)*RecNumNr+1:xx*RecNumNr       ,           (i-1)*RecNumNr+1:i*RecNumNr) ...
                                                *  RR( (nnnn-1)*RecNumNr+1:nnnn*RecNumNr     ,           (x-1)*RecNumNr+1:x*RecNumNr) )...
                                                +     trace(  RR( (x-1)*RecNumNr+1:x*RecNumNr     ,       (nnn-1)*RecNumNr+1:nnn*RecNumNr)  ...
                                                *  RR( (i-1)*RecNumNr+1:i*RecNumNr       ,        (x-1)*RecNumNr+1:x*RecNumNr))  ...
                                                *     trace(  RR( (nnnn-1)*RecNumNr+1:nnnn*RecNumNr     ,(xx-1)*RecNumNr+1:xx*RecNumNr) ...
                                                *  RR( (xx-1)*RecNumNr+1:xx*RecNumNr      ,      (i-1)*RecNumNr+1:i*RecNumNr) )             ;

                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


for m = 1 : M
    for k = 1 : K
        for l = 1 : K
            for mm = 1 : M

                if (m~=mm) && (l~=k)

%                     AB(:,:,k) = AB(:,:,k) + A(:,:,m,k,mm,l);

                elseif (m~=mm) && (l==k)

                    AB(:,:,k) = AB(:,:,k) + B(:,:,m,k,mm) ;

                elseif (m==mm) && (l~=k)

                    AB(:,:,k) = AB(:,:,k) + C(:,:,m,k,l)  ;

                elseif (m==mm) && (l==k)

                    AB(:,:,k) = AB(:,:,k) + D(:,:,m,k) ;

                end

            end
        end
    end
end

for k = 1 : K

    Gamma(:,:,k) =  AB(:,:,k);

end

for m = 1 : M
    for k = 1 : K

        numerator(:,:,k) = numerator(:,:,k) + sqrt(p)*Z(:,:,m,k) ;

    end
end

for k = 1 : K

    noize(:,:,k) =noize(:,:,k) +  var_noise*Z(:,:,m,k) ;

end


for k = 1 : K

    interfer(:,:,k) =p* Gamma(:,:,k) - numerator(:,:,k)*numerator(:,:,k)' + noize(:,:,k);

end

for k = 1 : K

    SE_MR_th_Level_2(k) = real(   log2( det( eye(TraNumNs) + numerator(:,:,k)' * (pinv( interfer(:,:,k)) ) * numerator(:,:,k) )  )    );

end







end


