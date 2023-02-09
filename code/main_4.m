clc
clear ;

nbrOfRealizations=400;
nbrOfSetups = 2 ;

SNR=5;
var_noise=10^(-0.1*SNR);
Pmax = 0.2;

lambda=0.03; % wavelength

UsersNum=8;
K=UsersNum;
UserSpacing=100*lambda;

Nr_X=12; Nr_Y=12;
Ns_X=12; Ns_Y=12;

RecNumNr=Nr_X*Nr_Y;
TraNumNs=Ns_X*Ns_Y;

%% Prepare

SE_MR_Level_2 = zeros(UsersNum,nbrOfSetups );
SE_MR_Level_1 = zeros(UsersNum,nbrOfSetups );

%% Go through all setups

for n = 1 %  1 : nbrOfSetups
for APS = 20 %  15:1:20
    M=APS;

    if n == 1

        TraSpacing_num=3;
        RecSpacing_num=6;

    elseif n == 2

        TraSpacing_num=3;
        RecSpacing_num=9;

    elseif n == 3

        TraSpacing_num=3;
        RecSpacing_num=9;

    elseif n == 4

        TraSpacing_num=6;
        RecSpacing_num=3;

    elseif n == 5

        TraSpacing_num=6;
        RecSpacing_num=6;   

    elseif n == 6

        TraSpacing_num=6;
        RecSpacing_num=9;

    elseif n == 7

        TraSpacing_num=9;
        RecSpacing_num=3;

    elseif n == 8

        TraSpacing_num=9;
        RecSpacing_num=6;

    elseif n == 9

        TraSpacing_num=9;
        RecSpacing_num=9;

    end

    fprintf(' M = %u, n = %u\n',M,n)

    RecSpacing=lambda/RecSpacing_num;
    TraSpacing=lambda/TraSpacing_num;%Dleta

    [RecVarianceVec,TraVarianceVec,RecResVector,TraResVector,ns,nr] = generateSetup(M,K,Nr_X,Nr_Y,RecSpacing,Ns_X,Ns_Y,TraSpacing,lambda);

    [Channel] =  functionChannelGeneration(  RecVarianceVec,TraVarianceVec,RecResVector,TraResVector,M,K,ns,nr,RecNumNr,TraNumNs,nbrOfRealizations   );
        
    [ SE_MR_Level_2(:,n)] = functionComputeMonteCarlo_SE_UL_Level_2(Channel,var_noise,M,K,TraNumNs,Pmax,nbrOfRealizations);

%     [ SE_MR_Level_1(:,n,M-14)] = functionComputeMonteCarlo_SE_UL_Level_1(Channel,var_noise,M,K,TraNumNs,Pmax,nbrOfRealizations);


end
end


%% draw


SE_Level_2=sum(SE_MR_th_Level_2);
SE_33=reshape(SE_Level_2(1,1,:),[8,1]);
SE_36=reshape(SE_Level_2(1,2,:),[8,1]);
SE_39=reshape(SE_Level_2(1,3,:),[8,1]);
SE_63=reshape(SE_Level_2(1,4,:),[8,1]);
SE_66=reshape(SE_Level_2(1,5,:),[8,1]);
SE_69=reshape(SE_Level_2(1,6,:),[8,1]);
SE_93=reshape(SE_Level_2(1,7,:),[8,1]);
SE_96=reshape(SE_Level_2(1,8,:),[8,1]);
SE_99=reshape(SE_Level_2(1,9,:),[8,1]);

figure;
hold on; box on;
plot(linspace(3,10,8),(SE_33),'d b - ','LineWidth',2);
plot(linspace(3,10,8),(SE_36),'d b --','LineWidth',2);
% plot(linspace(3,10,8),(SE_39),'s b -.','LineWidth',2);
plot(linspace(3,10,8),(SE_63),'^ r -','LineWidth',2);
plot(linspace(3,10,8),(SE_66),'^ r -.','LineWidth',2);
plot(linspace(3,10,8),(SE_69),'^ r --','LineWidth',2);
% plot(linspace(3,10,8),(SE_93),'x k -','LineWidth',2);
plot(linspace(3,10,8),(SE_96),'s k -','LineWidth',2);
plot(linspace(3,10,8),(SE_99),'s k --','LineWidth',2);
% legend('33','36','39','63','66','69','93','96','99','Interpreter','latex' )
legend('$\Delta_s$=$\lambda/3$,$\Delta_r$=$\lambda/3$','$\Delta_s$=$\lambda/3$,$\Delta_r$=$\lambda/6$',...
    '$\Delta_s$=$\lambda/6$,$\Delta_r$=$\lambda/3$','$\Delta_s$=$\lambda/6$,$\Delta_r$=$\lambda/3$','$\Delta_s$=$\lambda/6$,$\Delta_r$=$\lambda/9$',...
    '$\Delta_s$=$\lambda/9$,$\Delta_r$=$\lambda/6$','$\Delta_s$=$\lambda/9$,$\Delta_r$=$\lambda/9$','Interpreter','latex' )
xlabel('Number of BSs ','Interpreter','latex')
ylabel('Achievable sum SE[bit/s/Hz]','Interpreter','latex')
xticks(3:1:8); 
grid on



