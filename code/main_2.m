clc
clear ;

nbrOfRealizations=400;
nbrOfSetups =2;

SNR=5;
var_noise=10^(-0.1*SNR);
Pmax = 0.2;

lambda=0.03; % wavelength

UsersNum=8;
K=UsersNum;
UserSpacing=100*lambda;

Nr_X=12; Nr_Y=12;
Ns_X=6; Ns_Y=6;

RecNumNr=Nr_X*Nr_Y;
TraNumNs=Ns_X*Ns_Y;

%% Prepare

SE_MR_Level_2 = zeros(UsersNum,nbrOfSetups,6);
SE_MR_Level_1 = zeros(UsersNum,nbrOfSetups,6);


%% Go through all setups

for n = 1 : nbrOfSetups
    for APS = 15:1:20
        M=APS;

        if n == 1

            RecSpacing_num=3;
            TraSpacing_num=3;

        elseif n == 2

            RecSpacing_num=6;
            TraSpacing_num=6;

        end

        fprintf(' M = %u, n = %u\n',M,n)

        RecSpacing=lambda/RecSpacing_num;
        TraSpacing=lambda/TraSpacing_num;%Dleta

        [RecVarianceVec,TraVarianceVec,RecResVector,TraResVector,ns,nr] = generateSetup(M,K,Nr_X,Nr_Y,RecSpacing,Ns_X,Ns_Y,TraSpacing,lambda);

        [Channel] =  functionChannelGeneration(  RecVarianceVec,TraVarianceVec,RecResVector,TraResVector,M,K,ns,nr,RecNumNr,TraNumNs,nbrOfRealizations   );

        [ SE_MR_Level_2(:,n,M-14)] = functionComputeMonteCarlo_SE_UL_Level_2(Channel,var_noise,M,K,TraNumNs,Pmax,nbrOfRealizations);

        [ SE_MR_Level_1(:,n,M-14)] = functionComputeMonteCarlo_SE_UL_Level_1(Channel,var_noise,M,K,TraNumNs,Pmax,nbrOfRealizations);



    end
end


%% draw

SE_Level_1=sum(SE_MR_Level_1);
SE_Level_2=sum(SE_MR_Level_2);
SE_13=reshape(SE_Level_1(1,1,:),[6,1]);
SE_16=reshape(SE_Level_1(1,2,:),[6,1]);
SE_23=reshape(SE_Level_2(1,1,:),[6,1]);
SE_26=reshape(SE_Level_2(1,2,:),[6,1]);


figure;
hold on; box on;
plot(linspace(15,20,6),(SE_23),'d r -','LineWidth',2);
plot(linspace(15,20,6),(SE_26),'d r -- ','LineWidth',2);
plot(linspace(15,20,6),(SE_13),'s b -','LineWidth',2);
plot(linspace(15,20,6),(SE_16),'s b --','LineWidth',2);
legend('Distributed,$\Delta$=$\lambda/3$','Distributed,$\Delta$=$\lambda/6$','Small-cell,$\Delta$=$\lambda/3$','Small-cell,$\Delta$=$\lambda/6$',...
    'Interpreter','latex' )
xlabel('Number of BSs ','Interpreter','latex')
ylabel('Achievable sum SE[bit/s/Hz]','Interpreter','latex')
xticks(15:1:20);
ylim([0,180])
grid on



