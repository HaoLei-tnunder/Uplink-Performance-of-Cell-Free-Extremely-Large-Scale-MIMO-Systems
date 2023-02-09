clc
clear ;

nbrOfRealizations=400;
nbrOfSetups =2;
FRAME=6;
double=0; % the double scattering environment
step=10;

SNR=5;
var_noise=10^(-0.1*SNR);
% Pmax = 1/var_noise;
Pmax = 0.2 ;

lambda=0.03; % wavelength

APsNum=20;
M=APsNum;

Ns_X=6; Ns_Y=6;
TraNumNs=Ns_X*Ns_Y;
Nr_X=12; Nr_Y=12;
RecNumNr=Nr_X*Nr_Y;

%% Prepare

SE_MR_Level_2 = zeros(nbrOfSetups ,FRAME );
SE_MR_Level_1 = zeros(nbrOfSetups ,FRAME );

%% Go through all setups

for i = 1 : FRAME

    UsersNum=i+4;
    K=UsersNum;


    for n =  1 : nbrOfSetups

        if n == 1

            RecSpacing_num=3;
            TraSpacing_num=3;

        elseif n == 2

            RecSpacing_num=6;
            TraSpacing_num=6;

        end

        RecSpacing=lambda/RecSpacing_num;
        TraSpacing=lambda/TraSpacing_num;%Dleta

        fprintf('n = %u,i = %u\n ',n,i)

        [RecVarianceVec,TraVarianceVec,RecResVector,TraResVector,ns,nr] = generateSetup(M,K,Nr_X,Nr_Y,RecSpacing,Ns_X,Ns_Y,TraSpacing,lambda);

        [Channel] =  functionChannelGeneration(  RecVarianceVec,TraVarianceVec,RecResVector,TraResVector,M,K,ns,nr,RecNumNr,TraNumNs,nbrOfRealizations   );

        
        [ SE_MR_Level_2_item] = functionComputeMonteCarlo_SE_UL_Level_2(Channel,var_noise,M,K,TraNumNs,Pmax,nbrOfRealizations);

        [ SE_MR_Level_1_item] = functionComputeMonteCarlo_SE_UL_Level_1(Channel,var_noise,M,K,TraNumNs,Pmax,nbrOfRealizations);


        SE_MR_Level_2(n,i)  = sum(SE_MR_Level_2_item)/K;
        SE_MR_Level_1(n,i) =  sum(SE_MR_Level_1_item)/K;

    end
end


%% draw


SE_13=SE_MR_Level_1(1,:);
SE_16=SE_MR_Level_1(2,:);
SE_23=SE_MR_Level_2(1,:);
SE_26=SE_MR_Level_2(2,:);
  

figure;
hold on; box on;
plot(linspace(5,10,FRAME),(SE_23(:)),'d - r','LineWidth',2);
plot(linspace(5,10,FRAME),(SE_26(:)),'d -- r','LineWidth',2);
plot(linspace(5,10,FRAME),(SE_13(:)),'s - b','LineWidth',2);
plot(linspace(5,10,FRAME),(SE_16(:)),'s -- b','LineWidth',2);
legend('Distributed,$\Delta$=$\lambda/3$','Distributed,$\Delta$=$\lambda/6$','Small-cell,$\Delta$=$\lambda/3$','Small-cell,$\Delta$=$\lambda/6$',...
    'Interpreter','latex' )
xlabel('Number of UEs ','Interpreter','latex')
ylabel('Achievable average SE[bit/s/Hz]','Interpreter','latex')
xticks(5:1:10); 
grid on


