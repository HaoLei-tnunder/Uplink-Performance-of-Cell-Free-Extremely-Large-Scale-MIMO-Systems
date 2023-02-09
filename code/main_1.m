clc
clear ;

nbrOfRealizations=500;
nbrOfSetups =1;
setups = 1;

SNR=5;
var_noise=10^(-0.1*SNR);

Pmax = 0.2;
% Pmax = 1/var_noise ;

lambda=0.03; % wavelength

APsNum=1;
M=APsNum;
UsersNum=1;
K=UsersNum;

RecSpacing_num=3;
TraSpacing_num=3;
RecSpacing=lambda/RecSpacing_num;
TraSpacing=lambda/TraSpacing_num;%Dleta

Ns_X=6; Ns_Y=6;
TraNumNs=Ns_X*Ns_Y;

%% Prepare

SE_MR_Level_2 = zeros(UsersNum,nbrOfSetups,setups);
SE_MR_Level_1 = zeros(UsersNum,nbrOfSetups,setups );

%% Go through all setups
for jj = 1 : setups
for n =  1 : nbrOfSetups

    Nr_X=n+6; Nr_Y=n+6;
    RecNumNr=Nr_X*Nr_Y;

    fprintf('jj=%u,n=%u\n',jj,n)

    [RecVarianceVec,TraVarianceVec,RecResVector,TraResVector,ns,nr] = generateSetup(M,K,Nr_X,Nr_Y,RecSpacing,Ns_X,Ns_Y,TraSpacing,lambda);


%     x = length(RecVarianceVec(:,1));
%     for xx = 1 :  x
%         if RecVarianceVec(xx,1) > 100
%             RecVarianceVec(xx,1) =  RecVarianceVec(xx-3,1)/2+ RecVarianceVec(xx+3,1)/2;
%             for iii = 2 : M
%                 RecVarianceVec(xx,iii) =RecVarianceVec(xx,1);
%             end
%         end
%     end

    [Channel] =  functionChannelGeneration(  RecVarianceVec,TraVarianceVec,RecResVector,TraResVector,M,K,ns,nr,RecNumNr,TraNumNs,nbrOfRealizations  );

    [ SE_MR_Level_2(:,n,jj) ] = functionComputeMonteCarlo_SE_UL_Level_2(Channel,var_noise,M,K,TraNumNs,Pmax,nbrOfRealizations);

    [ SE_MR_Level_1(:,n,jj) ] = functionComputeMonteCarlo_SE_UL_Level_1(Channel,var_noise,M,K,TraNumNs,Pmax,nbrOfRealizations);



end
end

%% draw

% SE_1=sum(SE_MR_Level_1);
% SE_2=sum(SE_MR_Level_2);
% 
% 
% figure;
% hold on; box on;
% plot(linspace(4,12,9),(SE_2(1:9)),'d - r','LineWidth',2);
% plot(linspace(4,12,9),( SE_1(1:9)  ),'S - b','LineWidth',2);
% plot(linspace(4,12,9),( SE_2(1:9)  ),'o  k','LineWidth',2);
% legend('Cell-Free','Small-cell','Analytical' ,'Interpreter','latex' )
% xlabel('$N_{Hr}$=$N_{Vr}$ ','Interpreter','latex')
% ylabel('Achievable sum SE [bit/s/Hz]','Interpreter','latex')
% xticks(4:1:12);
% ylim([0,160])
% grid on




