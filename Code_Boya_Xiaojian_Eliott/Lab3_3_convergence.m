clear all; close all; clc;
addpath(genpath('libs'));
addpath(genpath('myFunctions'));
%%
clear all; close all;
addpath(genpath('myFunctions'));
%% Tx signal parameters

bit_number = 2*10000;              % length of stream of bits
bps = 2e6;
Nbps = 2;                                % #Bits per symbol
FC = 2e9; 
T = 1/(bps/Nbps);                         % symbol interval 1e-6;         
fs = 1/T;
df = 30.*(T*FC); 
phi = 3.14;

%% downsampling / upsampling factors
M = 100;  % Factor
t0 = 10; % time offset sequence
    
%% Root raised cosine Filter parameters
beta = 0.3;     % roll-off factor
RRCtaps=8*M+1;     % taps of rrc filter (root raised cosine)

%% Initinalization for plot

    if Nbps==1                  
        modulation='pam'; %BPSK
    else
        modulation='qam'; 
    end
N = bit_number/Nbps;        %length of Rx symbol
Realiations = 50;
error = zeros(N,3,Realiations); %for mean the error
symb_rx = zeros(size(error));   
for R=1:Realiations    
%% Transmiter Side
        bitStreamTx = fix(rand(bit_number,1)*2);% bitStream
        symb_tx = mapping(bitStreamTx,Nbps,modulation);
        
%         symb_R = real(symb_tx);
%         symb_I = imag(symb_tx);
%         figure(1)
%         scatter(symb_R,symb_I);
%         title('Constellation of modulated signal');  
        
        symb_tx_upsampled = upsample(symb_tx,M);
        [H_RC,freqGrid,h_rrc,h_rc,timeGrid] = HalfNyquistFilter(beta,T,M,RRCtaps);
        Tx_signal = conv(h_rrc,symb_tx_upsampled);


%% Reciever Side
        %BaseBand equivalent of an ideal channel
        EbN0_db = 400;
        noised_Rx_Signal = addAWGN(Tx_signal,EbN0_db,bit_number,M/T);
        Rx_signal = conv(h_rrc, Tx_signal);
        Rx_signal = Rx_signal(RRCtaps:end-(RRCtaps-1)); % usefull part from t=0
        %Rx_Signal = Tx_signal(RRCtaps/2-0.5:end-(RRCtaps/2-0.5)-1);
        
        %time shift
        %Rx_signal_shift =  [zeros(t0,1); Rx_signal];
        %Rx_signal_new = Rx_signal_new(1:length(Rx_signal)); %cutoff or use circshift(s,[]);
        Rx_signal_shift = circshift(Rx_signal, t0);
        
        %gardener
        K = [0.05 0.005 0.0005];
        symb_rx = zeros(length(Rx_signal_shift)/M, length(K));
        for index_k = 1:length(K)
            k = K(index_k);
            [symb_rx_k, error_k] = gardnerz(Rx_signal_shift,k,M);
            error(:,index_k,R) = error_k;
            symb_rx(:,index_k) = symb_rx_k;
        end
end
%% preprocess
Merror1=mean(error(:,1,:),3);
Std1=std(error(:,1,:),1,3);
L1 = Merror1+Std1;
L2 = Merror1-Std1;

Merror2=mean(error(:,2,:),3);
Std2=std(error(:,2,:),1,3);
L3 = Merror2+Std2;
L4 = Merror2-Std2;

Merror3=mean(error(:,3,:),3);
Std3=std(error(:,3,:),1,3);
L5 = Merror3+Std3;
L6 = Merror3-Std3;

%% plot
h=[];
p1 =1:100:length(Merror1);
p =1:300:length(Merror1);
p =1:300:length(Merror1);
figure(1)
h(1) = plot(p1,Merror1(p1),'-x');
hold on
plot(p1,L1(p1),'-','color',[0.85 0.85 0.85]);
plot(p1,L2(p1),'-','color',[0.85 0.85 0.85]);

h(2) = plot(p,Merror2(p),'->');
plot(p,L3(p),'--','color',[0.7 0.7 0.7]);
plot(p,L4(p),'--','color',[0.7 0.7 0.7]);

h(3) = plot(p,Merror3(p),'-o');
plot(p,L5(p),'--','color',[0.4 0.4 0.4]);
plot(p,L6(p),'--','color',[0.4 0.4 0.4]);



% plot(p,Merror3(p),'-o','color','k',L5(p),L6(p),'--','color','k');
xlabel('Symbols');
ylabel('Time error(mean\pmstdv)');
title('Convergence as a function of K');
legend(h,{'K=0.05','K=0.005','K=0.0005'});
        

