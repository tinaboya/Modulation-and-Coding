clear all; 
%close all; clc;
addpath(genpath('libs'));
addpath(genpath('myFunctions'));
%%
%clear all; close all;
addpath(genpath('myFunctions'));
%% Tx signal parameters

bit_number = 2*10000;              % length of stream of bits
bps = 2e6;
Nbps = 2;                                % #Bits per symbol
FC = 2e9; 
T = 1/(bps/Nbps);                   % symbol interval 1e-6;         
fs = 1/T;                           %symbol frequency
k = 0.0005;
df = 0.*(T*FC);                    %CFO
phi = 3.14;

%% downsampling / upsampling factors
M = 100;  % Up/Down sampling Factor
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
error = zeros(N,Realiations); %for mean the error
symb_rx = zeros(size(error));   
%% Transmitter Side

 for R=1:Realiations            % implement 10 times realizations to get the mean
        bitStreamTx = fix(rand(bit_number,1)*2);% bitStream
        symb_tx = mapping(bitStreamTx,Nbps,modulation);      
        symb_tx_upsampled = upsample(symb_tx,M);
        [H_RC,freqGrid,h_rrc,h_rc,timeGrid] = HalfNyquistFilter(beta,T,M,RRCtaps);
        Tx_signal = conv(h_rrc,symb_tx_upsampled);

%% Reciever Side
        %BaseBand equivalent of an ideal channel
        %EbN0_db = 40;
        %noised_Rx_Signal = addAWGN(Tx_signal,EbN0_db,bit_number,M/T);
        noised_Rx_Signal = Tx_signal;  %no noise
        
%% Add CFO
       
        num = length(noised_Rx_Signal);
        % time vector
        t = (0:length(noised_Rx_Signal)-1)./(fs*M);
        % Add distortion
        temp = exp(1j.* (2*pi*df .* t + phi));
        
             for i = ((RRCtaps-1)/2+1):num
             CFO_Rx(i)= noised_Rx_Signal(i) * temp(i-(RRCtaps-1)/2);   
             end
             
%% RX side RRC
        Rx_signal = conv(h_rrc, CFO_Rx);
        Rx_signal = Rx_signal(RRCtaps:end-(RRCtaps-1)); % usefull part from t=0
        
%% Compensate CFO
%         num = length(Rx_signal);
%         % time vector
%         t = (0:length(Rx_signal)-1)./fs;
%         % Add compemsation
%         temp = exp(-1j.* (2*pi*df .* t + phi));
%             for i = 1:num
%             Rx_signal(i)= Rx_signal(i) * temp(i);
%             end
%% time shift
        Rx_signal_shift = circshift(Rx_signal,t0);
%% gardener
        
        [symb_rx(:,R), error(:,R)] = gardnerz(Rx_signal_shift,k,M);        
 end

%% plot

Merror=mean(error,2);
Std=std(error,1,2);
L1 = Merror+Std;
L2 = Merror-Std;
p =1:200:length(Merror);
figure(1)
plot(p,Merror(p),'-o','color','k');
xlabel('Symbols');
ylabel('Time error');
title('Gardener Robustness to CFO');
hold on
plot(p,L1(p),'--','color','k');
plot(p,L2(p),'--','color','k');
%%
legend('CFO=10ppm','','','CFO=50ppm','','','no CFO','','');
ylabel('Time error(mean\pmstdv)');
%%
%  figure(2)
%          plot(error(:,1));
%         hold on
%         plot(error(:,2));
%         plot(error(:,3));      
% xlabel('Symbols');
% ylabel('Time error');
% title('robustness');
% legend('no CFO','CFO=10ppm','CFO=50ppm');
       
