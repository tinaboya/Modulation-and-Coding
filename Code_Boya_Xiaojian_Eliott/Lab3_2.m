clear all; close all;
addpath(genpath('myFunctions'));
%% Tx signal parameters
bit_number = 6*100000;                   % length of stream of bits
T = 1e-6;                                % symbol interval   
fs = 1/T;
Nbps = 2;                                % #Bits per symbol
FC = 2e9; 
ppm = 30;
df = ppm.*(T*FC); 
phi = 3.14;

%% downsampling / upsampling
M = 20;  % Factor

for t0 = 0:1:2 % time offset sequence
    
    % Root raised cosine Filter parameters
    beta = 0.3;     % roll-off factor
    RRCtaps=101;     % taps of rrc filter (root raised cosine)

    % Initinalization for plot
    EbN0Grid = -4:1:10;                   % SNR per bit grid 
        if Nbps==1                  
            modulation='pam'; %BPSK
        else
            modulation='qam'; 
        end
    for index = 1:length(EbN0Grid)
        EbN0_db = EbN0Grid(index);
        %Transmiter Side
        bitStreamTx = fix(rand(bit_number,1)*2);% bitStream
        symb_tx = mapping(bitStreamTx,Nbps,modulation);
        symb_tx_upsampled = upsample(symb_tx,M);
        [H_RC,freqGrid,h_rrc,h_rc,timeGrid] = HalfNyquistFilter(beta,T,M,RRCtaps);
        Tx_signal = conv(h_rrc,symb_tx_upsampled);

        %BaseBand equivalent of an ideal channel
        noised_Rx_Signal = addAWGN(Tx_signal,EbN0_db,bit_number,M/T);

        %without CFO
        Rx_signal0 = conv(h_rrc, noised_Rx_Signal);
        Rx_signal0 = Rx_signal0(RRCtaps:end-(RRCtaps-1)); % usefull part from t=0

    %         % Add CFO
    %         num = length(noised_Rx_Signal);
    %         % time vector
    %         t = (0:length(noised_Rx_Signal)-1)./fs;
    %         % Add distortion
    %         temp = exp(1j.* (2*pi*df .* t + phi));
    %         for i = ((RRCtaps-1)/2+1):num
    %         noised_Rx_Signal(i)= noised_Rx_Signal(i) * temp(i-(RRCtaps-1)/2);   
    %         end

        %Reciever Side
        Rx_signal = conv(h_rrc, noised_Rx_Signal);
        Rx_signal = Rx_signal(RRCtaps:end-(RRCtaps-1)); % usefull part from t=0
        Rx_signal_new =  [zeros(t0,1); Rx_signal];
        Rx_signal_new = Rx_signal_new(1:length(Rx_signal));
       
            
        %compensate CFO
    %     num = length(Rx_signal);
    %         % time vector
    %         t = (0:length(Rx_signal)-1)./fs;
    %         % Add distortion
    %         temp = exp(-1j.* (2*pi*df .* t + phi));
    %         for i = 1:num
    %         compRx_signal(i)= Rx_signal(i) * temp(i);
    %         end
    %     
    %     symb_rx = downsample(compRx_signal.',M);
        symb_rx = downsample(Rx_signal_new .',M);
        bitStreamRx = demapping(symb_rx.',Nbps,modulation); 

%         symb_rx0 = downsample(Rx_signal0,M);
%         bitStreamRx0 = demapping(symb_rx0,Nbps,modulation); 

%         %BER
%         N_errorBit = length(find(bitStreamTx - bitStreamRx0 ~= 0));
%         BER0 =  N_errorBit/ length(bitStreamTx);
%         BER_p0(index)= BER0;

        %BER
        N_errorBit = length(find(bitStreamTx - bitStreamRx ~= 0));
        BER =  N_errorBit/ length(bitStreamTx);
        BER_p(index)= BER;

    end
    % BER figure
    figure(1)
    semilogy(EbN0Grid,BER_p,'->');
    hold on;
end
grid;
xlabel('Eb/N_0[dB]');
ylabel('BER');
title('Impact of time shift on bit error rate')
ylim([10^(-4) 1]);
legend('no time offset','t_0=0.05 T','t_0=0.1 T')
hold off
