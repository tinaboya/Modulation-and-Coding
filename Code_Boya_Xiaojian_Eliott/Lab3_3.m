clear all; close all;
addpath(genpath('myFunctions'));
%% Tx signal parameters
bit_number = 6*10000;                    % length of stream of bits
T = 1e-6;                                % symbol interval   
fs = 1/T;
Nbps = 2;                                % #Bits per symbol
FC = 2e9; 
df = 30.*(T*FC); 
phi = 3.14;

%% downsampling / upsampling
M = 100;  % Factor
t0 = 10; % time offset sequence
    
    % Root raised cosine Filter parameters
    beta = 0.3;     % roll-off factor
    RRCtaps=8*M+1;     % taps of rrc filter (root raised cosine)

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

        %Reciever Side
        Rx_signal = conv(h_rrc, noised_Rx_Signal);
        Rx_signal = Rx_signal(RRCtaps:end-(RRCtaps-1)); % usefull part from t=0
        
        %time shift
        Rx_signal_new =  [zeros(t0,1); Rx_signal];
        Rx_signal_new = Rx_signal_new(1:length(Rx_signal)); %cutoff or use circshift(s,[]);
        
        
        %gardner
        k=0.01;
        [symb_rx, error] = gardnerz(Rx_signal_new,k,M);
        bitStreamRx = demapping(symb_rx,Nbps,modulation); 

        symb_rx0 = downsample(Rx_signal_new,M); % rx without gardner
        bitStreamRx0 = demapping(symb_rx0,Nbps,modulation); 

        %BER no gardner
        N_errorBit = length(find(bitStreamTx - bitStreamRx0 ~= 0));
        BER0 =  N_errorBit/ length(bitStreamTx);
        BER_p0(index)= BER0;

        %BER
        N_errorBit = length(find(bitStreamTx - bitStreamRx ~= 0));
        BER =  N_errorBit/ length(bitStreamTx);
        BER_p(index)= BER;

    end
%%    % BER figure
    figure(1)
    semilogy(EbN0Grid,BER_p,'->');
    hold on;
    semilogy(EbN0Grid,BER_p0,'-o');
    grid;
    xlabel('Eb/N_0[dB]');
    ylabel('BER');
    title('The impact of the Gardner algorithm ')
    ylim([10^(-4) 1]);
    legend('after Gardner','time shift = 0.1T')
    hold off
  