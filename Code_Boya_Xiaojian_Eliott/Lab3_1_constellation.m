clear all; close all;
addpath(genpath('myFunctions'));

% Tx signal parameters
bit_number = 1000;   % length of stream of bits
T = 1e-6;  % symbol interval   
fs = 1/T;
Nbps = 2;  % #Bits per symbol
FC = 2e9; 
df = 10.*(T*FC); 
phi = 3.14;
M = 2;  % up/downsampling factor

%Root raised cosine Filter parameters
beta = 0.3;     % roll-off factor
RRCtaps=101;    % taps of rrc filter (root raised cosine)

    EbN0 = 10;                   % SNR per bit grid 
    if Nbps==1                  
        modulation='pam'; %BPSK
    else
        modulation='qam'; 
    end
    
    %Tx
    bitStreamTx = fix(rand(bit_number,1)*2); %bitStream
    symb_tx = mapping(bitStreamTx,Nbps,modulation);
    symb_tx_upsampled = upsample(symb_tx,M);
    [Hrc,freqGrid,h_rrc,h_rc,timeGrid] = HalfNyquistFilter(beta,T,M,RRCtaps);
    Tx_signal = conv(h_rrc,symb_tx_upsampled);

    %BaseBand equivalent of an ideal channel
    noised_Rx_Signal = addAWGN(Tx_signal,EbN0,bit_number,M/T);

    %without CFO
    Rx_signal0 = conv(h_rrc, noised_Rx_Signal);
    Rx_signal0 = Rx_signal0(RRCtaps:end-(RRCtaps-1)); % usefull part from t=0
    
        % Add CFO
        num = length(noised_Rx_Signal);
        % time vector
        t = (0:length(noised_Rx_Signal)-1)./(fs*M);
        % Add distortion
        temp = exp(1j.* (2*pi*df .* t + phi));
        for i = ((RRCtaps-1)/2+1):num
        noised_Rx_Signal(i)= noised_Rx_Signal(i) * temp(i-(RRCtaps-1)/2);   
        end
        
    %Reciever Side
    Rx_signal = conv(h_rrc, noised_Rx_Signal);
    Rx_signal = Rx_signal(RRCtaps:end-(RRCtaps-1)); % usefull part from t=0
    
    symb_rx = downsample(Rx_signal.',M);
    symb_rx0 = downsample(Rx_signal0,M);
    
    symb_R = real(symb_rx0);
    symb_I = imag(symb_rx0);

    symb_R1 = real(symb_rx);
    symb_I1 = imag(symb_rx);
    
    symb_R2 = real(symb_tx);
    symb_I2 = imag(symb_tx);
    
    figure(1)
    scatter(symb_R1,symb_I1,'ro','MarkerFaceColor','r','LineWidth',1);
    hold on;
    scatter(symb_R,symb_I,'bo','MarkerFaceColor','b','LineWidth',1);
    hold on;
    scatter(symb_R2,symb_I2,'x','y','LineWidth',10);
    legend('After CFO','After AWGN','TX Symbols');
    xlabel('Real part');
    ylabel('Imaginary part');
    xlim([-2 2]);
    ylim([-2 2]);
    grid;
    title('CFO 10ppm, Fc = 2GHz, Eb/N0 = 10dB');