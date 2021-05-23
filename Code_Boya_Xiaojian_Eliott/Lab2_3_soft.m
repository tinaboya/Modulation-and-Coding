clear all; close all;
addpath(genpath('myFunctions'));

% Tx signal parameters
Nbps = 1; % #Bits per symbol only BPSK here
K = 128;  % length of the block
N = 2*K;  % length after encoding
T = 1e-6; % symbol interval   
M = 3;  % up and down sampling factor
modulation='pam'; %BPSK
% Root raised cosine filter parameters 
beta = 0.3;     % roll-off factor
RRCtaps = 101;     % taps of rrc filter (root raised cosine)

% LDPC parameters
H = makeLdpc(K,N,0,1,3);    % Matrix Parity check H

EbN0 = -6:0.5:8;                 % SNR per bit 
Maxiterations = [1 2 3 4 5];

%% Transmitter 
NBits =500*K*Nbps;    % a multiple of Nbps to map successfully; a multiple of K to be divided into blocks
 

for index = 1:length(EbN0)
        EbN0_db = EbN0(index);
        
        % Tx
        bits = rand(NBits,1)>0.5;
        
        % LDPC encoding
        bitsTx_coded = [];
        for i=1:(length(bits)/K)
            block = bits(1 + (i-1)*K : i*K);
            [checkBits,H] = makeParityChk(block,H,0); % encode function
            bitsTx_coded = [bitsTx_coded; checkBits; block;];
        end
        
        symb_tx_coded = mapping(bitsTx_coded,Nbps,modulation);
        symb_tx_upsampled = upsample(symb_tx_coded,M);
        [Hrc,freqGrid,h_rrc,h_rc,timeGrid] = HalfNyquistFilter(beta,T,M,RRCtaps);
        Tx_signal_coded = conv(h_rrc,symb_tx_upsampled);
        
        % add noise
        [noised_Rx_Signal,N0] = addAWGN1(Tx_signal_coded,EbN0_db,NBits*size(H,2)/size(H,1),M/T);
        
        % Rx
        Rx_signal = conv(h_rrc,noised_Rx_Signal);
        symb_rx_upsampled = Rx_signal(RRCtaps:end-(RRCtaps-1)); % useful part from t=0
        symb_rx = downsample(symb_rx_upsampled,M);
        bitsRx = demapping(symb_rx,Nbps,modulation);
    for indexIt=1:length(Maxiterations)
        MIt = Maxiterations(indexIt);
        % LDPC decoding
        [bitsRx_decoded_soft]=softdecoder(symb_rx',H,MIt,N0);   
        [bitsRx_decoded_hard]=LDPCDecode(bitsRx',H,5);

        bitsRx_decoded_soft = bitsRx_decoded_soft(:);
        bitsRx_decoded_hard = bitsRx_decoded_hard(:);
        
        N_errorBit = length(find(bitsTx_coded - bitsRx ~= 0));
        BER_uncoded =  N_errorBit/ length(bitsTx_coded);

        N_errorBit_ldpc = length(find(bits - bitsRx_decoded_soft ~= 0));%soft
        BER_encoded =  N_errorBit_ldpc/ length(bits);
        
        N_errorBit_hard = length(find(bits - bitsRx_decoded_hard ~= 0));%hard
        BER_encoded_hard =  N_errorBit_hard/ length(bits);

        BER_uncoded_matrix(index) = BER_uncoded;
        BER_hard_decoded_matrix(index) = BER_encoded_hard;
        BER_encoded_matrix(indexIt,index) = BER_encoded;
    end 
end
%%
figure(1)

for i=[1 2 3 4]
    semilogy(EbN0,BER_encoded_matrix(i,:),'-s','MarkerSize',5);
    hold on;
end
hold on;
semilogy(EbN0,BER_hard_decoded_matrix,'-o','MarkerSize',5);
semilogy(EbN0,BER_uncoded_matrix,'->','MarkerSize',5);
% legend('encoded','uncoded');
legend('Soft decoded MaxIt=1','Soft decoded MaxIt=2','Soft decoded MaxIt=3','Soft decoded MaxIt=4','Hard decoded Max iteration=5','Uncoded');
xlim([0 8]);
ylim([2e-4 1])
grid on;
xlabel('E_b/N_0(dB)');
ylabel('BER');
title('Soft decoded, Hard decoded and Uncoded BER (BPSK)')
%%
% t=1:1:50;
% y =sin(t);
% figure
% plot(y,'o','MarkerSize',9,'MarkerFaceColor','auto');