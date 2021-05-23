clear all; close all; 
addpath(genpath('myFunctions'));

% Tx signal parameters
Nbps_vec = [1, 2, 4, 8]; % #Bits per symbol
K = 128;  % length of the block
N = 2*K;  % length after encoding
T = 1e-6; % bit interval   
M = 3;  % up and down sampling factor

% Root raised cosine filter parameters 
beta = 0.3;     % roll-off factor
RRCtaps = 53;     % taps of rrc filter (root raised cosine)

% LDPC parameters
H = makeLdpc(K,N,0,1,3);    % Matrix Parity check H

EbN0 = 0:10;                 % SNR per bit 
Maxiterations = [1 5 10];
for Nbps_loopindex=1:4
    Nbps = Nbps_vec(Nbps_loopindex);
    NBits = 1000*K*Nbps;    % a multiple of Nbps to map successfully
                          % a multiple of K to be divided into blocks
    if Nbps==1                  
        modulation='pam'; %BPSK
    else
        modulation='qam'; 
    end
for indexIt=1:length(Maxiterations)
    MIt = Maxiterations(indexIt);
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
        noised_Rx_Signal = addAWGN(Tx_signal_coded,EbN0_db,NBits*size(H,2)/size(H,1),M/T);
        
        % Rx
        Rx_signal = conv(h_rrc,noised_Rx_Signal);
        symb_rx_upsampled = Rx_signal(RRCtaps:end-(RRCtaps-1)); % useful part from t=0
        symb_rx = downsample(symb_rx_upsampled,M);
        bitsRx = demapping(symb_rx,Nbps,modulation);
        
        % LDPC decoding
        [bitsRx_decoded]=LDPCDecode(bitsRx',H,MIt);
        
        bitsRx_decoded = bitsRx_decoded(:);
        
        N_errorBit = length(find(bitsTx_coded - bitsRx ~= 0));
        BER_uncoded =  N_errorBit/ length(bitsTx_coded);

        N_errorBit_ldpc = length(find(bits - bitsRx_decoded ~= 0));
        BER_encoded =  N_errorBit_ldpc/ length(bits);

        BER_uncoded_matrix(index) = BER_uncoded;
        BER_encoded_matrix(indexIt,index) = BER_encoded;
    end 
end

figure(1)
subplot(2,2,Nbps_loopindex);
for i=1:length(Maxiterations)
    semilogy(EbN0,BER_encoded_matrix(i,:));
    hold on;
end
hold on;
semilogy(EbN0,BER_uncoded_matrix);
legend('MaxIt=1','MaxIt=5','MaxIt=10','Uncoded');
xlim([0 6]);
grid on;
xlabel('SNR per bit(dB)');
ylabel('BER');
if Nbps_loopindex == 1
title('Max iterations gain (BPSK)')
else if Nbps_loopindex == 2
        title('Max iterations gain (4QAM)');
    else if Nbps_loopindex == 3
            title('Max iterations gain (16QAM)');
        else
            title('Max iterations gain (64QAM)');
        end
    end
end
end
