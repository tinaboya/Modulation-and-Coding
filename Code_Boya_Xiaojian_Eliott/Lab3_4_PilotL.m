clear all; close all;
addpath(genpath('myFunctions'));

%Tx parameter
bit_number = 100000;                   % length of stream of bits
T = 1e-6;                                % symbol interval   
fs = 1/T;
Nbps = 2;                                % #Bits per symbol
FC = 2e9; 
K = 8;
M = 100;  % Factor

% Root raised cosine Filter parameters
    beta = 0.3;     % roll-off factor
    RRCtaps=8*M+1;     % taps of rrc filter (root raised cosine)
    % Initinalization for plot
    EbN0Grid = 0:15;                   % SNR per bit grid 
    ToA = zeros(length(EbN0Grid),3);
    CfoError = zeros(length(EbN0Grid),3);
        if Nbps==1                  
            modulation='pam'; %BPSK
        else
            modulation='qam'; 
        end
 Ngrid = [10,20,40];
 
 for iN = 1:length(Ngrid)
    for index = 1:length(EbN0Grid)
        EbN0_db = EbN0Grid(index);
        Ntest = 50;
        for test=1:Ntest
        % Transmiter Side
        bitStreamTx = fix(rand(bit_number,1)*2);% bitStream
        symb_tx = mapping(bitStreamTx,Nbps,modulation);
        
        % pilots
        pilotsNum = 50;
        pilotL = Ngrid(iN);
        pilotBits = fix(rand(pilotL*Nbps,1)*2); 
        pilotSym = mapping(pilotBits, Nbps, modulation); %modulation
        pilotSymDis = floor(length(symb_tx)/pilotsNum); %length of signals between each pilot
        
        symb_txWithPilots = [];
        for i = 1:pilotsNum
                symb_txWithPilots = [symb_txWithPilots; %previous symb_txWithPilots
                symb_tx((i-1)*pilotSymDis+1:i*pilotSymDis); %the bitstreams between the pilots
                pilotSym %the pilot symbol
                ];
        end
        
        if length(symb_txWithPilots)<length(symb_tx) + pilotL*pilotsNum
            symb_txWithPilots = [symb_txWithPilots;
                symb_tx(i*pilotSymDis+1:end);%add the signals that are not included between two pilot
                ];                              %in another way, the rest of the bitstream are added to the end of last pilot
        end
        
        symb_tx_upsampled = upsample(symb_txWithPilots,M);
        [Hrc,freqGrid,h_rrc,h_rc,timeGrid] = HalfNyquistFilter(beta,T,M,RRCtaps); 
        Tx_signal = conv(h_rrc,symb_tx_upsampled);
        
        noised_Rx_Signal = addAWGN(Tx_signal,EbN0_db,length(symb_txWithPilots)*Nbps ,M/T);
                        
    %Reciever Side
    Rx_signal = conv(h_rrc, noised_Rx_Signal);
    Rx_signal = Rx_signal(RRCtaps:end-(RRCtaps-1)); % usefull part from t=0
    symb_rx = downsample(Rx_signal,M); 
        
        block = pilotL + pilotSymDis; %the length of each block
        Nt= []; PPMt = [];
        for ip = 1:floor((length(symb_rx)-pilotL)/block)
            ninit = pilotL + 1 + (ip-1)*block;
            nfinal = pilotL + (ip)*block;
            rx = symb_rx(ninit:nfinal);
            nt = pilotToAesti(length(rx)-pilotL,rx,pilotSym, K);
            Nt(ip) = nt; %save the value of ntilda to NTilda
            
            CFOestimated = - CFOestimate(rx, pilotSym,n_tilda,K,T);
            PPMt(ip) = CFOestimated/FC/T; %save the value of CFOestimated to PPMTilda
        end
        
        %save ToA 
        ToA(index,iN) = ToA(index,iN) + (1/Ntest)*std(Nt);
        CfoError(index,iN) = CfoError(index,iN) + (1/Ntest)*std(PPMt);
        end
    end
 end  
 
% plot
figure(1)
for i = 1:3
plot(EbN0Grid,ToA(:,i));
hold on
end
grid;
legend('N = 10, K =8','N = 20, K =8','N = 40, K =8');
xlabel('EbN0');
ylabel('Time Error stdev');
title('CFO - Random Qpsk pilot symbols, no CFO, no time shift');

figure(2)
for i = 1:3
plot(EbN0Grid,CfoError(:,i));
hold on
end
grid;
legend('N = 10, K =8','N = 20, K =8','N = 40, K =8');
xlabel('EbN0');
ylabel('Frequency Error stdev[ppm]');
ylim([0,10]);
title('CFO - Random Qpsk pilot symbols, no CFO, no time shift');