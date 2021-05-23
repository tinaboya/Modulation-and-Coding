clear all; close all;
addpath(genpath('myFunctions'));

% Tx signal parameters
bit_number = 100000;                   % length of stream of bits
T = 1e-6;                                % symbol interval   
fs = 1/T;
Nbps = 2;                                % #Bits per symbol
FC = 2e9; 
ppm_vec = [0 10];
t0_vec = [0 2];
phi = 3.14;
K = 8;
M = 100;  % Factor
    
    % Root raised cosine Filter parameters
    beta = 0.3;     % roll-off factor
    RRCtaps = 8*M+1;     % taps of rrc filter (root raised cosine)
    % Initinalization for plot
    EbN0Grid = 0:15;                   % SNR per bit grid 
    ToA = zeros(length(EbN0Grid),2);
    CfoError = zeros(length(EbN0Grid),2);
        if Nbps==1                  
            modulation='pam'; %BPSK
        else
            modulation='qam'; 
        end
                
    for ppmi = 1:length(ppm_vec)
            ppm = ppm_vec(ppmi);
            t0 = t0_vec(ppmi);
            df = ppm *(T*FC); 
    for index = 1:length(EbN0Grid)
        EbN0_db = EbN0Grid(index);
        Ntest = 30;
        for test=1:Ntest
        % Transmiter Side
        bitStreamTx = fix(rand(bit_number,1)*2);% bitStream
        symb_tx = mapping(bitStreamTx,Nbps,modulation);
        
        % pilots
        pilotsNum = 50;
        pilotL = 40;
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
        
        %compensate CFO
        num = length(Rx_signal);
        % time vector
        t = (0:length(Rx_signal)-1)./(fs*M);
        % Add distortion
        temp = exp(-1j.* (2*pi*df .* t + phi));
        for i = 1:num
        Rx_signal(i)= Rx_signal(i) * temp(i);
        end
        
        Rx_signal = circshift(Rx_signal.',t0);
        symb_rx = downsample(Rx_signal,M); 

        block = pilotL + pilotSymDis; %the length of each block
        NT= []; PPMT = [];
        for ip = 1:floor((length(symb_rx)-pilotL)/block)
            ninit = pilotL + 1 + (ip-1)*block;
            nfinal = pilotL + (ip)*block;
            rx = symb_rx(ninit:nfinal);
            nt = pilotToAesti(length(rx)-pilotL,rx,pilotSym, K);
            NT(ip) = nt; %save the value of ntilda to NTilda
            
            CFOestimated = - CFOestimate(rx, pilotSym,nt,K,T);
            PPMT(ip) = CFOestimated/FC/T; %save the value of CFOestimated to PPMTilda
        end
        
            %save ToA 
            ToA(index,ppmi) = ToA(index,ppmi) + (1/Ntest)*std(NT);
            CfoError(index,ppmi) = CfoError(index,ppmi) + (1/Ntest)*std(PPMT);
        end
        end
    end
        % plot
        figure(1)
        plot(EbN0Grid,ToA(:,1),'-o','MarkerIndices',1:2:length(EbN0Grid));
        hold on
        plot(EbN0Grid,ToA(:,2),'-x','MarkerIndices',1:2:length(EbN0Grid));
        grid;
        legend('no CF0, t=0','10 ppm CFO, t=0.02');
        %ylim([0 1]);
        xlabel('EbN0');
        ylabel('Time Error stdev');
        title('Random QPSK pilot symbols, N=40, K=8');
     
        figure(2)
        plot(EbN0Grid,CfoError(:,1),'-o','MarkerIndices',1:2:length(EbN0Grid));
        hold on
        plot(EbN0Grid,CfoError(:,2),'-x','MarkerIndices',1:2:length(EbN0Grid));
        grid;
        legend('no CF0, t=0','10 ppm CFO, t=0.02');
        ylim([0 1]);
        xlabel('EbN0');
        ylabel('Frequency Error stdev');
        title('Random QPSK pilot symbols, N=40, K=8');