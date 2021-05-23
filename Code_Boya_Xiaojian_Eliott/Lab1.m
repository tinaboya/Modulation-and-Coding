clear all; close all;
addpath(genpath('myFunctions'));

%Tx signal parameters
NBits = 6*1e7;                          % length of stream of bits
T = 1e-6;                                % bit interval   
M = 2;                                   % up/down sampling factor

%Root raised cosine Filter parameters
beta = 0.3;     % roll-off factor
RRCtaps=101;    % taps of rrc filter (root raised cosine)

%Tx and Rx
EbN0Grid = 0:25;
error = zeros(4,length(EbN0Grid)); 
Nbps_vec = [1, 2, 4, 8]; % different modulation

for Nbps_loopindex=1:4
    Nbps = Nbps_vec(Nbps_loopindex);
    if Nbps==1                  
        modulation='pam'; %BPSK
    else
        modulation='qam'; 
    end
for index = 1:length(EbN0Grid)
    EbN0 = EbN0Grid(index);
    
    %Tx
    bits = rand(NBits,1)>0.5; % bitStream
    txSymb = mapping(bits,Nbps,modulation);
    upsampledTxSymb = upsample(txSymb,M);
    [Hrc,freqGrid,hrrc,hrc,timeGrid] = HalfNyquistFilter(beta,T,M,RRCtaps);
    TxSignal = conv(hrrc,upsampledTxSymb);
    
    %add noise
    noisedRxSignal = addAWGN(TxSignal,EbN0,NBits,M/T); %add noise

    %Rx
    rxSignal = conv(hrrc,noisedRxSignal);
    upsampledSymbRx = rxSignal(RRCtaps:end-(RRCtaps-1)); % usefull part from t=0
    rxSymb = downsample(upsampledSymbRx,M);
    bitsRx = demapping(rxSymb,Nbps,modulation); 

    %BER
    errorBitN = length(find(bits - bitsRx ~= 0));
    BER =  errorBitN/ length(bits);
  
    error(Nbps_loopindex,index) = BER;
end
end

%limited communication bandwidth
figureA = figure;
plot(freqGrid/1e6,Hrc);
xlabel('Frequency[MHz]');
title('HRC Spectrum')

%cancellation of ISI
figureB = figure;
plot(timeGrid*1e6, hrc); hold on;
plot(timeGrid(1:M:RRCtaps)*1e6, hrc(1:M:RRCtaps),'o');
legend('RC filter impulse response','Sampling at symbol frequency');
grid;
xlabel('time[us]');
ylabel('Normalized amplitude');
title('Illustration of the cancellation of ISI by the Nyquist root raised cosine filter');
set(figureB,'position',[400,0,400,400]);

% BER figure
figureC = figure;
for i=1:4
    semilogy(EbN0Grid,error(i,:));
    hold on;
end
legend('BPSK','4QAM','16QAM','64QAM');
grid;
xlabel('Eb/N0');
ylabel('BER');
title('BER for different modulations')
ylim([10^(-4) 1]);