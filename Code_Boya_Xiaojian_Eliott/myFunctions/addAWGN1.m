function [noisedSignal,N0]=addAWGN1(signalTx,EbN0,NBits,fs)
% INPUTS:
%   signalTx    : add noise to this signal
%   EbN0        : SNR
%   NBits       : number of bits
%   fs          : sampling frequency

EbN0 = 10^(EbN0/10);
SignalEnergy = (trapz(abs(signalTx).^2)) /fs; 
Eb = SignalEnergy/NBits;        
Eb = Eb/2;
N0 = Eb/EbN0;       
NoisePower = 2*N0*fs;

% noise
if isreal(signalTx)  % if pam modulation, only real part
    % Half power in the real part, half power in the imaginary part
    noise = sqrt(NoisePower/2)*(randn(length(signalTx),1)); 
else %if qam, with real part and im part
    noise = sqrt(NoisePower/2)*(randn(length(signalTx),1) + 1i*randn(length(signalTx),1));
end

noisedSignal = signalTx + noise;
end