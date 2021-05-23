function [Hrc,freqGrid,h_rrc,h_rc,timeGrid]=HalfNyquistFilter(beta,T,M,RRCtaps)
% INPUTS:
%   beta    : Roll-off factor.
%   T       : Symbol interval.
%   M       : upsampling factor.
%   RRCtaps : numbers of taps

fs = M*1/T; % sample frequency
Hrc = zeros(RRCtaps,1) ;     

stepOffset = (1/RRCtaps)*fs;
highestFreq = stepOffset*(RRCtaps-1)/2; 
freqGrid = linspace(-highestFreq,highestFreq,RRCtaps);
timeGrid = (-(RRCtaps-1)/2 : (RRCtaps-1)/2)*(1/fs);

%Generation of H_RC
for i=1:1:RRCtaps
     if abs(freqGrid(i)) <= (1-beta)/(2*T) 
         Hrc(i) = T;
     elseif (abs(freqGrid(i)) > (1-beta)/(2*T)) && (abs(freqGrid(i))<=(1+beta)/(2*T))
         Hrc(i) = (T/2)*(1+cos((pi*T/beta)*(abs(freqGrid(i))-(1-beta)/(2*T))));
     else
         Hrc(i) = 0;
     end
end

%Normalization of Hrc
hrc = fftshift(ifft(ifftshift(Hrc),'symmetric'));
hrc0 = hrc(((RRCtaps-1)/2)+1);
HrcNormalized = Hrc/hrc0;

h_rc =fftshift(ifft(ifftshift(HrcNormalized),'symmetric'));

%Hrrc
Hrrc = sqrt(HrcNormalized);  
h_rrc = ifft(ifftshift(Hrrc),'symmetric');
h_rrc = fftshift(h_rrc);
end