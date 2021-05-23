function [x] = CFOestimate(yn,a,nEst,K,T)
x = 0;
for k=1:K
    x = x - 1/(K*2*pi*k*T)*angle(differentialCorr(yn,a,nEst,k));
end
end
