function [ynew, error] = gardnerz(y,k,M)
%y:recieved symbols
%k:error weight
%M:upsampling factor
error = zeros(length(y)/M,1);
error(2) =0.2;
T = 1e-6;
ynew = zeros(size(error));
ynew(1) = y(1);
for n=2:length(error)
    x= (0:M).*T;
    v= y(M*(n-2)+1:M*(n-1)+1);
    xq= [M/2-error(n)*M, M-error(n)*M].*T;
    ymid = interp1(x,v,xq,'linear','extrap');   %interpolation
    ynew(n) = ymid(2);
    error(n+1) = error(n) + 2*k*real(ymid(1)*(ynew(n)'- ynew(n-1)'));
end
 error=error(2:end).*T;%remove the extra 

