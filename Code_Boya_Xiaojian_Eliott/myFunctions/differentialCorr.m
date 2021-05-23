function x =  differentialCorr(y,a,n,k)
    N = length(a);
    sum = 0;
    for l=k:N-1
        sum = sum + (conj(y(n+l+1))*a(l+1))*(conj(conj(y(n+l-k+1))*a(l-k+1)));
    end
    x = (1/(N-k))*sum;
end