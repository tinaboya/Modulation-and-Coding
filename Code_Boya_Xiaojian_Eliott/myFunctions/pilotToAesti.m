function n_esti  = pilotToAesti(nl,y,pilot_symbol,K)

    maxSum= -Inf;
    n_esti = 0;
    for n = 1:nl
        sum = 0;
        for k=1:K
            sum = sum + abs(differentialCorr(y,pilot_symbol,n,k));
        end
        if maxSum < sum  %argmax
            n_esti = n;
            maxSum = sum;
        end        
    end
    
end
