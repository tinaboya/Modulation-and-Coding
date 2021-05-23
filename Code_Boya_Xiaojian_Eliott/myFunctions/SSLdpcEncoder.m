function [encodedBits,H]=SSLdpcEncoder(inBits,H)
    % put inBits to vector.
    inBits = inBits';
    inBits = inBits(:); 
    [K,N] = size(H);    
    % H = [I P']
    H=mod(rref(H),2);     
    % G2 = [P I];
    G = gen2par(H);       
    encodedBits = zeros((length(inBits)/K),N);
    % encode bits
    for i=1:(length(inBits)/K)
        block = inBits(1+(i-1)*K:i*K);
        encodedBits(i,:) = mod(block'*G,2);
    end
end

