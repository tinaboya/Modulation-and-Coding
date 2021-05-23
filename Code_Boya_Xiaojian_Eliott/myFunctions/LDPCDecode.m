function bits_i = LDPCDecode(coded_bits, H, maxit) 

% parameters
[R_H, C_H] = size(H);
r = C_H/R_H;
bits_l = length(coded_bits);
bits_i = zeros(1, bits_l/r);
% bloc by bloc
for k = 1:C_H:bits_l
    block = coded_bits(k:k+C_H-1);
    newblock = block;
    iter = 0;
    syndrom = mod(H*newblock',2);
    % Iterate on tanner graph
    while  any(syndrom) && (iter < maxit)
        iter = iter + 1;
        newblock = (sum(xor(syndrom,newblock).*H) + block) > sum(H)/2; 
        syndrom = mod(H*newblock',2);
    end
    % bits with information
    bits_i(1+(k-1)/r:(k-1)/r+R_H) = newblock(end-R_H+1:end);
end
