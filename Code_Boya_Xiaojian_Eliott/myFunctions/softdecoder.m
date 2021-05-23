function bits_i = softdecoder(coded_bits, H, maxit,N0) 
    %maxit: max iterations
    %N0/2:psd of noise
    %% initialize variables
    [Nc, Nv] = size(H);
    r = Nv/Nc;
    bits_l = length(coded_bits);
    bits_i = zeros(1,bits_l/r);  
for k = 1:Nv:bits_l
    block = coded_bits(k:k+Nv-1);       
%     Lq = zeros(size(H));            %initial Lq(check nodes table)
    for j = 1:1:Nc
        for i = 1:1:Nv
            if H(j,i)==1
               Lq(i,j) = -2*block(i)/(N0/2);
            end
        end
    end
    
    blocknew = block;
    iter = 0;
    syndrom = mod(H*blocknew',2);
    
    
    %% iterations
    while  any(syndrom) && (iter < maxit) 
      %% step 1 update variables nodes
      
      Lr = zeros(size(H));
        for j = 1:1:Nc
            R = Lq(find(H(j,:)~=0),j);  %set R: all variable nodes communicate to check nodes i
            lrr = zeros(1,length(R));
            for i = 1:length(R)
                    R_exclude = [R(1:i-1);R(i+1:end)]; %set R_exclude: all variable nodes in R except variable node j
                    Gamma = sign(R_exclude);
                    Alpha = abs(R_exclude);
                    lrr(i) = prod(Gamma)*min(Alpha);
            end
              Lr(j,find(H(j,:)~=0))=lrr;
       end

 
        %% step 2 update check nodes
        
       Lc = -2*blocknew/(N0/2);      %update Lc
       Lq = zeros(size(Lq));
        for i = 1:1:Nv
            Q = Lr(find(H(:,i)~=0),i); %set Q: all check nodes communicate to variable nodes j
            LrSum=zeros(1,length(Q));
            Q_exclude=zeros(length(Q)-1,1);
            for j = 1:length(Q)
                    Q_exclude = [Q(1:j-1);Q(j+1:end)]; %set Q_exclude: all check nodes in Q except check node i
                    LrSum(j)=sum(Q_exclude);
            end
            Lq(i,find(H(:,i)~=0))=(Lc(i)+LrSum).'; 
        end

        %% decision
        LQ=zeros(size(block));
        for i=1:1:Nv
        LQ(i) = Lc(i)+ sum(Lr(find(H(:,i)~=0),i));
%         LQ(i) = Lc(i)+ sum(Lr(:,i));
            if LQ(i)<0
                blocknew(i)=1;
            else
                blocknew(i)=-1;
            end
        end
        iter = iter + 1;
        syndrom = mod(H*blocknew',2);
        
    end
    
    for i=1:1:length(blocknew)
        if blocknew(1,i)==-1
            blocknew(1,i)=0;
        else
            blocknew(1,i)=1;
        end
    end
    
    bits_i(1+(k-1)/r:(k-1)/r+Nc) = blocknew(end-Nc+1:end);
end
end