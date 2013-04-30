function [y] = binconvenc(u,S,O)
%UNTITLED Encodes the input binary vector with the convolutional 
%encoder represented by maps S,O
%   [y] = binconvenc(u,S,O)
    
    mu = length(u);
    %nu = log(length(S(:,1)))/log(2);
    y = zeros(2,mu);
    s = 0;
    
    for i=1:mu
        y(:,i) = de2bi(O(s+1,u(i)+1),2,'left-msb').';
        s = S(s+1,u(i)+1);
    end  

end

