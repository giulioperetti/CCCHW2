function [Y] = nbinconvenc(N,D,u)
%NBINCONVENC Compute the output of a rate 1/n binary convolutional encoder
%  [Y] = convenc(N,D,u)
%   N = nunmerators matrix by rows [z x y] = z+xD+yD^2 
%   D = denunmerators matrix by rows [1 x y] = 1+xD+yD^2 
%   u = input vector 
%   Y = output matrix by rows

n = size(N,1); 
mu = length(u);
Y = gf(zeros(n,mu));

if (n~=size(D,1))
    Y = 0;
else
    for i=1:n
        Y(i,:) = binconvenc(N(i,:),D(i,:),u);
    end
end
end

