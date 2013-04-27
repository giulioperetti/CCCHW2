function [S,O,N] = binconvmaps(n)
%BINCONVENC Build a binary convolution encoder from the generator polynomials in binary vector form.
%[S,O,N] = binconvec(n) 
% S(s,u) is the state update function matrix
% O(s,u) is the output function matrix
% N(s) is the neghbour vector
    
    n = gf(n);
    nu = length(n)-1; % memory of the system
    %s = gf(zeros(1,nu)); % state of the finite state machine
    
    O = gf(zeros(2^nu,2)); % initialization of the output map
    for i=1:2^nu
        O(i,1) = n*gf([0 de2bi(i-1,nu,'left-msb')]).';
        O(i,2) = n*gf([1 de2bi(i-1,nu,'left-msb')]).';
    end
    
    S = zeros(2^nu,2); % initialization of the state update map
    for i=1:2^nu
        temp = de2bi(i-1,nu,'left-msb');
        S(i,1) = bi2de([0 temp(1:end-1)],'left-msb');
        S(i,2) = bi2de([1 temp(1:end-1)],'left-msb');
    end
   
    N = zeros(2^nu,2); % initialization of the state update map
    for i=1:2^nu
        [temp,temp2] = find(S==i-1);
        N(i,:) = bi2de(de2bi(temp-1,2,'left-msb'),'left-msb');   
    end


end

