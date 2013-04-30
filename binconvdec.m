function [u] = binconvdec(r,S,O,N)
%BINCONVDEC Viterbi decoder
%   Detailed explanation goes here

    mu = length(r(1,:));
    nu = log(length(S(:,1)))/log(2);
    
    C = zeros(2^nu,1); % path cost vector
    Caux = C; % auxiliary cost vector
    
    % initial condition
    C = C - 10^3; 
    C(1) = 0;
    
    U = zeros(2^nu,mu); % survivors matrix
    Uaux = zeros(2^nu,5*nu); % auxiliary matrix
    
    for i=1:mu    
        
        for j=1:2^nu
            
            % analyzing state j-1
            
            % predecessors
            pred = N(j,:); 
            
            % transitions from predecessors
            trans = [find(S(pred(1)+1,:)==j-1)-1 
                find(S(pred(2)+1,:)==j-1)-1];
            
            % cost function
            temp = [C(pred(1)+1) + r(:,i).'*pamap(de2bi(O(pred(1)+1,trans(1)+1),2,'left-msb')).' 
                C(pred(2)+1) + r(:,i).'*pamap(de2bi(O(pred(2)+1,trans(2)+1),2,'left-msb')).'];
            
            % select the max value
            [mcost ind] = max(temp);
            
            % path and cost update
            Caux(j) = mcost;
            Uaux(j,1:min(i,5*nu)) = U(pred(ind)+1,max(i-5*nu+1,1):i);
            Uaux(j,min(i,5*nu)) = trans(ind);

        end
        
        C = Caux - max(Caux);
        U(:,max(i-5*nu+1,1):i) = Uaux(:,1:min(i,5*nu));
        
    end
    
    u = U(find(C == max(C)),:);
    
end

