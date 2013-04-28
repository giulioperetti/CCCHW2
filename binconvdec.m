function [u] = binconvdec(r,S,O,N)
%BINCONVDEC Viterbi decoder
%   Detailed explanation goes here

    mu = length(r);
    nu = log(length(S(:,1)))/log(2);
    
    C = zeros(2^nu,1); % path cost vector
    Caux = C; % auxiliary cost vector
    
    % initial condition
    C = C - 10^3; 
    C(1) = 0;
    
    U = zeros(2^nu,mu); % survivors matrix
    
    for i=1:mu
        
        Uaux = U;
        
        for j=1:2^nu
            
            % analyzing state j-1
            
            % predecessors
            pred = N(j,:); 
            
            % transitions from predecessors
            trans = [find(S(pred(1)+1,:)==j-1)-1 
                find(S(pred(2)+1,:)==j-1)-1];
            
            % cost function
            temp = [C(pred(1)+1) + r(i)*pamap(O(pred(1)+1,trans(1)+1)) 
                C(pred(2)+1) + r(i)*pamap(O(pred(2)+1,trans(2)+1))];
            
            % select the max value
            [mcost ind] = max(temp);
            
            % path and cost update
            
            Uaux(j,1:i) = U(pred(ind)+1,1:i);
            Uaux(j,i) = trans(ind); 
            Caux(j) = mcost;
                        
        end
        
        C = Caux - max(Caux);
        U = Uaux;
    end
    
    winner = find(C == max(C));
    u = U(winner,:);
    
end

