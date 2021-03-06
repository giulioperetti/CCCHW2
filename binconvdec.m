function [u] = binconvdec(r,S,O,N)

    mu = length(r(1,:));
    nu = log(length(S(:,1)))/log(2);
    
    C = zeros(2^nu,1); % path cost vector
    Caux = C; % auxiliary cost vector
    
    % initial condition
    C = C - 10^5; 
    C(1) = 0;
    
    U = zeros(2^nu,mu); % survivors matrix
    Uaux = zeros(2^nu,10*nu); % auxiliary matrix
    
    for i=1:mu    
        
        for j=1:2^nu
            
            % analyzing state j-1
            
            % predecessors
            pred = N(j,:); 
            
            % transitions from predecessors
            trans = floor((j-1)/2);
            
            % cost function
            temp = [C(pred(1)+1) + r(:,i).'*(pamap(de2bi(O(pred(1)+1,trans+1),2,'left-msb')).') 
                C(pred(2)+1) + r(:,i).'*(pamap(de2bi(O(pred(2)+1,trans+1),2,'left-msb')).')];
            
            % select the max value
            [mcost ind] = max(temp);
            
            % path and cost update
            Caux(j) = mcost;
            Uaux(j,1:min(i,10*nu)) = U(pred(ind)+1,max(i-10*nu+1,1):i);
            Uaux(j,min(i,10*nu)) = trans;

        end
        
        C = Caux - max(Caux);
        U(:,max(i-10*nu+1,1):i) = Uaux(:,1:min(i,10*nu));
        
    end
    
    u = U(find(C == max(C)),:);
    
end

