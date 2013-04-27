%% 2.1
%Devise a MatLab function which, given the generator polynomials in binary vector form g1 = [101], g2 = [111], 
%builds three maps for: 
%a) the state update function s'+1 = S(s', u'); 
%b) the output function y' = O(s', u');
%c) the neighbors sets N(a)

%%

Alf = [1+1i 1-1i -1+1i -1-1i];

L1 = 0;
L2 = 2;
M = 4; % alfabeto quaternario
N_s = M^(L1+L2); % numero stati

r_DATI = 16; r_TS = 0;

eta = Psi_a(N1_a+D_a+1-L1:N1_a+D_a+1+L2); % risposta complessiva dopo c da -L1 a L2

% Pe_d = zeros(1,10);
% 
% for u=1:10

[a_k,a_n,r_c,w] = tx(r_TS,r_DATI,q_c,Q_0,M_wc); % 2^(12-1) elementi
x = filter(g_mf,1,r_c);

m_max = 30;% ritardo massimo del canale
r_xa=abs(xcorr(x,a_n,m_max,'biased'));
r_xa=r_xa((length(r_xa)+1)/2:end);

t0 = find(r_xa==max(r_xa))-1;

x_T = x(t0+1:Q_0:end);

rho = filter(c_a,1,x_T);
rho = rho(1+D_a:end);

comb = combn(Alf,L1+L2);

% matrice di possibili transizioni, da i posso raggiungere j se p_ij=1 

pp = zeros(N_s,M);

for i=1:N_s
    n=1;
    for j=1:N_s
        if comb(i,2:L1+L2)==comb(j,1:L1+L2-1)
            pp(i,n)=j; % i raggiungibile da j
            n=n+1;
        end
    end
end

K = length(rho);

% seq di training = sigma_1 = comb(1,:)

rho = [comb(1,:) rho];
a_k = [comb(1,:) a_k];

C = inf*ones(N_s,1);
C_pre = C; C_pre(1)=0; % seq training costo zero
L = zeros(N_s,200); % buffer surv seq
L_pre = L;
L_def = zeros(1,K+L1+L2); l=1; % surv seq comuni

for i=1:N_s
    L_pre(i,1:L1+L2)=comb(1,:);
end

tt=L1+L2; % contatore ottimizzazione
for t=1:K
    t/K
    tt=tt+1;
    cost=zeros(N_s,M); % path metric
    for i=1:N_s
        for k=1:M
            % 8.202
            cost(i,k) = C_pre(pp(i,k))+abs(rho(t+L2)-eta*[comb(i,:) comb(pp(i,k),L1+L2)].')^2;
        end
        winner = find(cost(i,:) == min(cost(i,:)));
        % se i vincitori son pi� di uno a parimerito scelgo il primo
        j = winner(1);
        % aggiorno i costi e le survivor seq
        C(i) = cost(i,j);
        L(i,1:tt) = [L_pre(pp(i,j),1:tt-1) comb(i,1)];
    end
    
    % ottimizzazine memoria
    if mod(t,50)==0
        i=0;
        while any((L(:,i+1)==L(1,i+1)*ones(N_s,1))-ones(N_s,1))==0
            i=i+1;
        end
        % ho i colonne uguali
        L_def(l:l+i-1)=L(1,1:i);
        l=l+i; tt=tt-i;
        L = [L(:,i+1:end) zeros(N_s,200-(i+1))];
    end
    L_pre = L;
    C_pre = C;
end

winner = find(C == min(C));
seq_stimata = [L_def(1:l-1) L(winner,1:tt)];