% HW2

clear;

% binary code
g1 = [1 0 1];
g2 = [1 1 1];

% input
ord = 6;
mu = 10^ord; 
u = round(rand(1,mu));

% maps
[SS,OO,NN] = binconvmaps(g1,g2);

% encoding
y = cbinconvenc(u,SS,OO);

% modulation
s = pamap(y);

% channel simulation
Eb_N0_min = 0;
Eb_N0_max = 11;
Eb_N0_dB = Eb_N0_min:Eb_N0_max;

errors = zeros(1,Eb_N0_max-Eb_N0_min);

x = 50; 
e = zeros(1,x);

for i=1:Eb_N0_max-Eb_N0_min+1
    
    for j=1:x
        
        r = awgn(s,Eb_N0_dB(i),'measured');
    
        % viterbi decoding
        uhat = cbinconvdec(r.',SS.',OO.',NN.');

        e(j) = sum(u~=uhat);
        
    end
    errors(i) = sum(e)/(mu*x);

    fprintf('%d % \n',i);
end
%%
Pb_PAM = qfunc(sqrt(2*10.^(Eb_N0_dB/10)));
Pb_vit = qfunc(sqrt(5*10.^(Eb_N0_dB/10)));

close all;

semilogy(Eb_N0_dB,errors,'-+r');
hold on;
grid on;
semilogy(Eb_N0_dB,Pb_PAM,'b');
semilogy(Eb_N0_dB,Pb_vit,'k');
axis([Eb_N0_dB(1) Eb_N0_dB(end) 1e-7 1e-0]);
legend('Simulation','Uncoded','Approximation');
title('BER/SNR');

