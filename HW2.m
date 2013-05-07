% Exercice 1 
% (MatLab implementation of Viterbi algorithm) We wish to implement in 
% MatLab a convolutional encoder and decoder (Viterbi algorithm) for 
% binary code
%
%   g_can = [1+D^2;1+D+D^2]
%
% and to evaluate its performance.
%
clear;

% binary code
g1 = [1 0 1];
g2 = [1 1 1];

ord = 6;

% input
mu = 10^ord; 
u = round(rand(1,mu));

% maps
[SS,OO,NN] = binconvmaps(g1,g2);

% encoding
y = cbinconvenc(u,SS,OO);
%y = binconvenc(u,SS,OO);

% modulation
s = pamap(y);

% channel simulation
Eb_N0_min = 0;
Eb_N0_max = 11;
Eb_N0_dB = Eb_N0_min:Eb_N0_max;
snr_dB = Eb_N0_dB + 10*log10(2);
Ec_N0_dB = Eb_N0_dB - 10*log10(2);
errors = zeros(1,Eb_N0_max-Eb_N0_min);
%errors2 = zeros(1,Eb_N0_max-Eb_N0_min);

for i=1:Eb_N0_max-Eb_N0_min+1

    %r = awgn(s,snr_dB(i),'measured');
    n = 1/sqrt(2)*[randn(size(s)) + 1i*randn(size(s))]; % white gaussian noise, 0dB variance
    r = real(s + 10^(-Ec_N0_dB(i)/20)*n); % additive white gaussian noise
    
    % viterbi decoding
    uhat = cbinconvdec(r.',SS.',OO.',NN.');
 %   uhat2 = binconvdec(r,SS,OO,NN);

    errors(i) = sum(u~=uhat)/mu;
  %  errors2(i) = sum(u~=uhat2)/mu;
    
    fprintf('%d',i);
end
%%
Pb_PAM = qfunc(sqrt(2*10.^(Eb_N0_dB/10)));
Pb_vit = qfunc(sqrt(5*10.^(Eb_N0_dB/10)));

close all;

semilogy(Eb_N0_dB,errors,'-+r');
hold on;
grid on;
%semilogy(Eb_N0_dB,errors2,'-sg');
semilogy(Eb_N0_dB,Pb_PAM,'b');
semilogy(Eb_N0_dB,Pb_vit,'k');
axis([Eb_N0_dB(1) Eb_N0_dB(end) 1e-7 1e-0]);

