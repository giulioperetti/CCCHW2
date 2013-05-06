% Exercice 1 
% (MatLab implementation of Viterbi algorithm) We wish to implement in 
% MatLab a convolutional encoder and decoder (Viterbi algorithm) for 
% binary code
%
%   g_can = [1+D?2;1+D+D^2]
%
% and to evaluate its performance.
%
clear;

% binary code
g1 = [1 0 1];
g2 = [1 1 1];

% input
mu = 10^7; 
u = round(rand(1,mu));

% maps
[SS,OO,NN] = binconvmaps(g1,g2);

% encoding
y = cbinconvenc(u,SS,OO);
%y = binconvenc(u,SS,OO);

% modulation
s = pamap(y);

% channel simulation
Eb_N0_max = 8;
Eb_N0_dB = 0:Eb_N0_max;
snr_dB = Eb_N0_dB - 10*log10(0.5);
errors = zeros(1,Eb_N0_max);

for i=1:Eb_N0_max+1

    r = awgn(s,snr_dB(i),'measured');

    % viterbi decoding
    uhat = cbinconvdec(r,SS,OO,NN);
    %uhat = binconvdec(r,SS,OO,NN);

    errors(i) = sum(abs(u-uhat))/mu;
    
end

semilogy(Eb_N0_dB,errors);
grid on;
