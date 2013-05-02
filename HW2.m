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
mu = 20; 
u = round(rand(1,mu));
%u = [1 0 1];

% maps
[SS,OO,NN] = binconvmaps(g1,g2);

% encoding
y = cbinconvenc(u,SS,OO);
%y = binconvenc(u,SS,OO);

% modulation
r = pamap(y);

% viterbi decoding
uhat = binconvdec(r,SS,OO,NN);

