% Lab2_1
clear all; close all; clc;
addpath(genpath('.\myFunctions'));

% Initial Parameters
% initial parity matrix
H = [1 1 0 1 1 0 0 1 0 0;           
     0 1 1 0 1 1 1 0 0 0;
     0 0 0 1 0 0 0 1 1 1;
     1 1 0 0 0 1 1 0 1 0;
     0 0 1 0 0 1 0 1 0 1];
 
% bits to encode     
bits = repmat([1 0 0 1 1],10,1);

bits_v = bits';
bits_v = bits_v(:)'; 
% maximum number of iterations    
maxIterations = 10;    

% Encoder
[encodedBits,H] = SSLdpcEncoder(bits,H);

% Errors due to noisy the channel
error =  eye(10);            % introducing one error in a different position.
noisyBits = mod (encodedBits + error,2);
    noisyBits = noisyBits';
    noisyBits = noisyBits(:)'; 

% Decoder
decodedBits = LDPCDecode(noisyBits,H,maxIterations);

N_errorBit = length(find(decodedBits - bits_v ~= 0));
BER =  N_errorBit/ length(bits_v)
    