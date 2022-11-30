clc; clear;
close all;
rng('default')
addpath(genpath('utils'));
var_struct = load('container.mat');
name_cell = fieldnames(var_struct);
Input = double(getfield(var_struct,char(name_cell)));
Input = Normalize(Input);
dim = size(Input);
%% Sampling
SamRate = 0.1; % this can be tuned
A = PermuteWHT_partitioned(dim(1)*dim(2),dim(3),SamRate);
y = A*Input(:);
m = length(y);
NoiseLevel = 0;
noise = NoiseLevel * randn(m,1);
y = y + noise;
%% Solving by TCTV-DFT
fprintf('=========== TCTV-DFT ============\n');
opts.MaxIter = 200;
opts.tol = 1e-8;
opts.shift_dim = [1,3,2];
opts.X = Input;
opts.dim = dim;
opts.mu = 1e-3;
opts.rho = 1.1;
transform.L = @fft; transform.l = dim(2); transform.inverseL = @ifft;
% transform.L = dctmtx(dim(2)); transform.l = 1;
% transform.L = RandOrthMat(dim(2)); transform.l = 1;
opts.transform = transform;
lambda = 1e-1; % this can be tuned
xrec_tctv = tctv_cs(A, y, lambda, opts);




