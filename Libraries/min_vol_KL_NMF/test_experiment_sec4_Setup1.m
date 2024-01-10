clearvars;close all;clc

% Piano sample "Mary had a little lamb" with factorization rank=4

% Click on Run Button

% Parameter definition
[options]=parameters;
options.fileselection=2;
options.timecut=0;
options.RM=2;
options.WINDOWSIZE =512;  
options.hopsize=256;
options.winselection=1;
options.K=4; 
options.MAXITER=200;
options.algo={'disjointconstraint_minvol_KLNMF'};
options.sparsity=1;
options.beta=1;
options.obj_compute = 1;
% options.init=1; replaced by options.display
options.display = 1;
options.epsi=10^-6;
options.alpha=0;
options.lambda_tilde=0.3;
options.delta=1;
options.BM=0;
options.BMfun={'sparse','baseline'};
options.focus=0;
options.metric_eval=0; 
options.timemax = 5;

% Call of the main function
MyProject_NMF(options)