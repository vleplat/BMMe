clearvars;close all;clc

% "Prelude" from J.S. Bach with factorization rank=16 

% Click on Run Button

% Parameter definition
[options]=parameters;
options.fileselection=3;
options.timecut=0;
options.RM=2;
options.WINDOWSIZE =1024;  
options.hopsize=512;
options.winselection=1;
options.K=16; 
options.MAXITER=300;
options.obj_compute = 1;
options.algo={'disjointconstraint_minvol_KLNMF'};
options.sparsity=1;
options.beta=1;
options.display = 1;
options.epsi=10^-6;
options.alpha=0;
options.lambda_tilde=0.04;
options.delta=1;
options.BM=0;
options.BMfun={'sparse','baseline'};
options.focus=0;
options.metric_eval=0; 

% Call of the main function
MyProject_NMF(options)