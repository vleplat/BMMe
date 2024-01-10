clearvars;close all;clc

% Piano sample "Mary had a little lamb" with factorization rank=3

% Click on Run Button

% Parameter definition
[options]=parameters;
options.fileselection=2;
options.timecut=0;
options.RM=2;
options.WINDOWSIZE =512;  
options.hopsize=256;
options.winselection=1;
options.K=3; 
options.MAXITER=200;
options.sparsity=1;
options.beta=1;
options.init=1;
options.epsi=10^-6;
options.alpha=0;
options.lambda_tilde=0.1;
options.delta=1;
options.gamma=0.5;
options.BM=0;
options.BMfun={'sparse','baseline'};
options.focus=0;
options.metric_eval=0; 



J_ini=1;
%% Experiment 1
tablelossfun_exp1=zeros(J_ini,4);
for j=1:J_ini
    options.algo={'minvol'};
    MyProject_NMF(options)
    load variables.mat;
    tablelossfun_exp1(j,1)=betaDiv(V+eps,W*H+eps,options.beta);
    tablelossfun_exp1(j,2)=lossfunction(end);
    close all;
    clc;
    options.algo={'disjointconstraint_minvol_KLNMF'};
    MyProject_NMF(options)
    load variables.mat;
    tablelossfun_exp1(j,3)=betaDiv(V+eps,W*H+eps,options.beta);
    tablelossfun_exp1(j,4)=lossfunction(end);
    close all;
    clc;
end
% Mean and std computation 
tablelossfun_exp1_MU=mean(tablelossfun_exp1,1);
tablelossfun_exp1_std=std(tablelossfun_exp1,0,1);


%% Experiment 2
% Parameter definition
[options]=parameters;
options.fileselection=2;
options.timecut=0;
options.RM=2;
options.WINDOWSIZE =512;  
options.hopsize=256;
options.winselection=1;
options.K=7; 
options.MAXITER=200;
options.sparsity=1;
options.beta=1;
options.init=1;
options.epsi=10^-6;
options.alpha=0;
options.lambda_tilde=0.1;
options.delta=1;
options.BM=0;
options.BMfun={'sparse','baseline'};
options.focus=0;
options.metric_eval=0; 
tablelossfun_exp2=zeros(J_ini,4);
for j=1:J_ini
    options.algo={'minvol'};
    MyProject_NMF(options)
    load variables.mat;
    tablelossfun_exp2(j,1)=betaDiv(V+eps,W*H+eps,options.beta);
    tablelossfun_exp2(j,2)=lossfunction(end);
    close all;
    clc;
    options.algo={'disjointconstraint_minvol_KLNMF'};
    MyProject_NMF(options)
    load variables.mat;
    tablelossfun_exp2(j,3)=betaDiv(V+eps,W*H+eps,options.beta);
    tablelossfun_exp2(j,4)=lossfunction(end);
    close all;
    clc;
end
% Mean and std computation 
tablelossfun_exp2_MU=mean(tablelossfun_exp2,1);
tablelossfun_exp2_std=std(tablelossfun_exp2,0,1);


%% Experiment 3
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
options.sparsity=1;
options.beta=1;
options.init=1;
options.epsi=10^-6;
options.alpha=0;
options.lambda_tilde=0.022;
options.delta=1;
options.BM=0;
options.BMfun={'sparse','baseline'};
options.focus=0;
options.metric_eval=0; 
tablelossfun_exp3=zeros(J_ini,4);
for j=1:J_ini
    options.algo={'minvol'};
    MyProject_NMF(options)
    load variables.mat;
    tablelossfun_exp3(j,1)=betaDiv(V+eps,W*H+eps,options.beta);
    tablelossfun_exp3(j,2)=lossfunction(end);
    close all;
    clc;
    options.algo={'disjointconstraint_minvol_KLNMF'};
    MyProject_NMF(options)
    load variables.mat;
    tablelossfun_exp3(j,3)=betaDiv(V+eps,W*H+eps,options.beta);
    tablelossfun_exp3(j,4)=lossfunction(end);
    close all;
    clc;
end
% Mean and std computation 
tablelossfun_exp3_MU=mean(tablelossfun_exp3,1);
tablelossfun_exp3_std=std(tablelossfun_exp3,0,1);
%% Table : Average performance of the algorithms 

clear cond1
clear cond2
clear cond3
clear cputime_b
clear cputime_logdetKL
clear cputime_s
clear extensionFile
clear freq
clear H
clear H_b
clear H_s
clear lossfunction
clear lossfunctionb
clear myFolder
clear objective
clear options
clear params
clear T
clear time
clear time_frame_max
clear time_frame_min
clear tmax
clear tmin
clear V
clear W
clear W_b
clear W_s
clear x
clear X
save variables_lossfunContest.mat


% Display table
columnHeaders = {'Experiment #1-betaDiv','Experiment #1-F', 'Experiment #2-betaDiv','Experiment #2-F','Experiment #3-betaDiv','Experiment #3-F'};
rowHeaders{1} = sprintf('minvol KL-NMF');
tableData{1,1} = sprintf('%.2f %s %.2f', tablelossfun_exp1_MU(1), 177, tablelossfun_exp1_std(1));
tableData{1,2} = sprintf('%.2f %s %.2f', tablelossfun_exp1_MU(2), 177, tablelossfun_exp1_std(2));
tableData{1,3} = sprintf('%.2f %s %.2f', tablelossfun_exp2_MU(1), 177, tablelossfun_exp2_std(1));
tableData{1,4} = sprintf('%.2f %s %.2f', tablelossfun_exp2_MU(2), 177, tablelossfun_exp2_std(2));
tableData{1,5} = sprintf('%.2f %s %.2f', tablelossfun_exp3_MU(1), 177, tablelossfun_exp3_std(1));
tableData{1,6} = sprintf('%.2f %s %.2f', tablelossfun_exp3_MU(2), 177, tablelossfun_exp3_std(2));

rowHeaders{2} = sprintf('Algorithm 4');
tableData{2,1} = sprintf('%.2f %s %.2f', tablelossfun_exp1_MU(3), 177, tablelossfun_exp1_std(3));
tableData{2,2} = sprintf('%.2f %s %.2f', tablelossfun_exp1_MU(4), 177, tablelossfun_exp1_std(4));
tableData{2,3} = sprintf('%.2f %s %.2f', tablelossfun_exp2_MU(3), 177, tablelossfun_exp2_std(3));
tableData{2,4} = sprintf('%.2f %s %.2f', tablelossfun_exp2_MU(4), 177, tablelossfun_exp2_std(4));
tableData{2,5} = sprintf('%.2f %s %.2f', tablelossfun_exp3_MU(3), 177, tablelossfun_exp3_std(3));
tableData{2,6} = sprintf('%.2f %s %.2f', tablelossfun_exp3_MU(4), 177, tablelossfun_exp3_std(4));


% Create the table and display it.
hTable = uitable();
% Apply the row and column headers.
set(hTable, 'RowName', rowHeaders);
set(hTable, 'ColumnName', columnHeaders);
% Display the table of values.
set(hTable, 'data', tableData);
% Size the table.
set(hTable, 'units', 'normalized');
set(hTable, 'Position', [.1 .1 .8 .8]);
set(hTable, 'ColumnWidth', {100, 100, 100, 100, 100, 100});
set(gcf,'name','Average Performance of the algorithms','numbertitle','off')
