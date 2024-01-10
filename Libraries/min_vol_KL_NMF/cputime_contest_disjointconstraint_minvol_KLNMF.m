close all;clc;clear all;
% cputime comparaison for algorithms
% ATTENTION: Before launching this file:
%                 - comment line 1 in files "test_experiment_sec4_Setup1.m", "test_experiment_sec4_Setup2.m" and "test_experiment_sec4_Setup3.m"
%                 - set options.BM=1; and options.init=1;
%                 - a second launch must be made to get the results with Algorithm min-vol
%                     KL-NMF from [1] by setting options.algo={'minvol'}; in each file .m called and replace line 96 by "rowHeaders{2} = sprintf('min-vol KL NMF');" 
%                     in this file
J_ini=20;
%% Experiment 1
tablecpu_exp1=zeros(J_ini,3);
for j=1:J_ini
    test_experiment_sec4_Setup1
    load variables.mat;
    tablecpu_exp1(j,1)=cputime_logdetKL;
    tablecpu_exp1(j,2)=cputime_s;
    tablecpu_exp1(j,3)=cputime_b;
    close all;
    clc;
end
% Mean and std computation 
tablecpu_exp1_MU=mean(tablecpu_exp1,1);
tablecpu_exp1_std=std(tablecpu_exp1,0,1);


%% Experiment 2
tablecpu_exp2=zeros(J_ini,3);
for j=1:J_ini
    test_experiment_sec4_Setup2
    load variables.mat;
    tablecpu_exp2(j,1)=cputime_logdetKL;
    tablecpu_exp2(j,2)=cputime_s;
    tablecpu_exp2(j,3)=cputime_b;
    close all;
    clc;
end
% Mean and std computation 
tablecpu_exp2_MU=mean(tablecpu_exp2,1);
tablecpu_exp2_std=std(tablecpu_exp2,0,1);


%% Experiment 3
tablecpu_exp3=zeros(J_ini,3);
for j=1:J_ini
    test_experiment_sec4_Setup3
    load variables.mat;
    tablecpu_exp3(j,1)=cputime_logdetKL;
    tablecpu_exp3(j,2)=cputime_s;
    tablecpu_exp3(j,3)=cputime_b;
    close all;
    clc;
end
% Mean and std computation 
tablecpu_exp3_MU=mean(tablecpu_exp3,1);
tablecpu_exp3_std=std(tablecpu_exp3,0,1);
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
save variables_cpuContest.mat


% Display table
columnHeaders = {'Experiment #1','Experiment #2', 'Experiment #3'};
rowHeaders{1} = sprintf('baseline KL-NMF');
tableData{1,1} = sprintf('%.2f %s %.2f', tablecpu_exp1_MU(3), 177, tablecpu_exp1_std(3));
tableData{1,2} = sprintf('%.2f %s %.2f', tablecpu_exp2_MU(3), 177, tablecpu_exp2_std(3));
tableData{1,3} = sprintf('%.2f %s %.2f', tablecpu_exp3_MU(3), 177, tablecpu_exp3_std(3));

rowHeaders{2} = sprintf('Algorithm 4');
tableData{2,1} = sprintf('%.2f %s %.2f', tablecpu_exp1_MU(1), 177, tablecpu_exp1_std(1));
tableData{2,2} = sprintf('%.2f %s %.2f', tablecpu_exp2_MU(1), 177, tablecpu_exp2_std(1));
tableData{2,3} = sprintf('%.2f %s %.2f', tablecpu_exp3_MU(1), 177, tablecpu_exp3_std(1));

rowHeaders{3} = sprintf('sparse KL-NMF');
tableData{3,1} = sprintf('%.2f %s %.2f', tablecpu_exp1_MU(2), 177, tablecpu_exp1_std(2));
tableData{3,2} = sprintf('%.2f %s %.2f', tablecpu_exp2_MU(2), 177, tablecpu_exp2_std(2));
tableData{3,3} = sprintf('%.2f %s %.2f', tablecpu_exp3_MU(2), 177, tablecpu_exp3_std(2));

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
set(hTable, 'ColumnWidth', {100, 100, 100});
set(gcf,'name','Average Performance of the algorithms','numbertitle','off')
