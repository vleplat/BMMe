function [options]=parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Audio file selection
options.fileselection=4;
% Samples provided
%1: Synthetic Bass and Drum example 1       
%2: Piano sample "Mary had a little lamb"
%3: Prelude Jean-Sebastian Bach 
%4: Piano sample from Fevotte, Bertin and Durieux
% Your sample: set the parameter to 99, rename your wav file as "myfile"
% and put it in subdirecory "Mysample"

% parameter to select a time interval ("cut") from the full 
% input signal to be analyzed:
options.timecut=0;
%0: intervals for time cut not active
%1: intervals for time cut active
options.timecut_interval=[0 2];
%define the time interval [a b] with a<b expressed in seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Reconstruction method
options.RM=2;
%1: Reconstruction of sources by ISTFT by "Synthesis"
%2: Reconstruction of sources by ISTFT by "Wiener Filter / Masking coefficients "

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % STFT parameters
% windows length
options.WINDOWSIZE =1024;  
% hopsize 
options.hopsize=512;
% choose the window type
options.winselection=1;
%1: Hamming window       
%2: Hann window
%3: Blackman window             
%4: Sinebell window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % NMF parameters
% Number of source estimates
options.K=7; 
% Max Number of iterations for the factorization
options.MAXITER=500;
% Algorithm selection
options.algo={'minvol'};
% 'minvol','baseline','sparse'
% if not 'minvol', benchmarking is deactivated: minvol is the reference for
% the benchmarking
% if sparse is selected, specify the sparsity penalty weight for sparsity
% penalty term here-under:
options.sparsity=1;

% beta value
options.beta=1; %keep beta=1 
%1: KullBack Leibler Divergence 

% options.init=0;
%0: graphics (lossfunction evolution) in real time activated
%1: graphics (lossfunction evolution) in real time deactivated

options.gamma=0.5;
% 0<gamma<=1: initial value for relaxation weight for logdetMU in order to ensure decrease
% of the loss function 
% -1: Cancel gamma effect

options.alpha=0;
% paramter to accelerate  convergence of entries of H (activations)
% choose value between 0.01 and 0.1 (only available for beta-NMF baseline)

% penalty weight for min-volume 
options.lambda_tilde=0.03; % initial ratio of the min-vol penalty, choose 
%it between 0.001 and 0.2, depending on the data
options.delta=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % BenchMark options
options.BM=0;
%0: No BenchMarking test - BLOC 4 is not active
%1: BenchMarking test on - BLOC 3 not active, BLOC 4 active
options.BMfun={'sparse','baseline'};
% Specify the functions to consider for the benchmark vs minvol
% 'baseline' and/or 'sparse'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Post-Processing parameters

% parameter to define a time interval for the post-processing only, i.e.,
% the figures present results on this time interval only, but the full
% audio signal is considered for the analysis:
options.focus=0;
%0: intervals for post-processing not active
%1: intervals for post-processing active
options.time_interval=[2 3];
%define the time interval [a b] with a<b expressed in seconds

options.metric_eval=0; 
%0: deactivate bss metric evaluation for stand alone run (not benchmark)
%1: activate bss metric evaluation for stand alone run (not benchmark)

end%EOF
