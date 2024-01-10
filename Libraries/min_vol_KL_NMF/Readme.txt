%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.INTRODUCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This framework is dedicated to the blind audio source separation and implements the algorithms 
presented in the following references:
[1]: Valentin Leplat, Nicolas Gillis and Man Shun Ang, "Blind Audio Source Separation with Minimum-Volume Beta-Divergence NMF", 2019, IEEE TSP
[2]: Valentin Leplat, Nicolas Gillis and Jerôme Idier, "Multiplicative Updates for NMF with beta-Divergences under Disjoint Equality Constraints", preprint, 2020
In order to reproduce the results presented in [2] (section 4): launch the following MatLab files:
- test_experiment_sec4_Setup1: results obtained for the piano sample "Mary had a little lamb" with factorization rank=3
- test_experiment_sec4_Setup2: results obtained for the piano sample "Mary had a little lamb" with factorization rank=7
- test_experiment_sec4_Setup3: results obtained for the piano sample "Prelude" from J.S. Bach with factorization rank=16
- test_experiment_sec4_Setup4: additional results obtained for the piano sample from Févotte, Bertin and Durieux with factorization rank=7
- cputime_contest_disjointconstraint_minvol_KLNMF.m: results for Table 3 from [2], a second launch must be made to get the results with Algorithm min-vol KL-NMF from [1] by setting options.algo={'minvol'}; in each .m called by the batchtest
       please read carefully the instructions at the beginning of the file
- lossfun_contest_disjointconstraint_minvol_KLNMF.m: results for Table 4 from [2]        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
2.Global Structure of the framework
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The framework contains the following (most important) files:
- parameters.m : used to define all the parameters for the audio source separation, refer to section 3 or 
 open the file and see the comments for the description of the different parameters  that have to be defined to launch an analysis.
- MyProject_NMF.m: main file to be launched (no modification to do inside this file)
- betaNMF_baseline.m: algorithm to solve beta-NMF wihtout penalty term 
- betaNMF_logdetKL.m: algorithm min-vol KL-NMF described in [1]
- disjointconstraint_minvol_KLNMF= New algorithm to tackle min-vol KL-NMF [2]


In the directory "Utils", several utilitary matlab files are present such as the computation of the STFT, the STFT inverse, 
audio sample selection, etc.. 

You have the possibility to launch two kinds of analysis:
- A standalone analysis, i.e., one audio source separation performed with one algorithm. See section 3 for the description of the parameters
- A benchmark analysis: a compartive study between three different algorithms: min-vol NMF/disjointconstraint_minvol_KLNMF, baseline NMF and sparse NMF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
3.Launch a standalone analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
- Open file "parameters.m" and set the parameters with the values you want. Here-under you can find a description of the parameters:
  - "options.fileselection": you can whether choose a sample considered in [1] (set the parameter to 1,2 or 3)
   or you can input your own audio file in wav format => rename the file as "myfile", put it in subdirectory "Mysample"
   and set the parameter "options.fileselection" to 99
  - "options.timecut": you can cut the full audio signal and select one chunk (defined by the time interval in the 
   parameter "options.timecut_interval") to be processed.
  - "options.RM": two techniques to reconstruct the source estimates spectrograms S_i from the NMF results (W and H):
      -> Synthesis: S_i = W(:,i)*H(i,:)
      -> Wiener Filter / Masking coefficients: S_i= (W(:,i)*H(i,:)).(W*H) * V
  - "options.WINDOWSIZE": size of the window to compute the STFT, usual values: 256,512 and 1024
  - "options.hopsize": hop size for two consecutive windows (number of points between two windows): in order to get
   50% overlapping (usual value) between successive windows, just set this to half value of options.WINDOWSIZE
  - "options.winselection": chose the window type (Hamming, Hann, Blackman, Sinebell)
  - "options.K": the factorization rank for the NMF
  - "options.MAXITER": the maximum number of iterations for the optimization loop used to "solve" NMF problems
  - "options.algo": the algorithm, you can choose between min-vol [1], disjointconstraint_minvol_KLNMF [2], baseline or sparse. 
  - "options.sparsity": sparsity penalty weight for sparse-NMF
  - "options.beta": value for beta. Here beta = 1, the cost function considered will correspond to Kullback-Leibler (KL) divergence
                    
  - "options.init": activate (set to 0) of deactivate (set to 1) some real-time plots that give the evolution of the objective functions and gamma parameter
  - "options.gamma": initial value parameter for line-search procedure, see [1]. We advice to use 0.5
  - "options.alpha": a parameter used to slightly accelerate the convergence of baseline NMF and enforce a little bit
  of sparsity, set the value in the interval 0.001 - 0.01
  - "options.lambda_tilde": initial ratio of the logdet term over the total objective function: choose a value in the
  intervale 0.001 and 0.2, depending on the data: several tries are in order to get a good estimation for this parameter
  
  - "options.BM": MUST BE SET TO ZERO for a standalone analysis
  - "options.focus": a post-processing parameter for the figures to zoom on a specific time interval (defined with 
  the parameter "options.time_interval")
  - "options.metric_eval": If you have the true sources and you want to compute the quality of the source separation using metrics SDR, SIR and SAR [3]
    , put the audio files of the true sources in the subdirectory './Metric_eval/true_sources' and set the parameter to 1
    the metrics are displayed on the Command Window of MatLab

=> an efficient way to build up your analysis, create a .m file dedicated and follow the structure of an example such as "test_experiment_sec4_Setup1.m"

The framework generates several files:
- All kinds of figures helpfull for an audio source separation (Input signals, input amplitude spectrogram, results for W
  and H, mask coefficients, etc) are generated to help the user to analyze the results of the source separation.
  The figures are saved in subdirectory "Graphs" in .fig and .png formats
- The audio files for the source estimates in wav format. The audio files are put in subdirectory "Extractions"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
4.Launch a benchmark analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Most of the parameters to define are described in the previous section. BUT, the following parameters require a 
particular attention:
- "options.BM" must be set to 1 
- "options.algo={'minvol'}"[1] or "options.algo={'disjointconstraint_minvol_KLNMF'}"[2] must be considered
- "options.BMfun" must be defined with available algorithms to be compared with min-vol NMF: sparse and/or baseline
 For a full comparative study, just set the parameter as follows: options.BMfun={'sparse','baseline'};
-"options.metric_eval" must be set to zero (not yet available for benchmark tests)


The figures generated are stored in subdiretory "BenchMark", the audio files for the source estimates obtained with each
algorithm are stored in subdirectories: "BenchMark/Extractions/minvol", "BenchMark/Extractions/sparse" and/or "BenchMark/Extractions/baseline"


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Additional References
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[3]: BSS eval tool box: http://bass-db.gforge.inria.fr/bss_eval/
