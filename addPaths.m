clear all;close all;clc;

%%-------------------------------------------------------------------------
% Add paths to directories
%--------------------------------------------------------------------------
addpath("Datasets\");   % directory with data sets
addpath("Methods\");    % directory with Proposed Algorithms
addpath("Utils\");      % directory with utiliterian functions

% below: libraries for benchmarked methods
addpath("Libraries\group_robust_NMF\") ;
addpath("Libraries\min_vol_KL_NMF\");
addpath("Libraries\min_vol_KL_NMF\Utils_lib\");
