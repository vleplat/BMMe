%   This code solves the betaNMF problem
%
%          min_{W,H >= options.epsilon} D_{options.beta}(X,WH)
%
%   where options.beta is the beta divergence used to measure the error,
%   and options.epsilon is a small nonnegative constant (we recommend to
%   use the machine epsilon, which is the default value; cf. the discussion
%   about the zero locking phenomenon of MU in the book). 
%   See https://gitlab.com/ngillis/nmfbook/ and 
%       https://sites.google.com/site/nicolasgillis/book
%
% ****** Input ******
% X     :  the nonnegative input matrix pair
% r     :  the rank of the sought approximation
%
% ---Options---
% .extrapol   : choice of the extrapolation sequence (see the paper for more details) 
%               'nesterov' --> see the paper. Default for beta in [1,2]. 
%               'ptsengv1' --> (t-1)/t.
%               'ptsengv2' --> t/(t+1).
%               'noextrap' --> 0 (no extrapolation). Default for beta NOT in [1,2]. 
% .maxiter    : the maximum number of iterations performed by the algorithm
%             -default = 500.
% .timemax   : the maximum time in seconds alloted to the algorithm
%             -default = 60.
% .beta       : beta divergence considered
% .epsilon    : lower bound on the entries of W and H to ensure convergence
%             -default: Matlab machine precision 2^-52
% .accuracy   : stop iterations if the relative error does not decrease by
%               this value between two iterations (default: 1e-6)
% .W and .H    : initial values for W and H.
%                default: rand(m,r) and rand(r,n)
% .display in {0,1} : =1 displays the evolution of the iterations (default),
%                     =0 otherwise.
%
% ****** Output ******
% (W,H)     : W>=0 and H>=0, and WH approximates X according the the
%               criterion described above.
% e         : e gives the values of the objective functions during the
%             iterative process, that is, D_beta(X,WH)
%
% Code modified from https://gitlab.com/ngillis/nmfbook/

function [W,H,e,t] = betaNMF(X,r,options);

time0 = cputime;
[m,n] = size(X);
if nargin <= 2
    options = [];
end
if ~isfield(options,'maxiter')
    options.maxiter = 500;
end
if ~isfield(options,'timemax')
    options.timemax = 60;
end
if ~isfield(options,'beta')
    warning('Beta not specified, default selected: beta=1 => KL divergence.');
    options.beta = 1;
end
if ~isfield(options,'extrapol')
    if options.beta >= 1 && options.beta <= 2
        options.extrapol = 'nesterov';
    else
        options.extrapol = 'noextrap';
    end
end
if ~isfield(options,'epsilon')
    options.epsilon = eps;
end
if ~isfield(options,'accuracy')
    options.accuracy = 1e-4;
end
if ~isfield(options,'W')
    W = rand(m,r);
else
    W = max(options.epsilon,options.W);
end
if ~isfield(options,'H')
    H = rand(r,n);
else
    H = max(options.epsilon,options.H);
end
if ~isfield(options,'display')
    options.display = 1;
end
if options.beta == 2
    warning('Since beta=2, you might want to use more efficient algorithms; see FroNMF.m');
end
i = 1;
cpuinit = cputime;
% Scaling: the maximum entry in each column of W is 1
for k = 1 : r
    mxk = max( W(:,k) );
    W(:,k) = W(:,k)/mxk;
    H(k,:) = H(k,:)*mxk;
end
% Keep previous iterate in memory for extrapolation
Hp = H; Wp = W;
cparam = 10^30;
nutprev = 1;
if options.display == 1
    disp('Iterations:');
    cntdis = 0;
    mintime = 0.1; % display parameters
end
while i <= options.maxiter ...
        && cputime <= cpuinit+options.timemax...
        && (i <= 12 || abs(e(i-1)-e(i-11)) > options.accuracy*abs(e(i-2)))
    % Compute the extrapolated points
    stepH = max(0,H-Hp);
    if options.extrapol == 'nesterov'
        nut = 0.5 * ( 1 + sqrt(1+4*nutprev^2) );
        extrapolparam = (nutprev-1)/nut;
        nutprev = nut;
    elseif options.extrapol == 'ptsengv1'
        extrapolparam = (i-1)/i;
    elseif options.extrapol == 'ptsengv2'
        extrapolparam = i/(i+1);
    elseif options.extrapol == 'noextrap'
        extrapolparam = 0;
    end
    normstepH = norm(stepH,'fro');
    if normstepH == 0
        extrapH = extrapolparam;
    else
        extrapH = min( extrapolparam, cparam / i^(1.5/2) / normstepH);
    end
    He = H + extrapH * stepH;
    stepW = max(0,W-Wp);
    normstepW = norm(stepW,'fro');
    if normstepW == 0
        extrapW = extrapolparam;
    else
        extrapW = min(extrapolparam, cparam / i^(1.5/2) / normstepW);
    end
    We = W + extrapW * stepW;
    % Keep previous iterate in memory for extrapolation
    Hp = H; Wp = W;
    % Update of W and H with MU with extrapolation
    H = MUbeta(X,W,He,options.beta,options.epsilon);
    W = MUbeta(X',H',We',options.beta,options.epsilon);
    W = W';
    % Error: this time should not be taken into account in the
    % computational cost of the method
    if nargout >= 3
        timeei = cputime;
        e(i) = betadiv(X,W*H,options.beta);
        timeei = cputime-timeei; % do not take the computation of the error in the cost
        if i == 1
            t(i) = cputime-time0-timeei;
        else
            t(i) = t(i-1) + (cputime-timei)-timeei;
        end
    end
    % Display evolution of the iterations
    if options.display == 1
        if cputime-time0 >= mintime
            fprintf('%1.0d...',i);
            mintime = mintime*2;
            cntdis = cntdis+1;
            if mod(cntdis,10) == 0
                fprintf('\n');
            end
        end
    end
    i = i+1;
    timei = cputime;
end
if options.display == 1
    fprintf('\n');
end
