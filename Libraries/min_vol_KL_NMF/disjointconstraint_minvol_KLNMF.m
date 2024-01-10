function [W, H, lossfun, t] = disjointconstraint_minvol_KLNMF(V,options)

% This function tackles the NMF problem with volume constraints and normalization on the 
% factor W using the KL divergence: 
% 
% min_{W,H >= 0} D_KL(V | WH) + lambda * logdet(W^TW + delta I), 
% 
% where the paramters (lambda, delta) are specified via options
% and the columns of W must lie on the unit simplex: W(:,k) \in \Delta^F
% for all k. 
% Note that lambda is specified via lambda_tilde: lambda will
% be equal to lambda_tilde * D_KL(V | WH) / logdet(W^TW + delta I) in
% order to balance the two terms properly in the objective, where (W,H) is 
% the initialization.  

%INPUTS:
%       V:                 Spectrogram of x (input audio signal) 
%       options.K:         Number of sources (determined prior to the
%       factorization)
%       options.MAXITER:   Maximum number of iterations for updates of W and H
%       options.lambda_tilde:   Relative Weight for Volume minimization
%       options.delta:     Parameter within the logdet function
%       options.epsi:      Stopping criterion for Newton-Raphson
%                          method


%OUTPUTS:
%       W: Matrix dictionnaries
%       H: Matrix activations
%       lossfun: value of the loss function f(W,H) after MAXITER iterations
%       t: time in seconds to return a solution

% % Loading parameters

K=options.K;
beta=1;
delta=options.delta;
alpha=options.alpha;
epsi=options.epsi;

% % Algorithm Beta-NMF with logdet(W'W+delta*I) regularization
% disp(' ->disjoint-constraint min-vol KL-NMF')
tic
F = size(V,1);
N = size(V,2);
addpath('Libraries\min_vol_KL_NMF\Utils_lib\');

% %  Initialization for W and H
if ~isfield(options,'init')
    % disp(' ->Random Initialization for W and H')
    % Random initialization
    W = 1+rand(F, K);
    H = 1+rand(K, N);
else
    % disp(' ->Initialization given for W and H')
    W = options.init.W; 
    H = options.init.H; 
end

if ~isfield(options,'MAXITER')
    options.MAXITER = 200; 
end
if ~isfield(options,'timemax')
    options.timemax = 5; 
end

% % Initialization for loop parameters
Y = inv(W'*W + delta*eye(K));

% % array to save the value of the loss function
% lossfunsave = zeros(options.MAXITER,1);
lossfunsave = [];

% % initialization for lambda 
lambda=options.lambda_tilde*betaDiv(V+eps,W*H+eps,beta)/abs(log10(det(W'*W+delta*eye(K))));
if(options.display==1)
    fprintf(' -> Initial value for betadivergence : %0.2f \n', betaDiv(V+eps,W*H+eps,beta));
    fprintf(' -> Value for the penalty weight : %0.2f \n', lambda);
    fprintf(' -> Initial value for penalty term : %0.2f \n', lambda * log10(det(W'*W+delta*eye(K))));
    fprintf(' -> Initial ratio of terms : %0.2f \n', lambda * log10(det(W'*W+delta*eye(K)))/betaDiv(V+eps,W*H+eps,beta));
end

% % initialization for mu (Lagrangian multipliers) 
mu = zeros(K,1);
% Others parameters/variables
ONES = ones(F,N);
JF1  = ones(F,1);
if(options.display==1)
    ani1=animatedline('Color','r','Marker','o');
    title(['Evolution of objective function - $\beta$ = ' num2str(beta)],'FontSize',12, 'Interpreter','latex')
    xlabel('iteration','FontSize',12, 'Interpreter','latex')
end

% % Optimization loop
% iter = 1;
% t(1) = toc(cputime0);
% time1=tic;
% WtW = W'*W;
% lossfunsave(iter) = betaDiv(V+eps,W*H+eps,beta) + lambda * log10(det(WtW+delta*eye(K)));
% timeerr=toc(time1);

cputime0 = tic;
iter = 1; 
t(1) = toc(cputime0);

timeerr=0;
if options.obj_compute==1
    time1=tic;
    lossfunsave(iter) = betaDiv(V+eps,W*H+eps,beta) + lambda * log10(det(W'*W+delta*eye(K)));
    timeerr=toc(time1); % to remove the time of finding the objective function
end


while iter<=options.MAXITER && t(iter) < options.timemax
    
    % % update matrix  H Coefficients ("activations")
    H = (H .* (W'*(((W*H).^(beta-2)).*V))./(W'*(W*H).^(beta-1)+eps)).^(1+alpha);
    H = max(H,eps);
    
    %  % update mu (lagrangian multipliers)
    Y_plus=max(0,Y);
    Y_minus=max(0,-Y);
    C=ONES*H'-4*lambda*W*Y_minus;
    S=8*lambda*(W*(Y_plus+Y_minus)).*((V./(W*H+eps))*H');
    D=4*lambda*W*(Y_plus+Y_minus);
    mu=updatemu(C,S,D,W,mu,epsi);
    %  % update maxtrix W ("dictionaries")
    W = W .*(((C+JF1*mu').^2+S).^(1/2)-(C+JF1*mu'))./(D+eps);
    W = max(W,eps);

    % % Update iteration counter
    iter = iter + 1; 
    
    % % Compute the loss function = Beta-Divergence + Penalty term 
    time1 = tic; 
    WtW = W'*W;
    lossfunsave(iter) = betaDiv(V+eps,W*H+eps,beta) + lambda * log10(det(WtW+delta*eye(K)));
    timeerr=timeerr + toc(time1);
    t(iter)= toc(cputime0)-timeerr;
    
    % % Update Y
    Y = inv(W'*W + delta*eye(K));
    
    % Drawing
    if(options.display==1)
        if(options.delta<1)
            addpoints(ani1,iter,lossfunsave(iter));
        else
            addpoints(ani1,iter,log10(lossfunsave(iter))) ;
        end
        drawnow
    end

    
end


% % Save lossfun
lossfun=lossfunsave;
% % Final Value for loss function
if(options.display==1)
    fprintf(' -> Final value for betadivergence : %0.2f \n', betaDiv(V+eps,W*H+eps,beta));
    fprintf(' -> Final value for penalty term : %0.2f \n', lambda * log10(det(W'*W+delta*eye(K))));
    fprintf(' -> Ratio of terms : %0.2f \n', lambda * log10(det(W'*W+delta*eye(K)))/betaDiv(V+eps,W*H+eps,beta));
end
end %EOF