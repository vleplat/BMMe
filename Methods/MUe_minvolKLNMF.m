
%%-------------------------------------------------------------------------
% Multiplicative Updates with extrapolation for solving min-vol KL NMF, 
% see Sections 4.3 and 5.3 from the paper
%%-------------------------------------------------------------------------
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: 
%  - X      : an input matrix of size m by n
%  - r      : factorization rank
%  - struc options:  options include {display, init.W, init.H, maxiter, timemax, paraepsi,obj_compute}
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output: 
% - W, H    : the factors of the decomposition
% - e       : the objective sequence and 
% - t       : the sequence of cpu running time [sec]
%
% if options.obj_compute = 1 then the output (default value) then the output e would be 
% the sequence of the objective values;  otherwise, the output e would be [].
%

function [W,H,e,t] = MUe_minvolKLNMF(X,r,options) 

[m,n] = size(X); 
%% -------------------------------------------------------------------------
% Parameters of NMF algorithm
%%-------------------------------------------------------------------------
if nargin < 3
    options = [];
end
if ~isfield(options,'display')
    options.display = 1; 
end
if ~isfield(options,'init')
    W = 1 + rand(m,r); 
    H = 1 + rand(r,n); 
else
    W = options.init.W; 
    H = options.init.H; 
end
if ~isfield(options,'maxiter')
    options.maxiter = 200; 
end
if ~isfield(options,'no_ex')
    options.no_ex = 0; 
end
if ~isfield(options,'timemax')
    options.timemax = 5; 
end

if ~isfield(options,'paramepsi')
    options.paramepsi = eps; 
end

if ~isfield(options,'delta')
    options.delta = 1; 
end

if ~isfield(options,'accuracy_bisMethod')
    options.accuracy_bisMethod = 10^-8; 
end

if ~isfield(options,'lambda_tilde')
    options.lambda_tilde = 0.05;
end

if ~isfield(options,'projsplx')
    options.projsplx = 0;
end

if ~isfield(options,'obj_compute')
    options.obj_compute = 1; % if it is 0 then MU does not evaluate the objective, so e = []
end

% computation of min-vol penalty weight
options.lambda1 = options.lambda_tilde*(betaDiv(X+eps,W*H+eps,1))/(abs(log10(det(W'*W+options.delta*eye(r)))));
if options.display == 1
    fprintf(' -> Initial value for betadivergence = %0.2f \n', betaDiv(X+eps,W*H+eps,1));
    fprintf(' -> Value for the penalty weight = %0.2f \n', options.lambda1);
    fprintf(' -> Initial value for penalty term = %0.2f \n', options.lambda1 * log10(det(W'*W+options.delta*eye(r))));
    fprintf(' -> Initial ratio for logdet term = %d \n',(options.lambda1*(abs(log10(det(W'*W+options.delta*eye(r))))))/betaDiv(X+eps,W*H+eps,1));

end
%% -------------------------------------------------------------------------
% Initialization
%%-------------------------------------------------------------------------
W_prev=W;
H_prev=H;

mu=zeros(r,1);
e_m1 = ones(m,1);
e_n1 = ones(n,1);

cputime0 = tic; 
i = 1; 
t(1) = toc(cputime0);
e=[];
timeerr=0;
if options.obj_compute==1
    time1=tic;
    e(1)=  betaDiv(X+eps,W*H+eps,1) + options.lambda1*log10(det(W'*W+options.delta*eye(r)));
    timeerr=toc(time1); % to remove the time of finding the objective function
end
t1=1;
const=1e16;
% const=1e30;

%% -------------------------------------------------------------------------
% Main Optimization Loop
%%-------------------------------------------------------------------------

while i <= options.maxiter && t(i) < options.timemax  

    %%% update extrapolation coeff
    t2=1/2*(1+sqrt(1+4*t1*t1));
    ex_coef=(t1-1)/t2;
    t1=t2;
    i2=i^(1.5/2);

    %%% update H
    Wt=W';
    rj=sum(W)';
    rjc=repmat(rj,1,n);
    H_bar=max(H-H_prev,0);
    ex_coef1=min(ex_coef,const*1/i2/norm(H_bar,'fro'));
    if options.no_ex==1
        ex_coef1 = 0;
    end
    H_ex=H + ex_coef1*H_bar;
    bAx=X./(W*H_ex+eps);
    cj=Wt*bAx; 
    H_prev=H;
    H=max(options.paramepsi,(H_ex.*cj)./(rjc+eps));
     
    %%% update of W
    % extrapolation
    W_bar=max(W-W_prev,0);
    ex_coef1=min(ex_coef,const*1/i2/norm(W_bar,'fro'));
    if options.no_ex==1
        ex_coef1 = 0;
    end
    W_ex=W+ex_coef1*(W_bar);
    WWex_deltaI=W_ex'*W_ex+options.delta*eye(r);
    % Computation of L-smoothness constant
    L_phi1 = 2*norm(inv(WWex_deltaI),2);

    % Computation of A matrix
    A = 2*W_ex/WWex_deltaI;

    % Computation of V_tilde and B1
    B1 = (((X)./(W_ex*H + options.paramepsi))*H').*W_ex;

    % Computation of \mu_k
    for k=1:r
        % Computation of \tilde{\mu}_k
        mu_tilde = 1/options.lambda1*(4*options.lambda1*L_phi1*B1(:,k)*m - 1/m - sum(H(k,:))*e_m1) + L_phi1*W_ex(:,k) - A(:,k);
        % Call of bisection method
        mu(k) = bisectionMethod(min(mu_tilde),max(mu_tilde),options.accuracy_bisMethod,B1,A,H,options.lambda1,L_phi1,W_ex,k);
    end

    % Computation of B2 and actual update of W
    B2 = e_m1*(e_n1'*H') + options.lambda1*(A - L_phi1*W_ex + e_m1*mu');
    W_prev = W;
    W = max(options.paramepsi,1/2*(-B2 + (B2.^2 + 4*options.lambda1*L_phi1*B1).^(1/2)));
   
    %%% Update iteration counter
    i=i+1; 

    %%% Computation of error functions and time evolution
    if options.obj_compute==1
        time1=tic; 
        e(i)= betaDiv(X+eps,W*H+eps,1) + options.lambda1*log10(det(W'*W+options.delta*eye(r)));
        timeerr = timeerr + toc(time1);
    end
    t(i)= toc(cputime0) - timeerr;

end
if(options.display==1)
    fprintf(' -> Final value for betadivergence : %0.2f \n', betaDiv(X+eps,W*H+eps,1));
    fprintf(' -> Final value for penalty term : %0.2f \n', options.lambda1 * log10(det(W'*W+options.delta*eye(r))));
    fprintf(' -> Ratio of terms : %0.2f \n', (options.lambda1*(abs(log10(det(W'*W+options.delta*eye(r))))))/betaDiv(X+eps,W*H+eps,1));
end
end


