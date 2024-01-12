%% Inertial Multiplicative Update for solving KL NMF
% Input: X, r, options
% options include {display, init.W, init.H, maxiter, timemax, paraepsi,obj_compute}
% Output: W, H (the factors), e (the objective sequence) and t (the sequence of
% running time)
%
% if options.obj_compute = 1 then the output (default value) then the output e would be 
% the sequence of the objective values;  otherwise, the output e would be [].
%
% this version uses t_Nesterov in the extrapolation sequences. 
% written by LTK Hien
% Latest update: December 2023
function [W,H,e,t] = MUe_KLNMF(X,r,options) 
cputime0 = tic; 
[m,n] = size(X); 
%% Parameters of NMF algorithm
if nargin < 3
    options = [];
end
if ~isfield(options,'display')
    options.display = 1; 
end
if ~isfield(options,'init')
    W = rand(m,r); 
    H = rand(r,n); 
else
    W = options.init.W; 
    H = options.init.H; 
end

if ~isfield(options,'const')
    options.const = 1e30; 
end

if ~isfield(options,'maxiter')
    options.maxiter = 200; 
end
if ~isfield(options,'timemax')
    options.timemax = 5; 
end

if ~isfield(options,'paramepsi')
    options.paramepsi = eps; 
end

if ~isfield(options,'obj_compute')
    options.obj_compute = 1; % if it is 0 then MU does not evaluate the objective, so e = []
end

const=options.const;
W_prev=W;
H_prev=H;

i = 1; 
t(1) = toc(cputime0);
e=[];
timeerr=0;
if options.obj_compute==1
    time1=tic;
    e(1)= KLobj(X,W,H);
    timeerr=toc(time1); % to remove the time of finding the objective function
end
t1=1;
kdis = 0; 
while i <= options.maxiter && t(i) < options.timemax  
    t2=1/2*(1+sqrt(1+4*t1*t1));
    ex_coef=(t1-1)/t2;
    t1=t2;

    i2=i^(1.5/2);
    %update H
    Wt=W';
    rj=sum(W)';
    rjc=repmat(rj,1,n);
    H_bar=max(H-H_prev,0);
    ex_coef1=min(ex_coef,const*1/i2/norm(H_bar,'fro'));
    H_ex=H + ex_coef1*H_bar;
    bAx=X./(W*H_ex+eps);
    cj=Wt*bAx; 
    H_prev=H;
    H=max(options.paramepsi,(H_ex.*cj)./(rjc+eps));
     
    %update W
   
    rj=sum(H,2)';
    rjc=repmat(rj,m,1); 
    W_bar=max(W-W_prev,0);
    ex_coef1=min(ex_coef,const*1/i2/norm(W_bar,'fro'));
    W_ex=W+ex_coef1*(W_bar);

    bAX=X./(W_ex*H+eps);
    cj=bAX*(H');
    W_prev=W;
    W=max(options.paramepsi,(W_ex.*cj)./(rjc+eps));
   
    i=i+1; 
    if options.obj_compute==1
        time1=tic; 
        e(i)= KLobj(X,W,H);
        timeerr=timeerr + toc(time1);
    end
    t(i)= toc(cputime0)-timeerr;
    if options.display ==1 && options.obj_compute==1
        % only display at iteration i s.t. i = 2^k for some k
        if i-1 == 2^kdis
            fprintf('MU: iteration %4d fitting error: %1.4e \n',i-1,e(i));
            kdis = kdis+1;
        end
    elseif options.display ==1
        if i-1 == 2^kdis
            fprintf('MU: iteration %4d:',i-1);
            kdis = kdis+1;
        end
    end
end
end 
