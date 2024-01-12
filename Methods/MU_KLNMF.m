%% Multiplicative Update for solving KL NMF
% Input: X, r, options
% options include {display, init.W, init.H, maxiter, timemax, paraepsi,obj_compute}
% Output: W, H (the factors), e (the objective sequence) and t (the sequence of
% running time)
%
% if options.obj_compute = 1 then the output (default value) then the output e would be 
% the sequence of the objective values;  otherwise, the output e would be [].
%
% written by LTK Hien
% Latest update: December 2023
function [W,H,e,t] = MU(X,r,options) 
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

W=options.init.W ;
H=options.init.H ;

i = 1; 
t(1) = toc(cputime0);
e=[];
timeerr=0;
if options.obj_compute==1
    time1=tic;
    e(1)= KLobj(X,W,H);
    timeerr=toc(time1); % to remove the time of finding the objective function
end
while i <= options.maxiter && t(i) < options.timemax  
   
    %update H
    Wt=W';
    rj=sum(W)';
    rjc=repmat(rj,1,n);
    bAx=X./(W*H+eps);
    cj=Wt*bAx; 
    H=max(options.paramepsi,(H.*cj)./(rjc+eps));
     
    %update W
   
    rj=sum(H,2)';
    rjc=repmat(rj,m,1); 
    bAX=X./(W*H+eps);
    cj=bAX*(H');
    W=max(options.paramepsi,(W.*cj)./(rjc+eps));
   
    i=i+1; 
    if options.obj_compute==1
        time1=tic; 
        e(i)= KLobj(X,W,H);
        timeerr=timeerr + toc(time1);
    end
    t(i)= toc(cputime0)-timeerr;
    if  options.display ==1 && options.obj_compute==1
          fprintf('MU: iteration %4d fitting error: %1.4e \n',i,e(i));     
    end
end
end


