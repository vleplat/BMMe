function [W, H, lossfun, t] = betaNMF_baseline(V,options)

% This function tackles the NMF problem using the beta-divergence
% (specified in options) via multiplicative updates. 

%INPUTS:
%       V:                 Spectrogram of x (input audio signal) 
%       options.K:         Number of sources (determined prior to the
%       factorization)
%       options.MAXITER:   Maximum number of iterations for updates of W and H
%       options.beta:      Choice of Beta-Divergence computation
%       options.alpha:     Parameter of the MU to control the sparsity of H
%       options.init       Use this algorithm for initialisation if
%                          options.init=1 and do not display anything, 
%                          options.init==0 otherwise. 


%OUTPUTS:
%       W: Matrix dictionnaries
%       H: Matrix activations
%       lossfun: value of the loss function f(W,H) after MAXITER iterations
%       t: time in seconds to return a solution

% % Loading parameters
K=options.K;
MAXITER=options.MAXITER;
beta=options.beta;
alpha=options.alpha;
init=options.init;

% % Algorithm Beta-NMF
disp(' ->baseline-MU NMF algorithm')
tic
F = size(V,1);
T = size(V,2);
addpath('./Utils');

% %  Initialization for W and H

% Random initialization

if(init==0)
    disp(' ->Random Initialization for W and H')
end
W = 1+rand(F, K);
H = 1+rand(K, T);



% % array to save the value of the loss function
lossfunsave = zeros(MAXITER,1);
if(init==0)
    fprintf(' ->The initial value for betadivergence is %0.2f \n', betaDiv(V+eps,W*H+eps,beta));
end

if(init==0)
    ani1=animatedline('Color','r','Marker','o');
    title(['Evolution of objective function - $\beta$ = ' num2str(beta)],'FontSize',12, 'Interpreter','latex')
    xlabel('iteration','FontSize',12, 'Interpreter','latex')
end

% % Optimization loop
for iter=1:MAXITER 
    
    %  % update maxtrix W ("dictionaries")
    if(beta==0)
        W = W .* (((((W*H).^(beta-2)).*V)*H')./((W*H).^(beta-1)*H'+eps)).^(1/2);
    else
        W = W .* ((((W*H).^(beta-2)).*V)*H')./((W*H).^(beta-1)*H'+eps);
    end
    
    W = max(W,eps);

    
    % % update matrix  H Coefficients ("activations")
    if(beta==0)
        H = (H .* ((W'*(((W*H).^(beta-2)).*V))./(W'*(W*H).^(beta-1)+eps)).^(1/2)).^(1+alpha);
    else
        H = (H .* (W'*(((W*H).^(beta-2)).*V))./(W'*(W*H).^(beta-1)+eps)).^(1+alpha);
    end
    H = max(H,eps);

    % % In practice, we impose ||Wk||1 = 1 
    sumW = sum(W)+eps;
    W = W*diag(1./sumW);
    H = diag(sumW)*H;

    % % Compute the loss function = Beta-Divergence 
    lossfunsave(iter) = betaDiv(V+eps,W*H+eps,beta);

    % Drawing
    if(init==0)
        addpoints(ani1,iter,log10(lossfunsave(iter))) ;
        drawnow
    end
    
end%End of optimization loop

% % Loss function graph saving
if (options.init==0)
    figure1=figure;
    semilogy(lossfunsave,'-o');
    title(['Evolution of objective function - $\beta$ = ' num2str(beta)],'FontSize',12, 'Interpreter','latex');
    xlabel('iteration','FontSize',12, 'Interpreter','latex');
    saveas(figure1,'./Graphs/lossfun.png')
    close figure 1
end

% % Save lossfun
lossfun=lossfunsave;
% % Final Value for loss function
if(init==0)
    fprintf(' ->The final value for betadivergence is %0.2f \n', betaDiv(V+eps,W*H+eps,beta));
end
t = toc;
end %EOF
