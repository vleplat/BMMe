% numerical experiment for MUe
close all; clear all; clc;
% Parameters for the experiment:
numdatasets = 4; % number of data sets tested
numinit = 10; % number of random initializations used
options.maxiter = 100; % number of iterations
options.beta = 1.5; % value of beeta tested
% Make sure 100 iterations are performed, do not use other stopping
% conditions:
options.accuracy = 0;
options.timemax = Inf;
% lower bound for the entries of W and J:
options.epsilon = 2^(-52);
fprintf('Running HSI experiment: %2.0f data sets, %2.0f initializations.\n',...
    numdatasets,numinit);
fprintf('This should take of the order of %2.0f minutes.\n', 2*numdatasets*options.maxiter*numinit/60);
for ida = 1 : numdatasets
    if ida == 1
        load SanDiego % ~ 3 second for 1 iteration
        r = 8;
    elseif ida == 2
        load Urban % ~ 2 second for 1 iteration
        r = 6;
    elseif ida == 3
        load Cuprite; % ~ 1 second for 1 iteration
        r = 20;
    elseif ida == 4
        load Pines; % ~ 1 second for 1 iteration
        r = 16;
    end
    [m,n] = size(X);
    rng(2023); % control random seed
    for i = 1 : numinit
        fprintf('Running data set %2.0f, init %2.0f...\n',ida,i);
        % Initialization
        W0 = max(options.epsilon,rand(m,r));
        H0 = max(options.epsilon,rand(r,n));
        %% Beta = 3/2 -- KL divergence
        options.H = MUbeta(X,W0,H0,options.beta,options.epsilon);
        options.W = MUbeta(X',options.H',W0',options.beta,options.epsilon);
        options.W = options.W';
        [W,H,e,t] = betaNMF(X,r,options);
        options.extrapol = @extrapol_nesterov;
        [We,He,ee,te] = betaNMFextra(X,r,options);
        etot(ida,i,:) = e;
        ttot(ida,i,:) = t;
        eetot(ida,i,:) = ee;
        tetot(ida,i,:) = te;
        disp('Done.');
    end
end
% Display results for each data set
set(0, 'DefaultAxesFontSize', 22);
set(0, 'DefaultLineLineWidth', 2);
for ida = 1 : numdatasets
    if ida == 1
        load SanDiego % ~ 3 second for 1 iteration
        r = 8;
    elseif ida == 2
        load Urban % ~ 2 second for 1 iteration
        r = 6;
    elseif ida == 3
        load Cuprite; % ~ 1 second for 1 iteration
        r = 20;
    elseif ida == 4
        load Pines; % ~ 1 second for 1 iteration
        r = 16;
    end
    ida
    e = reshape(etot(ida,:,:),numinit,options.maxiter);
    ee = reshape(eetot(ida,:,:),numinit,options.maxiter);
    T = (ee < e(end));
    R(ida,:) = 100-sum(T')
    emin = min(min(e(:)),min(ee(:)));
    figure;
    dist = betadiv(X,zeros(size(X)),options.beta);
    semilogy((median(e)-emin)/dist); hold on;
    semilogy((median(ee)-emin)/dist,'--');
    grid on;
    legend('MU', 'MUe', 'Interpreter', 'latex');
    xlabel('Iterations', 'Interpreter', 'latex','FontSize',22);
    ylabel('$\frac{D_{3/2}(X,WH)}{D_{3/2}(X,0)}-e_{min}$', 'Interpreter', 'latex','FontSize',22);
end
