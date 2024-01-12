clear all;close all;clc; 
rng(2023); % for reproducibility
%% ------------------------------------------------------------------------
% Loading data set
%%------------------------------------------------------------------------
 
% % Mary had a little lamb
% load('2_piano_Mary.mat') 
% options.WINDOWSIZE = 256;   % for Mary only
% options.fs = 16000;         % for Mary only
% length_x =  75200;          % for Mary only
% r = 4;                      % for Mary only
% options.lambda_tilde = 0.3; % for Mary only

% % Prelude from J.S. Bach
% load('3_prelude_JSB.mat')
% options.WINDOWSIZE = 256;     % for prelude
% options.fs = 11025;           % for prelude
% length_x = 330750;            % for prelude
% r = 16;                       % 16 for over-shoot
% options.lambda_tilde = 0.04;  % for prelude

% % Piano sample from Fevotte, Bertin and Durieux
load('1_Bertin.mat')
options.WINDOWSIZE =1024;   % for Bertin
options.fs = 22050;         % for Bertin
length_x = 345500;          % for Bertin
r = 7;                      % for Bertin
options.lambda_tilde=0.015; % for Bertin

% copy and size computation
X=V;
[m,n]=size(X);

% report relative errors: error divided by KL(X, all-one matrix*average(X))
colX = sum(X,2)/n+eps;
nX = X.*log(X./repmat(colX,1,n)+eps);
nX = sum(nX(:));


%% ------------------------------------------------------------------------
% Setup parameters
%%------------------------------------------------------------------------
% Common parameters
options.timemax = 1;
options.display = 0;
options.delta = 1;
options.accuracy_bisMethod = 10^-6;
options.maxiter = 10^4;
options.obj_compute = 1;

% parameters for disjoint KL NMF (Alg. from [Leplat et al., 2021])
options.K = r;
options.MAXITER=10^4;
options.epsi=10^-6; %check with options.accuracy_bisMethod 
options.alpha=0;
options.init_points = 1;
%% ------------------------------------------------------------------------
% Benchmarks
%%------------------------------------------------------------------------
e_min_min=1e16;
nb_Trials = 10;
compa_setup = 1; %0: fixed max number of iterations, 1: fixed max cpu time
MUe_struc = [];
MUe_noex_struc = [];
disj_mvkl_struc = [];
disj_mvklex_struc = [];
gen_struc = [];
for numtrial=1:nb_Trials
    fprintf('------------------------------------ \n');
    fprintf('Trial %4d over %4d running... \n',numtrial,nb_Trials); 

    % initialization for W and H
    W = rand(m,r); 
    H = rand(r,n);
    % scale initial point 
    WH=W*H;
    alpha=sqrt(sum(X(:))/sum(WH(:)));
    W=alpha*W; 
    H=alpha*H; 
    options.init.W=W; 
    options.init.H=H;
    
    % Call of MUe algorithm
    fprintf('  Call of MUe Alg. \n');
    options.no_ex = 0;
    [Wmv,Hmv,emv,tmv] = MUe_minvolKLNMF(X,r,options);
    % Relative errors:
    emv = emv/nX;
    MUe_struc.W{numtrial} = Wmv; MUe_struc.H{numtrial} = Hmv;
    MUe_struc.e{numtrial} = emv; MUe_struc.t{numtrial} = tmv;
    
    % Call of MU algorithm 
    fprintf('  Call of MU Alg. \n');
    options.no_ex = 1;
    [Wmv_noex,Hmv_noex,emv_noex,tmv_noex] = MUe_minvolKLNMF(X,r,options);
    % Relative errors:
    emv_noex = emv_noex/nX;
    MUe_noex_struc.W{numtrial} = Wmv_noex; MUe_noex_struc.H{numtrial} = Hmv_noex;
    MUe_noex_struc.e{numtrial} = emv_noex; MUe_noex_struc.t{numtrial} = tmv_noex;
    
    % Call ok disjoint min-vol KL NMF
    fprintf('  Call of Alg. from [Leplat et al., 2021]\n');
    [Wdismv, Hdismv, edismv, tdismv] = disjointconstraint_minvol_KLNMF(X,options);
    % Relative errors:
    edismv = edismv/nX;
    disj_mvkl_struc.W{numtrial} = Wdismv; disj_mvkl_struc.H{numtrial} = Hdismv;
    disj_mvkl_struc.e{numtrial} = edismv; disj_mvkl_struc.t{numtrial} = tdismv;

    % Call ok disjoint min-vol KL NMF with extrapolation
    fprintf('  Call of Alg. from [Leplat et al., 2021] with ex.\n');
    [Wdismv_ex, Hdismv_ex, edismv_ex, tdismv_ex] = disjointconstraint_minvol_KLNMFe(X,options);
    % Relative errors:
    edismv_ex = edismv_ex/nX;
    disj_mvklex_struc.W{numtrial} = Wdismv_ex; disj_mvklex_struc.H{numtrial} = Hdismv_ex;
    disj_mvklex_struc.e{numtrial} = edismv_ex; disj_mvklex_struc.t{numtrial} = tdismv_ex;

    
    e_min = min([min(emv),min(emv_noex),min(edismv),min(edismv_ex)]);
    gen_struc.e_mini{numtrial} = e_min;
    if e_min<e_min_min
        e_min_min = e_min;
    end

    fprintf('------------------------------------ \n');

end

if compa_setup == 0
    % computation of the averaged curves (for comparison w.r.t. iteration
    % number)
    MUe_e_mean = zeros(1,length(MUe_struc.e{1}));
    MUe_noex_e_mean = zeros(1,length(MUe_noex_struc.e{1}));
    disj_mvkl_e_mean = zeros(1,length(disj_mvkl_struc.e{1}));
    disj_mvklex_e_mean = zeros(1,length(disj_mvklex_struc.e{1}));
    for i=1:length(MUe_struc.e)
        MUe_e_mean = MUe_e_mean + (MUe_struc.e{i});
        MUe_noex_e_mean = MUe_noex_e_mean + (MUe_noex_struc.e{i});
        disj_mvkl_e_mean = disj_mvkl_e_mean + (disj_mvkl_struc.e{i});
        disj_mvklex_e_mean = disj_mvklex_e_mean + (disj_mvklex_struc.e{i});
    end
    MUe_e_mean = MUe_e_mean/length(MUe_struc.e);
    MUe_noex_e_mean = MUe_noex_e_mean/length(MUe_noex_struc.e);
    disj_mvkl_e_mean = disj_mvkl_e_mean/length(disj_mvkl_struc.e);
    disj_mvklex_e_mean = disj_mvklex_e_mean/length(disj_mvklex_struc.e);
    
    MUe_ttot = zeros(length(MUe_struc.t),1);
    MUe_noex_ttot = zeros(length(MUe_noex_struc.t),1);
    disj_mvkl_ttot = zeros(length(disj_mvkl_struc.t),1);
    disj_mvklex_ttot = zeros(length(disj_mvklex_struc.t),1);
    for i=1:length(MUe_struc.t)
        MUe_ttot(i) =  MUe_struc.t{i}(end);
        MUe_noex_ttot(i) =  MUe_noex_struc.t{i}(end);
        disj_mvkl_ttot(i) =  disj_mvkl_struc.t{i}(end);
        disj_mvklex_ttot(i) =  disj_mvklex_struc.t{i}(end);
    end

elseif compa_setup==1
    % % % Computation of the interpolations for each trial
    
    % find the first temporal point for each algo to ensure the good
    % interpolations
    maximum_t_1 = MUe_struc.t{1}(1); maximum_t_2 = MUe_noex_struc.t{1}(1);
    maximum_t_3 = disj_mvkl_struc.t{1}(1); maximum_t_4 = disj_mvklex_struc.t{1}(1);
    for ti=2:nb_Trials
        if MUe_struc.t{ti}(1)>maximum_t_1
            maximum_t_1 = MUe_struc.t{ti}(1);
        end
        if MUe_noex_struc.t{ti}(1)>maximum_t_2
            maximum_t_2 = MUe_noex_struc.t{ti}(1);
        end
        if disj_mvkl_struc.t{ti}(1)>maximum_t_3
            maximum_t_3 = disj_mvkl_struc.t{ti}(1);
        end
        if disj_mvklex_struc.t{ti}(1)>maximum_t_4
            maximum_t_4 = disj_mvklex_struc.t{ti}(1);
        end
    end
    % Create the query grids common to all trials for each algo
    temp_acc = 5*10^-3; % setup time step

    [Xq1,Yq1] = meshgrid(maximum_t_1:temp_acc:options.timemax);
    sizeXq1 = size(Xq1);
    tempo_grid1 = Xq1(1,:);

    [Xq2,Yq2] = meshgrid(maximum_t_2:temp_acc:options.timemax);
    sizeXq2 = size(Xq2);
    tempo_grid2 = Xq2(1,:);

    [Xq3,Yq3] = meshgrid(maximum_t_3:temp_acc:options.timemax);
    sizeXq3 = size(Xq3);
    tempo_grid3 = Xq3(1,:);

    [Xq4,Yq4] = meshgrid(maximum_t_4:temp_acc:options.timemax);
    sizeXq4 = size(Xq4);
    tempo_grid4 = Xq4(1,:);

    % Init. of arrays
    MUe_eInterpol = zeros(sizeXq1(1),nb_Trials);
    MUe_noex_eInterpol = zeros(sizeXq2(1),nb_Trials);
    disj_mvkl_eInterpol = zeros(sizeXq3(1),nb_Trials);
    disj_mvklex_eInterpol = zeros(sizeXq4(1),nb_Trials);
    
    % % Interpolation for all trials and for all algorithms
    for ti=1:nb_Trials
        % % MUe alg.
        % Create the meshgrid with existing point values
        [X,Y] = meshgrid(MUe_struc.t{ti});
        sizeGrid = size(X);
        V = repmat((MUe_struc.e{ti})',1,sizeGrid(2));
        % Compute the interpolation
        Vq = interp2(X,Y,V,Xq1,Yq1);
        MUe_eInterpol(:,ti) = Vq(:,1);

        % % MUe alg. without extrapolation
        [X,Y] = meshgrid(MUe_noex_struc.t{ti});
        sizeGrid = size(X);
        V = repmat((MUe_noex_struc.e{ti})',1,sizeGrid(2));
        % Compute the interpolation
        Vq = interp2(X,Y,V,Xq2,Yq2);
        MUe_noex_eInterpol(:,ti) = Vq(:,1);

        % % disjoint min-vol KL NMF
        [X,Y] = meshgrid(disj_mvkl_struc.t{ti});
        sizeGrid = size(X);
        V = repmat((disj_mvkl_struc.e{ti})',1,sizeGrid(2));
        % Compute the interpolation
        Vq = interp2(X,Y,V,Xq3,Yq3);
        disj_mvkl_eInterpol(:,ti) = Vq(:,1);

        % % disjoint min-vol KL NMF with extrapolation
        [X,Y] = meshgrid(disj_mvklex_struc.t{ti});
        sizeGrid = size(X);
        V = repmat((disj_mvklex_struc.e{ti})',1,sizeGrid(2));
        % Compute the interpolation
        Vq = interp2(X,Y,V,Xq4,Yq4);
        disj_mvklex_eInterpol(:,ti) = Vq(:,1);


    end
    
    % % Computation of the medians for interpolated curves
    MUe_eInterpol_median = median(MUe_eInterpol,2);
    MUe_noex_eInterpol_median = median(MUe_noex_eInterpol,2);
    disj_mvkl_eInterpol_median = median(disj_mvkl_eInterpol,2);
    disj_mvklex_eInterpol_median = median(disj_mvklex_eInterpol,2);

    % % Computation of the min for interpolated curves
    MUe_eInterpol_min = min(MUe_eInterpol,[],2);
    MUe_noex_eInterpol_min = min(MUe_noex_eInterpol,[],2);
    disj_mvkl_eInterpol_min = min(disj_mvkl_eInterpol,[],2);
    disj_mvklex_eInterpol_min = min(disj_mvklex_eInterpol,[],2);
    
end


%% ------------------------------------------------------------------------
% Post-processing
%%------------------------------------------------------------------------
close all
if compa_setup == 0
    figure;
    set(0, 'DefaultAxesFontSize', 18);
    set(0, 'DefaultLineLineWidth', 2);
    semilogy(MUe_e_mean-e_min_mins,'c','LineWidth',1.5);hold on; 
    semilogy(disj_mvkl_e_mean-e_min_mins,'r','LineWidth',1.5);hold on;
    semilogy(MUe_noex_e_mean-e_min_mins,'b','LineWidth',1.5);hold on; 
    semilogy(disj_mvklex_e_mean-e_min_mins,'k','LineWidth',1.5);hold on;
    ylabel('$e_{rel} - e_{rel,min}$','Interpreter','latex');
    xlabel('Iteration counter','Interpreter','latex');
    title([' Av. over ' num2str(nb_Trials)  ' runs - $\tilde{\lambda}=$ ' num2str(options.lambda_tilde)],'FontSize',18, 'Interpreter','latex')
    legend(['MUe - ' num2str(mean(MUe_ttot),'%2.2f') ' +/- ' num2str(std(MUe_ttot),'%2.2f') ' sec.' ],['[Leplat et al., 2021] - ' num2str(mean(disj_mvkl_ttot),'%2.2f') ' +/- ' num2str(std(disj_mvkl_ttot),'%2.2f') ' sec.' ],['MU - ' num2str(mean(MUe_noex_ttot),'%2.2f') ' +/- ' num2str(std(MUe_noex_ttot),'%2.2f') ' sec.'],['[Leplat et al., 2021]e - ' num2str(mean(disj_mvklex_ttot),'%2.2f') ' +/- ' num2str(std(disj_mvklex_ttot),'%2.2f') ' sec.' ]); 
    grid on;

 elseif compa_setup==1
    e_min_mins = 1*e_min_min;  %0.98 for Mary
    figure;
    set(0, 'DefaultAxesFontSize', 18);
    set(0, 'DefaultLineLineWidth', 2);
    MUe_curve_1 = MUe_eInterpol_median-e_min_mins;
    disj_mvkl_curve_1 = disj_mvkl_eInterpol_median-e_min_mins;
    MUe_noex_curve_1 = MUe_noex_eInterpol_median-e_min_mins;
    disj_mvklex_curve_1 = disj_mvklex_eInterpol_min-e_min_mins;
    star_idx = 2; % remove 2-3 first iterations to improve visibility
    semilogy(tempo_grid1(star_idx:end)',MUe_curve_1(star_idx:end),'c-.','LineWidth',1.5);hold on; 
    semilogy(tempo_grid3(star_idx:end)',disj_mvkl_curve_1(star_idx:end),'r--','LineWidth',1.5);hold on;
    semilogy(tempo_grid2(star_idx:end)',MUe_noex_curve_1(star_idx:end),'b','LineWidth',1.5);hold on; 
    semilogy(tempo_grid4(star_idx:end)',disj_mvklex_curve_1(star_idx:end),'k-x','LineWidth',1.5);hold on;
    ylabel('$e_{rel} - e_{rel,min}$','Interpreter','latex');
    xlabel('time (s.)','Interpreter','latex'); 
    title([' Median curves over ' num2str(nb_Trials)  ' runs - $\tilde{\lambda}=$ ' num2str(options.lambda_tilde)],'FontSize',18, 'Interpreter','latex')
    legend(['MUe'],['[Leplat et al., 2021] ' ],['MU' ],['[Leplat et al., 2021]e.' ]); 
    grid on;

    figure;
    set(0, 'DefaultAxesFontSize', 18);
    set(0, 'DefaultLineLineWidth', 2);
    MUe_curve_2 = MUe_eInterpol_min-e_min_mins;
    disj_mvkl_curve_2 = disj_mvkl_eInterpol_min-e_min_mins;
    MUe_noex_curve_2 = MUe_noex_eInterpol_min-e_min_mins;
    disj_mvklex_curve_2 = disj_mvklex_eInterpol_min-e_min_mins;
    star_idx = 2; % remove 2-3 first iterations to improve visibility
    semilogy(tempo_grid1(star_idx:end)',MUe_curve_2(star_idx:end),'c-.','LineWidth',1.5);hold on; 
    semilogy(tempo_grid3(star_idx:end)',disj_mvkl_curve_2(star_idx:end),'r--','LineWidth',1.5);hold on;
    semilogy(tempo_grid2(star_idx:end)',MUe_noex_curve_2(star_idx:end),'b','LineWidth',1.5);hold on; 
    semilogy(tempo_grid4(star_idx:end)',disj_mvklex_curve_2(star_idx:end),'k-x','LineWidth',1.5);hold on;
    ylabel('$e_{rel} - e_{rel,min}$','Interpreter','latex');
    xlabel('time (s.)', 'Interpreter','latex')
    title([' Min. curves over ' num2str(nb_Trials)  ' runs - $\tilde{\lambda}=$ ' num2str(options.lambda_tilde)],'FontSize',18, 'Interpreter','latex')
    legend(['MUe'],['[Leplat et al., 2021] ' ],['MU' ],['[Leplat et al., 2021]e.' ]); 
    grid on;
end 

%% ------------------------------------------------------------------------
% Visualize some audio sources
%%------------------------------------------------------------------------
% PLOT Basis Vectors

T = n;
set(0, 'DefaultAxesFontSize', 12);
time=linspace(0,length_x/options.fs,T);
freq=linspace(0,options.fs/2,options.WINDOWSIZE/2+1)/1000;
nCourbes = options.K;
couleurs = hsv(nCourbes);
% choose which W ?
W = Wdismv_ex; %Wdismv, Wmv, Wdismv_ex
figure2=figure;
for i=1:options.K
    plot((2*i-1)*max(max(W))+(1-W(:,i)),freq,'LineWidth', 1,'color', couleurs(i,:),'DisplayName',['source estimate ' num2str(i)]);
    hold on
end
ylabel('Frequency [kHz]','FontSize',14, 'Interpreter','latex')
xlabel('$W$ columns ','FontSize',14, 'Interpreter','latex')
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
% saveas(figure2,'Basis_Vectors_latex.fig')

% Plot activations
figure3=figure;
% choose which H ?
H = Hdismv_ex; % Hdismv, Hmv, Hdismv_ex
for i=1:options.K
     plot(time, (i-1)*max(max(H))+(H(i,:)),'LineWidth', 1,'color', couleurs(i,:),'DisplayName',['source estimate ' num2str(i)]);
    hold on
end
ylabel('Activations','FontSize',14, 'Interpreter','latex')
xlabel('Time[s]','FontSize',14, 'Interpreter','latex')
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
% saveas(figure3,'Activations_latex.fig')


%% ------------------------------------------------------------------------
% Backlog of codes
%%------------------------------------------------------------------------
% Sanity check
% figure; semilogy(X(1,:),V(:,1),'r','LineWidth',1.5);hold on;
% semilogy(Xq(1,:),Vq(:,1),'.k','LineWidth',1.5);hold on;
% legend('original','interpolated');
% grid on;
