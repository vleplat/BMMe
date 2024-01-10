function MyProject_NMF(options)

%% BLOC 0 : INITIALIZATION and Data Loading
disp('Start Data loading')
% % Load parameter
if nargin == 0
    [options]=parameters; % default parameters 
end
% % Load audio file
% Add functions from Utils directory
addpath('./Utils_lib');
[x,options.fs]=loadAudioFile(options);
if(options.timecut==1)
    x = time_cut(x,options);
    options.metric_eval=0;%disable metric eval
end

% % Audio signal Playing
%sound(x,options.fs);

% % Suppress previous sources extracted and graphs
myFolder = './Extractions/';
extensionFile='*.wav';
emptyDir(myFolder,extensionFile);
myFolder = './Graphs/';
extensionFile='*.png';
emptyDir(myFolder,extensionFile);
extensionFile='*.fig';
emptyDir(myFolder,extensionFile);
myFolder = './BenchMark/';
extensionFile='*.fig';
emptyDir(myFolder,extensionFile);
extensionFile='*.png';
emptyDir(myFolder,extensionFile);
myFolder = './BenchMark/Extractions/baseline/';
extensionFile='*.wav';
emptyDir(myFolder,extensionFile);
myFolder = './BenchMark/Extractions/minvol/';
extensionFile='*.wav';
emptyDir(myFolder,extensionFile);
myFolder = './BenchMark/Extractions/sparse/';
extensionFile='*.wav';
emptyDir(myFolder,extensionFile);

disp('End of data loading')
fprintf('\n');
%% BLOC 1: Short Time Fourier Transform

disp('Start STFT computation')
% windows generation
[options.windows]= windowsGeneration(options.winselection,options.WINDOWSIZE);
% STFT computation
X = STFT(x, options.WINDOWSIZE, options.windows, options.hopsize);
V=abs(X);
T = size(V,2);

disp('End of STFT computation')
fprintf('\n');

%% BLOC 2: NMF
disp('Start NMF computation')
if(strcmp(options.algo,'minvol'))
    if(options.beta==0)
        [W, H, lossfunction, cputime_logdetIS,grad] = betaNMF_logdetIS(V,options);
    elseif(options.beta==1)
        [W, H, lossfunction, cputime_logdetKL] = betaNMF_logdetKL(V,options);
    end
elseif(strcmp(options.algo,'disjointconstraint_minvol_KLNMF'))
    [W, H, lossfunction, cputime_logdetKL] = disjointconstraint_minvol_KLNMF(V,options);
elseif(strcmp(options.algo,'hypersphericstructured_minvol_KLNMF'))
    [W, H, lossfunction, cputime_logdetKL] = hypersphericstructured_minvol_KLNMF(V,options);
elseif(strcmp(options.algo,'baseline'))
    [W, H, lossfunction, cputime_b] = betaNMF_baseline(V,options);
    options.BM=0;
elseif(strcmp(options.algo,'sparse'))
    params.sparsity=options.sparsity;
    params.max_iter = options.MAXITER;
        if (options.beta==1)
            params.cf = 'kl';
        elseif(options.beta==0)
            params.cf = 'is';
        else
            params.cf = 'ed';
        end
    params.r=options.K;
    [W, H, objective, cputime_s] = sparse_nmf(V, params);
    lossfunction=objective.cost;
    options.BM=0;
end

% % Final Value for loss function
fprintf(' ->The final value for loss function (beta-divergence and penalty term) is %0.2f \n', lossfunction(end)); 

disp('End of NMF computation')
fprintf('\n');

%% BLOC 3: POSTPROCESSING
time=linspace(0,length(x)/options.fs,T);
freq=linspace(0,options.fs/2,options.WINDOWSIZE/2+1);
if (options.focus==1)
    
    % Check if the intervals boundary are correct
    if(options.time_interval(1)<time(1) || options.time_interval(2)>time(end))
        warning('Time interval specified not correct, the time selection is disabled');
        tmin=time(1);
        tmax=time(end);
        time_frame_min=1;
        time_frame_max=T;
  
    elseif(options.time_interval(1)>=options.time_interval(2))
        warning('Time interval specified not correct, the time selection is disabled');
        tmin=time(1);
        tmax=time(end);
        time_frame_min=1;
        time_frame_max=T;
    
    elseif(abs(options.time_interval(1)-options.time_interval(2))<=1/options.fs)
        warning('Time interval lower than the sampling period, the time selection is disabled');
        tmin=time(1);
        tmax=time(end);
        time_frame_min=1;
        time_frame_max=T;
        
    elseif(length(options.time_interval)>2 ||length(options.time_interval)<2 )
        warning('Time interval specified not correct, the time selection is disabled');
        tmin=time(1);
        tmax=time(end);
        time_frame_min=1;
        time_frame_max=T;
    else
        % Find the closest bounds
        [~,idx]=min(abs(time-options.time_interval(1)));
        tmin=time(idx);
        time_frame_min=idx;
        [~,idx]=min(abs(time-options.time_interval(2)));
        tmax=time(idx);
        time_frame_max=idx;
        if(time_frame_min==time_frame_max)
            if(time_frame_max<T)
                time_frame_max=time_frame_max+1;
                tmax=time(time_frame_max);
            else
                if(time_frame_min>1)
                    time_frame_min=time_frame_min-1;
                    tmin=time(time_frame_min);
                else
                    warning('Time interval specified not correct, the time selection is disabled');
                    tmin=time(1);
                    tmax=time(end);
                    time_frame_min=1;
                    time_frame_max=T;
                end
            end
        end
    end

else
    tmin=time(1);
    tmax=time(end);
    time_frame_min=1;
    time_frame_max=T;

end
if (options.BM==0)
    disp('Start POSTPROCESSING')
    % % Source estimates generation by ISTFT
    [x_Source,Mask] = temporal_Reconstruction(W,H,angle(X),options,V);

    % % Estimated Sources saving
    % Specify the destination folder.
    myFolder = './Extractions/';
    saveEstimates(x_Source,options,myFolder)

    % % Save variables from the workspace
    save variables.mat

   
    % % Graph generation - Latex Format
    disp(' ->Graphs generation');
    % Spectrograms and Oginal Time-domain signal
    length(x);
    t=(1:length(x))*1/options.fs;
    figure1=figure;
    subplot(2,1,1);
    plot(t,x)
    title('Input Audio Signal $x$','FontSize',14, 'Interpreter','latex')
    ylabel('[-]','FontSize',14, 'Interpreter','latex')
    xlabel('Time[s]','FontSize',14, 'Interpreter','latex')
    xlim([tmin tmax])

    subplot(2,1,2);
    imagesc(db(V(:,time_frame_min:time_frame_max)))%db()=10*log10()
    set(gca,'YDir','normal') %Flip ordinate
    title('Amplitude Spectrogram','FontSize',14, 'Interpreter','latex')
    ylabel('Frequency $f$','FontSize',14, 'Interpreter','latex')
    xlabel('Time $n$','FontSize',14, 'Interpreter','latex')
    saveas(figure1,'./Graphs/Input_latex.fig')
    
    % % Comparison Spectrograms (Orignal VS Reconstruction)
    figure2=figure;
    subplot(2,1,1);
    imagesc(db(V(:,time_frame_min:time_frame_max)))
    set(gca,'YDir','normal') %Flip ordinate
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
    title('Spectrogram: Representation Frequency vs Time of the Input Audio Signal','FontSize',14, 'Interpreter','latex')
    ylabel('Frequency $f$','FontSize',14, 'Interpreter','latex')
    xlabel('Time $n$','FontSize',14, 'Interpreter','latex')
    colormap bone

    subplot(2,1,2);
    V_chap=W*H;
    imagesc(db(V_chap(:,time_frame_min:time_frame_max)))
    set(gca,'YDir','normal') %Flip ordinate
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
    title('Spectrogram: Representation Frequency vs Time of the Input Audio Signal','FontSize',14, 'Interpreter','latex')
    ylabel('Frequency $f$','FontSize',14, 'Interpreter','latex')
    xlabel('Time $n$','FontSize',14, 'Interpreter','latex')
    saveas(figure2,'./Graphs/Spectro_compa.fig')

    
    % Objective function
    figure3=figure;
    semilogy(lossfunction,'-o');
    title(['Evolution of objective function - $\beta$ = ' num2str(options.beta)],'FontSize',14, 'Interpreter','latex');
    xlabel('iteration','FontSize',14, 'Interpreter','latex');
    saveas(figure3,'./Graphs/lossfun_latex.fig')

    % PLOT Basis Vectors
    set(0, 'DefaultAxesFontSize', 12);
    time=linspace(0,length(x)/options.fs,T);
    freq=linspace(0,options.fs/2,options.WINDOWSIZE/2+1)/1000;
    nCourbes = options.K;
    couleurs = hsv(nCourbes);
    figure4=figure;
    for i=1:options.K
        plot((2*i-1)*max(max(W))+(1-W(:,i)),freq,'LineWidth', 1,'color', couleurs(i,:),'DisplayName',['source estimate ' num2str(i)]);
        hold on
    end
    ylabel('Frequency [kHz]','FontSize',14, 'Interpreter','latex')
    xlabel('$W$ columns ','FontSize',14, 'Interpreter','latex')
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
    saveas(figure4,'./Graphs/Basis_Vectors_latex.fig')

    % Plot activations
    figure5=figure;
    for i=1:options.K
         plot(time, (i-1)*max(max(H))+(H(i,:)),'LineWidth', 1,'color', couleurs(i,:),'DisplayName',['source estimate ' num2str(i)]);
        hold on
    end
    ylabel('Activations','FontSize',14, 'Interpreter','latex')
    xlabel('Time[s]','FontSize',14, 'Interpreter','latex')
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
    xlim([tmin tmax])
    saveas(figure5,'./Graphs/Activations_latex.fig')
    
    % % Comparison Masking coefficients
    figure6=figure;
    for i=1:options.K
        subplot(options.K,1,i);
        clims = [0 1];
        imagesc(Mask(:,time_frame_min:time_frame_max,i),clims)
        title(['Masking coeffcients for source estimate ' num2str(i)],'FontSize',10, 'Interpreter','latex')
        set(gca,'YDir','normal') %Flip ordinate
        set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
        set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
        ylabel('Frequency $f$','FontSize',8, 'Interpreter','latex')
        xlabel('Time $n$','FontSize',10, 'Interpreter','latex')
        colormap bone
    end
    saveas(figure6,'./Graphs/Masks.fig')

%     figure7=figure;
%     for i=1:options.K
%         subplot(options.K,1,i);
%         clims = [0 1];
%         imagesc(Mask(:,time_frame_min:time_frame_max,i),clims)
%         set(gca,'YDir','normal') %Flip ordinate
%         set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
%         set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
%         colormap bone
%     end
%     saveas(figure7,'./Graphs/Masks_nolabels.fig')
    
    figure8=figure;
    for i=1:options.K
        subplot(options.K,1,i);
        plot(freq,log(W(:,i)));%title('Plot of W(:,13) in log scale')
        title(['Plot of W(:,' num2str(i),') in log scale'],'FontSize',10, 'Interpreter','latex')
    end
    saveas(figure8,'./Graphs/basis_vectors_logscale.fig')

    disp('End of POSTPROCESSING')
end

%% BLOC 4: BENCHMARK
if (options.BM==1)
    disp('Start BENCHMARK')
    %remove previous wav files 
    myFolder = './BenchMark/Extractions/minvol/';
    extensionFile='*.wav';
    emptyDir(myFolder,extensionFile);
    myFolder = './BenchMark/Extractions/baseline/';
    emptyDir(myFolder,extensionFile);
    myFolder = './BenchMark/Extractions/sparse/';
    emptyDir(myFolder,extensionFile);
    % % BASELINE_NMF
    cond1=length(options.BMfun)==2;
    cond2=length(options.BMfun)==1 && strcmp(options.BMfun(1),'baseline');
    cond3=length(options.BMfun)==1 && strcmp(options.BMfun(1),'sparse');
    if(cond1 || cond2)
        options.init=1;
        [W_b, H_b, lossfunctionb, cputime_b] = betaNMF_baseline(V,options);
    end
    % % SPARSE-NMF
    if(cond1 || cond3)
        % parameters setting
        params.max_iter = options.MAXITER;
        if (options.beta==1)
            params.cf = 'kl';
        elseif(options.beta==0)
            params.cf = 'is';
        else
            params.cf = 'ed';
        end
        params.sparsity=10;
        params.r=options.K;
        % Call of NMF
        [W_s, H_s, objective,cputime_s] = sparse_nmf(V, params);
    end
    % % Assignement problem 
    if(cond1)
        [W,H,W_b,H_b]=matching(W,H,W_b,H_b);
        [W,H,W_s,H_s]=matching(W,H,W_s,H_s);
    elseif(cond2)
        [W,H,W_b,H_b]=matching(W,H,W_b,H_b);
    elseif(cond3)
        [W,H,W_s,H_s]=matching(W,H,W_s,H_s);
    end
    
%     % % Save variables from the workspace
      save variables.mat
    
    
    % % Comparison of results
    % Objective functions
    if(cond1)
        ymin=min([min(lossfunction) min(objective.cost) min(lossfunctionb)]);
        ymax=max([max(lossfunction) max(objective.cost) max(lossfunctionb)]);
        
        figurelfmin=figure;
        semilogy(lossfunction,'-o');
        if(options.beta==1)
            title('min-vol KL-NMF objective function','FontSize',20, 'Interpreter','latex');
        elseif(options.beta==0)
            title('min-vol IS-NMF objective function','FontSize',20, 'Interpreter','latex');
        end
        xlabel('iteration','FontSize',20, 'Interpreter','latex');
        ylim([ymin ymax])
        saveas(figurelfmin,'./BenchMark/lossfun_minvol.fig')
        saveas(figurelfmin,'./BenchMark/lossfun_minvol.png')
        
        figurelfs=figure;
        semilogy(objective.cost,'-o');
        if(options.beta==1)
            title('sparse KL-NMF objective function','FontSize',20, 'Interpreter','latex');
        elseif(options.beta==0)
            title('sparse IS-NMF objective function','FontSize',20, 'Interpreter','latex');
        end
        xlabel('iteration','FontSize',20, 'Interpreter','latex');
        ylim([ymin ymax])
        saveas(figurelfs,'./BenchMark/lossfun_sparse.fig')
        saveas(figurelfs,'./BenchMark/lossfun_sparse.png')
        
        figurelfb=figure;
        semilogy(lossfunctionb,'-o');
        if(options.beta==1)
            title('baseline KL-NMF objective function','FontSize',20, 'Interpreter','latex');
        elseif(options.beta==0)
            title('baseline IS-NMF objective function','FontSize',20, 'Interpreter','latex');
        end
        xlabel('iteration','FontSize',20, 'Interpreter','latex');
        ylim([ymin ymax])
        saveas(figurelfb,'./BenchMark/lossfun_baseline.fig')
        saveas(figurelfb,'./BenchMark/lossfun_baseline.png')
        
    elseif(cond2)
        ymin=min([min(lossfunction) min(lossfunctionb)]);
        ymax=max([max(lossfunction) max(lossfunctionb)]);
        figurelfmin=figure;
        semilogy(lossfunction,'-o');
        if(options.beta==1)
            title('min-vol KL-NMF objective function','FontSize',20, 'Interpreter','latex');
        elseif(options.beta==0)
            title('min-vol IS-NMF objective function','FontSize',20, 'Interpreter','latex');
        end
        xlabel('iteration','FontSize',20, 'Interpreter','latex');
        ylim([ymin ymax])
        saveas(figurelfmin,'./BenchMark/lossfun_minvol.fig')
        saveas(figurelfmin,'./BenchMark/lossfun_minvol.png')
        
        figurelfb=figure;
        semilogy(lossfunctionb,'-o');
        if(options.beta==1)
            title('baseline KL-NMF objective function','FontSize',20, 'Interpreter','latex');
        elseif(options.beta==0)
            title('baseline IS-NMF objective function','FontSize',20, 'Interpreter','latex');
        end
        xlabel('iteration','FontSize',14, 'Interpreter','latex');
        ylim([ymin ymax])
        saveas(figurelfb,'./BenchMark/lossfun_baseline.fig')
        saveas(figurelfb,'./BenchMark/lossfun_baseline.png')
    elseif(cond3)
        ymin=min([min(lossfunction) min(objective.cost)]);
        ymax=max([max(lossfunction) max(objective.cost)]);
        
        figurelfmin=figure;
        semilogy(lossfunction,'-o');
        if(options.beta==1)
            title('min-vol KL-NMF objective function','FontSize',20, 'Interpreter','latex');
        elseif(options.beta==0)
            title('min-vol IS-NMF objective function','FontSize',20, 'Interpreter','latex');
        end
        xlabel('iteration','FontSize',20, 'Interpreter','latex');
        ylim([ymin ymax])
        saveas(figurelfmin,'./BenchMark/lossfun_minvol.fig')
        saveas(figurelfmin,'./BenchMark/lossfun_minvol.png')
        
        figurelfs=figure;
        semilogy(objective.cost,'-o');
        if(options.beta==1)
            title('sparse KL-NMF objective function','FontSize',14, 'Interpreter','latex');
        elseif(options.beta==0)
            title('sparse IS-NMF objective function','FontSize',14, 'Interpreter','latex');
        end
        xlabel('iteration','FontSize',20, 'Interpreter','latex');
        ylim([ymin ymax])
        saveas(figurelfs,'./BenchMark/lossfun_sparse.fig')
        saveas(figurelfs,'./BenchMark/lossfun_sparse.png')
    end
    font_s_b=60;
    % PLOT Basis Vectors
    set(0, 'DefaultAxesFontSize', 12);
    time=linspace(0,length(x)/options.fs,T);
    freq=linspace(0,options.fs/2,options.WINDOWSIZE/2+1)/1000;
    nCourbes = options.K;
    couleurs = hsv(nCourbes);
    figure1=figure;
    for i=1:options.K
        plot((2*i-1)*max(max(W))+(1-W(:,i)),freq,'LineWidth', 1,'color', couleurs(i,:),'DisplayName',['source estimate ' num2str(i)]);
        hold on
    end
    if(options.beta==1)
        title('min-vol KL-NMF','FontSize',font_s_b, 'Interpreter','latex')
    elseif(options.beta==0)
        title('min-vol IS-NMF','FontSize',font_s_b, 'Interpreter','latex')
    end
    ylabel('Frequency [kHz]','FontSize',font_s_b, 'Interpreter','latex')
    xlabel('$W$ columns ','FontSize',font_s_b, 'Interpreter','latex')
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
    saveas(figure1,'./BenchMark/W_col_minvol.fig')
    saveas(figure1,'./BenchMark/W_col_minvol.png')
    
    if(cond1 || cond3)
        figure2=figure;
        for i=1:options.K
            plot((2*i-1)*max(max(W_s))+(1-W_s(:,i)),freq,'LineWidth', 1,'color', couleurs(i,:),'DisplayName',['source estimate ' num2str(i)]);
            hold on
        end
        if(options.beta==1)
            title('sparse KL-NMF','FontSize',font_s_b, 'Interpreter','latex')
        elseif(options.beta==0)
            title('sparse IS-NMF','FontSize',font_s_b, 'Interpreter','latex')
        end
        ylabel('Frequency [kHz]','FontSize',font_s_b, 'Interpreter','latex')
        xlabel('$W$ columns ','FontSize',font_s_b, 'Interpreter','latex')
        set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
        saveas(figure2,'./BenchMark/W_col_sparse.fig')
        saveas(figure2,'./BenchMark/W_col_sparse.png')
    end
    if(cond1 || cond2)
        figure3=figure;
        for i=1:options.K
            plot((2*i-1)*max(max(W_b))+(1-W_b(:,i)),freq,'LineWidth', 1,'color', couleurs(i,:),'DisplayName',['source estimate ' num2str(i)]);
            hold on
        end
        if(options.beta==1)
            title('baseline KL-NMF','FontSize',font_s_b, 'Interpreter','latex')
        elseif(options.beta==0)
            title('baseline IS-NMF','FontSize',font_s_b, 'Interpreter','latex')
        end
        ylabel('Frequency [kHz]','FontSize',font_s_b, 'Interpreter','latex')
        xlabel('$W$ columns ','FontSize',font_s_b, 'Interpreter','latex')
        set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
        saveas(figure3,'./BenchMark/W_col_baseline.fig')
        saveas(figure3,'./BenchMark/W_col_baseline.png')
    end
    font_s_a=60;
    % Plot activations
    figure4=figure;
    for i=1:options.K
         plot(time, (i-1)*max(max(H))+(H(i,:)),'LineWidth', 1,'color', couleurs(i,:),'DisplayName',['source estimate ' num2str(i)]);
        hold on
    end
    if(options.beta==1)
        title('min-vol KL-NMF','FontSize',font_s_a, 'Interpreter','latex')
    elseif(options.beta==0)
        title('min-vol IS-NMF','FontSize',font_s_a, 'Interpreter','latex')
    end
    ylabel('Activations','FontSize',font_s_a, 'Interpreter','latex')
    xlabel('Time[s]','FontSize',font_s_a, 'Interpreter','latex')
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
    xlim([tmin tmax])
    saveas(figure4,'./BenchMark/H_rows_minvol.fig')
    saveas(figure4,'./BenchMark/H_rows_minvol.png')
    
    if(cond1 || cond3)
        figure5=figure;
        for i=1:options.K
             plot(time, (i-1)*max(max(H_s))+(H_s(i,:)),'LineWidth', 1,'color', couleurs(i,:),'DisplayName',['source estimate ' num2str(i)]);
            hold on
        end
        if(options.beta==1)
            title('sparse KL-NMF','FontSize',font_s_a, 'Interpreter','latex')
        elseif(options.beta==0)
            title('sparse IS-NMF','FontSize',font_s_a, 'Interpreter','latex')
        end
        ylabel('Activations','FontSize',font_s_a, 'Interpreter','latex')
        xlabel('Time[s]','FontSize',font_s_a, 'Interpreter','latex')
        set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
        xlim([tmin tmax])
        saveas(figure5,'./BenchMark/H_rows_sparse.fig')
        saveas(figure5,'./BenchMark/H_rows_sparse.png')
    end
    if(cond1 || cond2)
        figure6=figure;
        for i=1:options.K
             plot(time, (i-1)*max(max(H_b))+(H_b(i,:)),'LineWidth', 1,'color', couleurs(i,:),'DisplayName',['source estimate ' num2str(i)]);
            hold on
        end
        if(options.beta==1)
            title('baseline KL-NMF','FontSize',font_s_a, 'Interpreter','latex')
        elseif(options.beta==0)
            title('baseline IS-NMF','FontSize',font_s_a, 'Interpreter','latex')
        end
        ylabel('Activations','FontSize',font_s_a, 'Interpreter','latex')
        xlabel('Time[s]','FontSize',font_s_a, 'Interpreter','latex')
        set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
        xlim([tmin tmax])
        saveas(figure6,'./BenchMark/H_rows_baseline.fig')
        saveas(figure6,'./BenchMark/H_rows_baseline.png')
    end
    
    % % Mask coefficients
    % Source estimates generation by ISTFT
    [x_Source,Mask] = temporal_Reconstruction(W,H,angle(X),options,V);
    myFolder = './BenchMark/Extractions/minvol/';
    saveEstimates(x_Source,options,myFolder)
    if(cond1)
        [x_Source_s,Mask_s] = temporal_Reconstruction(W_s,H_s,angle(X),options,V);
        myFolder = './BenchMark/Extractions/sparse/';
        saveEstimates(x_Source_s,options,myFolder)
        [x_Source_b,Mask_b] = temporal_Reconstruction(W_b,H_b,angle(X),options,V);
        myFolder = './BenchMark/Extractions/baseline/';
        saveEstimates(x_Source_b,options,myFolder)
    elseif(cond2)
        [x_Source_b,Mask_b] = temporal_Reconstruction(W_b,H_b,angle(X),options,V);
        myFolder = './BenchMark/Extractions/baseline/';
        saveEstimates(x_Source_b,options,myFolder)
    elseif(cond3)
        [x_Source_s,Mask_s] = temporal_Reconstruction(W_s,H_s,angle(X),options,V);
        myFolder = './BenchMark/Extractions/sparse/';
        saveEstimates(x_Source_s,options,myFolder)
    end
   
    figure7=figure;
    for i=1:options.K
        subplot(options.K,1,i);
        clims = [0 1];
        imagesc(Mask(:,time_frame_min:time_frame_max,i),clims)
        set(gca,'YDir','normal') %Flip ordinate
        set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
        set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
        colormap bone
    end
    subplot(options.K,1,1)
    if(options.beta==1)
        title('min-vol KL-NMF','FontSize',16, 'Interpreter','latex')
    elseif(options.beta==0)
        title('min-vol IS-NMF','FontSize',16, 'Interpreter','latex')
    end
    saveas(figure7,'./BenchMark/Masks_minvol.fig')
    saveas(figure7,'./BenchMark/Masks_minvol.png')
    if(cond1 || cond3)
        figure8=figure;
        for i=1:options.K
            subplot(options.K,1,i);
            clims = [0 1];
            imagesc(Mask_s(:,time_frame_min:time_frame_max,i),clims)
            set(gca,'YDir','normal') %Flip ordinate
            set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
            set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
            colormap bone
        end
        subplot(options.K,1,1)
        if(options.beta==1)
            title('sparse KL-NMF','FontSize',16, 'Interpreter','latex')
        elseif(options.beta==0)
            title('sparse IS-NMF','FontSize',16, 'Interpreter','latex')
        end
        saveas(figure8,'./BenchMark/Masks_sparse.fig')
        saveas(figure8,'./BenchMark/Masks_sparse.png')
    end
    if(cond1 || cond2)
        figure9=figure;
        for i=1:options.K
            subplot(options.K,1,i);
            clims = [0 1];
            imagesc(Mask_b(:,time_frame_min:time_frame_max,i),clims)
            set(gca,'YDir','normal') %Flip ordinate
            set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
            set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
            colormap bone
        end
        subplot(options.K,1,1)
        if(options.beta==1)
            title('baseline KL-NMF','FontSize',16, 'Interpreter','latex')
        elseif(options.beta==0)
            title('baseline IS-NMF','FontSize',16, 'Interpreter','latex')
        end
        saveas(figure9,'./BenchMark/Masks_baseline.fig')
        saveas(figure9,'./BenchMark/Masks_baseline.png')
    end
    
    
    disp('End of BENCHMARK')
end
%% BLOC 5: BSS Evaluation metrics
if(options.metric_eval==1)
    % % Metrics computation 
    myFolder = './Extractions';
    extension='*.wav';
    [SDR,SIR,SAR,perm] = metric_Eval(myFolder,extension,options);
    table(perm,SDR,SIR,SAR)
end

%% SOURCES:
% Paper: Valentin Leplat, Nicolas Gillis, Man Shun Ang. "Blind Audio Source Separation with Minimum-Volume Beta-Divergence NMF"          
% BSS eval tool box: http://bass-db.gforge.inria.fr/bss_eval/
% Sparse NMF code:        https://www.merl.com/publications/docs/TR2015-023.pdf
%                         such copying is by permission of Mitsubishi ElectricResearch Laboratories, Inc.; an acknowledgment of the authors and individual contributions to the work; and allapplicable portions of the copyright notic
% Munkres (Hungarian) Algorithm for Linear Assignment Problem: http://csclab.murraystate.edu/~bob.pilgrim/445/munkres.html
% Music sample drum+bass: http://isse.sourceforge.net/demos.html
% Music Sample Prelude  from J.S Bach: https://www.youtube.com/watch?v=ZlbK5r5mBH4