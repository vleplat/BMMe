clear all; close all; clc; 
%% ------------------------------------------------------------------------
% Data set loading
%%-------------------------------------------------------------------------

% % Image - "CBCL" data set
load CBCL.mat

% % Image - "FREY" data set
% load freyDataset.mat
% X = M;
% % Image - "ORL" data set
% load imagesORL.mat
% X = M;
% clear Vh Wh M 
% % Image - "Umistim" data set
% load Umistim.mat
% X = M;
% clear M;

% % Document classification - "classic"
% load Datasets/classic.mat
% X = dtm;
% % Document classification - "HiTech"
% load Datasets/hitech.mat
% X = dtm;
% % Document classification - "Sports"
% load Datasets/sports.mat
% X = dtm;
% % Document classification - "La1"
% load Datasets/la1.mat
% X = dtm;

% size computation
[m,n]=size(X);

%% ------------------------------------------------------------------------
% Parameters setup
%%-------------------------------------------------------------------------
darky=[204,204,0]/255;
max_iter=100000; % for ccd only
max_time=100;
options.timemax =max_time;
options.maxiter=inf;

% factorization rank
r=10;
% normalization constant
colX=sum(X,2)/n+eps;
nX=X.*log(X./repmat(colX,1,n)+eps);
nX=sum(nX(:));

%% ------------------------------------------------------------------------
% Benchmark
%%-------------------------------------------------------------------------
numtrial=10;
seed=[];
for idx=1:numtrial
    seed=[seed rng(idx)];
    % initialization
    W= rand(m,r);
    H = rand(r,n);
    % scale initial point 
     WH=W*H;
     alpha=sqrt(sum(X(:))/sum(WH(:)));
     W=alpha*W; 
     H=alpha*H; 

     options.init.W=W; 
     options.init.H=H;
     % run 1 MU step
     options.maxiter=1;
     [W,H]=MU(X,r,options); 
     % initialize by MU
     options.init.W=W; 
     options.init.H=H;

     options.maxiter=inf;
   
     [w0, h0, e_CD, t_CD] = KLnmf(full(X),r,max_iter,max_time,W',H, 1);
     e_CD=e_CD/nX;
     if idx==1 
        CD_obj_save=e_CD;
        CD_time_save=t_CD;
        CD_last_obj=e_CD(end); % save the last objective
     else 
        [~,n1]=size(CD_obj_save);
        l=min(n1,length(e_CD));
        CD_obj_save=[CD_obj_save(:,1:l);e_CD(1:l)];
        CD_time_save=[CD_time_save(:,1:l);t_CD(1:l)];
        CD_last_obj=[CD_last_obj,e_CD(end)];
     end

     [Wmu,Hmu,e_MU,t_MU]=MU_KLNMF(X,r,options);
     e_MU=e_MU/nX;
     fprintf('... MU done, final error = %f \n',e_MU(end));
      if idx==1 
        MU_obj_save=e_MU;
        MU_time_save=t_MU;
        MU_last_obj=e_MU(end);
     else 
        [~,n1]=size(MU_obj_save);
        l=min(n1,length(e_MU));
        MU_obj_save=[MU_obj_save(:,1:l);e_MU(1:l)];
        MU_time_save=[MU_time_save(:,1:l);t_MU(1:l)];
        MU_last_obj=[MU_last_obj,e_MU(end)];
     end
     [W_iMU,H_iMU,e_iMU,t_iMU]=MUe_KLNMF(X,r,options);
     e_iMU=e_iMU/nX; 
     fprintf('... MUe done, final error = %f \n',e_iMU(end));
    
     if idx==1 
        iMU_obj_save=e_iMU;
        iMU_time_save=t_iMU;
        iMU_last_obj=e_iMU(end);
     else 
        [~,n1]=size(iMU_obj_save);
        l=min(n1,length(e_iMU));
        iMU_obj_save=[iMU_obj_save(:,1:l);e_iMU(1:l)];
        iMU_time_save=[iMU_time_save(:,1:l);t_iMU(1:l)];
        iMU_last_obj=[iMU_last_obj,e_iMU(end)];
     end
end

%% ------------------------------------------------------------------------
% Post-processing
%%-------------------------------------------------------------------------
% computation of e_min
e_min=min([min(CD_last_obj),min(MU_last_obj),min(iMU_last_obj)]);


% graph generation - Median Curves
figure;
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultLineLineWidth', 2);
if numtrial>1
    semilogy(median(MU_time_save),median(MU_obj_save)-e_min,'g-.','LineWidth',2);hold on; %MU
    semilogy(median(iMU_time_save),median(iMU_obj_save)-e_min,'r','LineWidth',1.5);hold on; %MUe
    semilogy(median(CD_time_save),median(CD_obj_save)-e_min,'Color',darky,'LineWidth',3);hold on; %MU
else
    semilogy(MU_time_save,MU_obj_save-e_min,'g-.','LineWidth',2);hold on; %MU
    semilogy(iMU_time_save,iMU_obj_save-e_min,'r','LineWidth',1.5);hold on; %MUe
    semilogy(CD_time_save,CD_obj_save-e_min,'Color',darky,'LineWidth',3);hold on; %MU
end
ylabel('D(X,WH)/D(X,0) - e_{min}');
xlabel('Time [sec.]')
legend('MU','MUe','CCD'); 
title([' Median curves over ' num2str(numtrial)  ' runs - $r=$ ' num2str(r)],'FontSize',18, 'Interpreter','latex')
grid on

% Min-curves
figure;
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultLineLineWidth', 2);
semilogy(median(MU_time_save),min(MU_obj_save)-e_min,'g-.','LineWidth',2);hold on; %MU
semilogy(median(iMU_time_save),min(iMU_obj_save)-e_min,'r','LineWidth',1.5);hold on; %MUe
semilogy(median(CD_time_save),min(CD_obj_save)-e_min,'Color',darky,'LineWidth',3);hold on; %MU
ylabel('D(X,WH)/D(X,0) - e_{min}');
xlabel('Time [sec.]')
legend('MU','MUe','CCD'); 
title([' Min curves over ' num2str(numtrial)  ' runs - $r=$ ' num2str(r)],'FontSize',18, 'Interpreter','latex')
grid on


