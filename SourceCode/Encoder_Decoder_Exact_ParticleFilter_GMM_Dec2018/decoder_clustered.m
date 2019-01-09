function decoder_clustered(MIfile,startind,endind)

% Model GMM & Partcile Method
% Updated March 22 AY
% adapted from ay_run_decoder_may_10

% This loop runs between startind and endind

ind_a = startind;
ind_b = endind;

%% Decode Setting
% particle number
PARTCILE_NO = 8000;
% for gradient calculation
hx  = 4*0.25;
hy  = 4*0.25;
% likelihood method
likelihood_method = 1;
% choose criteria
choose_criteria   = 2;
% number of mixture and number of iterations
N_MIX   = 15;
N_Iter  = 500; 
% complete resolution
dXY = 2;

%% Load Data For Intensity Calculation 
% load Intensity Calculation Components - this bring path as well
load(MIfile);
% Path, Cell, Kernel
Mark.Path   = CellPath;
Mark.Cell   = Cell;
Mark.Kernel = Kernel;
Mark.dT     = (CellPath.T(2)-CellPath.T(1));
Mark.dxy    = dXY;

%% Define State Transition Process
% here, we focus on 4-D model. 
TransProc.A = [1   0;
               0   1];
% I have checked the peack of dx, dy 
% The peak is devided by 1.5 and it defines std 
% diagonal variance is std^2
TransProc.S = [2.5 .25;                  %[5        0.5;
                .25 4.5];                    %0.5       9];

%% Build initial mixture model - posterior model
% % ---- Model A
% % model A is based on shrinkage model
 INIT_VAR  = 5;
% PostMixedModelA.n_mix       = 1;
% PostMixedModelA.Model{1}.M  = [TestPath.X(ind_a) TestPath.Y(ind_a)];
% PostMixedModelA.Model{1}.S  = INIT_VAR *eye(2,2);
% PostMixedModelA.Model{1}.W  = 1;
% % prior model
% PreMixedModelA = PostMixedModelA;
% 
% ----- Model B
% model B is based on mixture addition
PostMixedModelB.n_mix       = 1;
PostMixedModelB.Model{1}.M  = [TestPath.X(ind_a) TestPath.Y(ind_a)];
PostMixedModelB.Model{1}.S  = INIT_VAR *eye(2,2);
PostMixedModelB.Model{1}.W  = 1;
% prior model
PreMixedModelB = PostMixedModelB;
% 
% % ------- Single Mixture Model
% % Build initial mixture model - posterior model
% PostSingleModel.n_mix       = 1;
% PostSingleModel.Model{1}.M  = [TestPath.X(ind_a) TestPath.Y(ind_a)];
% PostSingleModel.Model{1}.S  = INIT_VAR *eye(2,2);
% PostSingleModel.Model{1}.W  = 1;
% % prior model
% PreSingledModel = PostSingleModel;

% ------- Complete Model
% Complete Model
[L,mx,my]    = ay_likelihood(zeros(length(Mark.Cell),1),Mark);
PostComplete = ones(length(mx),length(my));
for i=1:length(mx)
    for j=1:length(my)
         PostComplete(i,j)= mvnpdf([mx(i) my(j)],[TestPath.X(ind_a) TestPath.Y(ind_a)],INIT_VAR*eye(2,2));
    end
end

PV  = zeros(length(mx)*length(my),2);
ind = 1;
for i=1:length(mx)
    for j=1:length(my)
        PV(ind,:)=[mx(i) my(j)];
        ind = ind + 1;
    end
end
dE = [mx(2)-mx(1);my(2)-my(1)];



%% Main Loop which Runs Different Method
% Run the decoder-encoder loop
fig_ind       = 10000;
for t=ind_a:1:ind_b
    t
    %% Current Data
    num_time_step = 1;
    Mo = [MS(t,2:end-2) 0];
    %% check this in ay_run
    %Mo = [MS(t,2:end-2)];
    
    %% Complete Method
%     tic
    % calculate likelihood
    [L,mx,my]    = ay_likelihood(Mo,Mark,num_time_step);
    % run filter
    PostComplete = ay_filter(L,PostComplete,mx,my,TransProc,num_time_step);
    exact_solution_time = toc;
    
%     %% Single Model
%     tic
%     % run one step prediction
%     PreSingleModel = ay_prior(PostSingleModel,TransProc,num_time_step);
%     if sum(Mo)==0   % no spike time
%        % a more accurate solution is Gaussian Approximation Here
%        % mode=0, is just copying Prior to Posterior
%        % mode=1, is using the idea proposed in paper
%        PostSingleModel = ay_gmm_update(PreSingleModel,Mark,hx,hy,Mo,num_time_step);
%     else            % on spike time
%        %tempPreSingleModel=ay_gmm_update(1,likelihood_method,PreSingleModel,Mark,hx,hy,Mo,num_time_step); 
%        %[Ps,Pw] = ay_gen_particle(1,likelihood_method,PARTCILE_NO,Mo,tempPreSingleModel,Mark);
%        % we draw samples, and then try to re-estimate Guassian mixture
%        [Ps,Pw] = ay_gen_particle(0.025*PARTCILE_NO,Mo,PreSingleModel,Mark);
%        % run gaussian mixture process
%        PostSingleModel= ay_gmm_select(1,1,N_Iter,Ps,Pw,choose_criteria);
%        % for display
%        Ps_single = Ps;
%        Pw_single = Pw;
%     end
%     single_gmm_time = toc;
    
%     %% Mixed Model A - Addition
%     tic
%     % run one step prediction
%     PreMixedModelA = ay_prior(PostMixedModelA,TransProc,num_time_step);
%     if sum(Mo)==0   % no spike time
%         % a more accurate solution is Gaussian Approximation Here
%         % mode=0, is just copying Prior to Posterior
%         % mode=1, is using the idea proposed in paper
%         PostMixedModelA = ay_gmm_update(PreMixedModelA,Mark,hx,hy,Mo,num_time_step);
%     else            % on spike time
%         %tempPreMixedModelA=ay_gmm_update(1,likelihood_method,PreMixedModelA,Mark,hx,hy,Mo,num_time_step); 
%         %[Ps,Pw] = ay_gen_particle(1,likelihood_method,PARTCILE_NO,Mo,tempPreMixedModelA,Mark);
%         % we draw samples, and then try to re-estimate Guassian mixture
%         [Ps,Pw] = ay_gen_particle(0.025*PARTCILE_NO,Mo,PreMixedModelA,Mark);
%         % run gaussian mixture process
%         PostMixedModelA= ay_gmm_select(2,N_MIX,N_Iter,Ps,Pw,choose_criteria);
%         % for display
%         Ps_A = Ps;
%         Pw_A = Pw;
%     end
%     multiple_gmm_a_time = toc;
    
    %% Mixed Model B - Shrinkage
    tic
    % run one step prediction
    PreMixedModelB = ay_prior(PostMixedModelB,TransProc,num_time_step);
    if sum(Mo)==0   % no spike time
        % a more accurate solution is Gaussian Approximation Here
        % mode=0, is just copying Prior to Posterior
        % mode=1, is using the idea proposed in pape
        PostMixedModelB = ay_gmm_update(PreMixedModelB,Mark,hx,hy,Mo,num_time_step);
    else            % on spike time
        %tempPreMixedModelB=ay_gmm_update(1,likelihood_method,PreMixedModelB,Mark,hx,hy,Mo,num_time_step); 
        %[Ps,Pw] = ay_gen_particle(1,likelihood_method,PARTCILE_NO,Mo,tempPreMixedModelB,Mark);
        % we draw samples, and then try to re-estimate Guassian mixture
        [Ps,Pw] = ay_gen_particle(PARTCILE_NO,Mo,PreMixedModelB,Mark);
        % run gaussian mixture process
        PostMixedModelB= ay_gmm_select(1,N_MIX,N_Iter,Ps,Pw,choose_criteria);
        % for display
        Ps_B = Ps;
        Pw_B = Pw;
    end
    multiple_gmm_b_time = toc;
     
%     %% Performance Measure
%      xPostComplete = zeros(length(mx)*length(my),2);
%      ind = 1;
%      for xx=1:length(mx)
%          for yy=1:length(my)
%              xPostComplete(ind)=PostComplete(xx,yy);
%              ind = ind +1;
%          end
%      end
%     [inHPD_Complete,SErr_Complete,MSErr_Complete] = ay_performance_estimate_xx(1,[MS(t,end-1) MS(t,end)],xPostComplete,PV,dE,0.95);
%     
%     [inHPD_Mixed_A,SErr_Mixed_A,MSErr_Mixed_A]    = ay_performance_estimate_xx(2,[MS(t,end-1) MS(t,end)],PostMixedModelA,PV,dE,0.95);
%     [inHPD_Mixed_B,SErr_Mixed_B,MSErr_Mixed_B]    = ay_performance_estimate_xx(2,[MS(t,end-1) MS(t,end)],PostMixedModelB,PV,dE,0.95);
%     [inHPD_Single,SErr_Single,MSErr_Single]       = ay_performance_estimate_xx(2,[MS(t,end-1) MS(t,end)],PostSingleModel,PV,dE,0.95);
%     Result =[t MS(t,1 ) MS(t,end-1)   MS(t,end-2)   inHPD_Complete   inHPD_Single   inHPD_Mixed_A  inHPD_Mixed_B   SErr_Complete   SErr_Single   SErr_Mixed_A  SErr_Mixed_B   MSErr_Complete  MSErr_Single   MSErr_Mixed_A  MSErr_Mixed_B exact_solution_time single_gmm_time multiple_gmm_a_time multiple_gmm_b_time PostMixedModelA.n_mix PostMixedModelB.n_mix];
%     %Result =[t MS(t,1 ) MS(t,2)   MS(t,3)   MS(t,4)    inHPD_Single   inHPD_Mixed     SErr_Single   SErr_Mixed      MSErr_Single   MSErr_Mixed];
%     fid = fopen('Result.txt','at');
%     for z=1:length(Result)
%              fprintf(fid,'%f  ',Result(z));
%     end
%     fprintf(fid,'\r\n');
%     fclose(fid);

    %% Graph Result
    %figure(1)
    figure; set(gcf,'Position',[675 148 918 834]);
    % Path, Current Point, Spike Information
    %subplot(2,5,1)
    subplot(2,2,1); axis square
    plot(MS(:,end-1),MS(:,end),'.','LineWidth',2,'Color',[0.7 0.7 0.7]);
    hold on;
    plot(MS(t,end-1),MS(t,end),'ro','Markersize',8,'LineWidth',2);
    if sum(Mo)>=1
        plot(MS(t,end-1),MS(t,end),'g*','Markersize',4,'LineWidth',2);
    end
    hold off;
    xlabel('X')
    ylabel('Y')
    title_str =['Time:' num2str(t)];
    spk_ind   = find(MS(t,2:end-2));
    if ~isempty(spk_ind)
        title_str = strcat(title_str,'(');
        for zz=1:length(spk_ind)-1
            title_str = strcat(title_str,num2str(spk_ind(zz)));
            title_str = strcat(title_str,',');
        end
        title_str = strcat(title_str,num2str(spk_ind(end)));
        title_str = strcat(title_str,')');
    end
    title(title_str);
    set(gca,'YDir','Reverse');
    axis([min(Mark.Path.X-15)  max(Mark.Path.X+15)   min(Mark.Path.Y-15)  max(Mark.Path.Y+15)])
    axis square
    set(gca,'fontsize',12);

%     % Partciles Model Single
%     % Ps = Ps_single;  AKG
%     Ps = PostSingleModel.Model{1}.S;
%     
%     subplot(2,5,3)
%     imagesc(mx,my,log(1e-6+(L'/sum(sum(L)))));
%     axis square
%     set(gca,'fontsize',12);
%     colormap parula
%     alpha 0.8
%     hold on;
%     plot(Ps(:,1),Ps(:,2),'m.');
%     set(gca,'YDir','Reverse');
%     hold off;
%     xlabel('X')
%     title('Gaussian Particles')
        
%     % Partciles Model A
%     %Ps = Ps_A;   AKG
%     Ps = PostMixedModelA.Model{1}.S;
%     
%     subplot(2,5,4)
%     imagesc(mx,my,log(1e-6+(L'/sum(sum(L)))));
%     axis square
%     set(gca,'fontsize',12);
%     colormap parula
%     alpha 0.8
%     hold on;
%     plot(Ps(:,1),Ps(:,2),'m.');
%     set(gca,'YDir','Reverse');
%     hold off;
%     xlabel('X')
%     title('Greedy Particles')
    
    
%     % Partciles Model B
%     %Ps = Ps_B;   AKG
%     Ps = PostMixedModelB.Model{1}.S;
%     
%     subplot(2,5,5)
%     imagesc(mx,my,log(1e-6+(L'/sum(sum(L)))));
%     axis square
%     set(gca,'fontsize',12);
%     colormap parula
%     alpha 0.8
%     hold on;
%     plot(Ps(:,1),Ps(:,2),'m.');
%     set(gca,'YDir','Reverse');
%     hold off;
%     xlabel('X')
%     title('Hybrid Particles')
    
    %% Show Likelihood function
    %subplot(2,5,6)
    subplot(2,2,2); axis square
    imagesc(mx,my,log(1e-6+(L'/sum(sum(L)))));
    xlabel('X')
    ylabel('Y')
    title('Likelihood')
    axis square
    set(gca,'fontsize',12);
    colormap parula
%         
        
    subplot(2,2,3); axis square
    imagesc(mx,my,sqrt(PostComplete'/sum(sum(PostComplete))));
    hold on
    plot(MS(:,end-1),MS(:,end),'.','LineWidth',0.1,'Color',[0.9 0.9 0.9]);
    hold on
    imagesc(mx,my,sqrt(PostComplete'/sum(sum(PostComplete))));
    alpha 0.8
    plot(MS(t,end-1),MS(t,end),'ro','Markersize',8,'LineWidth',2);
    hold off
    colormap parula
    xlabel('X')
    %ylabel('Y')
    title('Exact Solution')
    axis square
    set(gca,'fontsize',12);


%     subplot(2,5,8)
%     mx = min(Mark.Path.X-15):2.5:max(Mark.Path.X+15);
%     my = min(Mark.Path.Y-15):2.5:max(Mark.Path.Y+15);
%     Lt=ay_gmm_plot(PostSingleModel,mx,my);
%     imagesc(mx,my,sqrt(Lt'/sum(sum(Lt))));
%     hold on
%     plot(MS(:,end-1),MS(:,end),'.','LineWidth',0.1,'Color',[0.9 0.9 0.9]);
%     imagesc(mx,my,sqrt(Lt'/sum(sum(Lt))));
%     alpha 0.8
%     plot(MS(t,end-1),MS(t,end),'ro','Markersize',8,'LineWidth',2);
%     for zz=1:PostSingleModel.n_mix
%         plot(PostSingleModel.Model{zz}.M(1),PostSingleModel.Model{zz}.M(2),'k+','Markersize',8,'LineWidth',1);
%     end
%     hold off
%     xlabel('X')
%     %ylabel('Y')
%     title('Gaussian')
%     colormap parula
%     axis square
%     set(gca,'fontsize',12);

       
%     subplot(2,5,9)
%     Lt=ay_gmm_plot(PostMixedModelA,mx,my);
%     imagesc(mx,my,sqrt(Lt'/sum(sum(Lt))));
%     hold on
%     plot(MS(:,end-1),MS(:,end),'.','LineWidth',0.1,'Color',[0.9 0.9 0.9]);
%     imagesc(mx,my,sqrt(Lt'/sum(sum(Lt))));
%     alpha 0.8
%     plot(MS(t,end-1),MS(t,end),'ro','Markersize',8,'LineWidth',2);
%     for zz=1:PostMixedModelA.n_mix
%         plot(PostMixedModelA.Model{zz}.M(1),PostMixedModelA.Model{zz}.M(2),'k+','Markersize',8,'LineWidth',1);
%     end
%     hold off
%     xlabel('X')
%     %ylabel('Y')
%     title(['Greedy #' num2str(PostMixedModelA.n_mix)])
%     axis square
%     set(gca,'fontsize',12);
%     colormap parula

    %subplot(2,5,10)
    subplot(2,2,4); axis square
    Lt=ay_gmm_plot(PostMixedModelB,mx,my);
    imagesc(mx,my,sqrt(Lt'/sum(sum(Lt))));
    hold on
    plot(MS(:,end-1),MS(:,end),'.','LineWidth',0.1,'Color',[0.9 0.9 0.9]);
    imagesc(mx,my,sqrt(Lt'/sum(sum(Lt))));
    alpha 0.8
    plot(MS(t,end-1),MS(t,end),'ro','Markersize',8,'LineWidth',2);
    for zz=1:PostMixedModelB.n_mix
          plot(PostMixedModelB.Model{zz}.M(1),PostMixedModelB.Model{zz}.M(2),'k+','Markersize',8,'LineWidth',1);
    end
    hold off
    xlabel('X')
    %ylabel('Y')
    title(['Hybrid #' num2str(PostMixedModelB.n_mix)])
    axis square
    set(gca,'fontsize',12);
    colormap parula

    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 60 30])
    savename = sprintf('decode%dto%d/result_%d.png',startind,endind,fig_ind);
    %print('-dpng',['images\result_' num2str(fig_ind) '.png']);
    print('-dpng',savename);
    fig_ind = fig_ind + 1;

    %% Run
    %drawnow
    %pause(0.1) 

   close
end




