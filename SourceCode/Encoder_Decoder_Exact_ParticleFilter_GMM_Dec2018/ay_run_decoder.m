function ay_run_decoder(encoder_file_name,which_method,res_file_to_save,base_file_to_save)
% define number of particles
PARTCILE_NO = 4000;

%% Decode Setting
% complete resolution
dXY   = 8;
dVXY  = 8;
%% Load Data For Intensity Calculation 
% load Intensity Calculation Components - this bring path as well
load(encoder_file_name);
% Path, Cell, Kernel
Mark.Path   = CellPath;
Mark.Cell   = PCell;
Mark.Kernel = Kernel;
%Mark.dT     = (CellPath.T(2)-CellPath.T(1));
Mark.dxy    = dXY;
Mark.dvxy   = dVXY;

%% Define State Transition Process
% here, we focus on 2-D model. 
Mark.dT = TestPath.T(2)-TestPath.T(1);
TransProc.A = [1    0;
               0    1];
TransProc.S = [6.0      0;
                0      6.0];
            
%% GMM or Single Gaussian Model (2D)
if which_method ==3 || which_method ==4 
    choose_criteria   = 2;
    %% maximum number of mixtures
    N_MIX = 15;
    IterA = 250;
    IterB = 50;
    hx  = 1;    hy  = 1;
    steps = [hx;hy];
    ind_a    = 1;
    INIT_VAR = 3;
    M        = [TestPath.X(ind_a) TestPath.Y(ind_a)];
    Q        = [INIT_VAR    0        ;
                0           INIT_VAR];
    % model A is based on mixture addition
    PostMixedModelA.n_mix       = 1;
    PostMixedModelA.Model{1}.M  = M;
    PostMixedModelA.Model{1}.S  = Q;       
    PostMixedModelA.Model{1}.W  = 1;
    % % ------- Complete Model
    EXTND = 10;
    mx = min(TestPath.X)-EXTND:2:max(TestPath.X)+EXTND;
    my = min(TestPath.Y)-EXTND:2:max(TestPath.Y)+EXTND;

    PV  = zeros(length(mx)*length(my),2);
    ind = 1;
    for i=1:length(mx)
        for j=1:length(my)
            PV(ind,:)=[mx(i) my(j)];
            ind = ind + 1;
        end
    end
    dE  = [mx(2)-mx(1);my(2)-my(1)];
    HPD = 0.95;
end

%% GMM or Single Gaussian Model (4D)
if which_method == 5 
    ind_a   =1;
    %% whcih method
    choose_criteria   = 2;
    %% maximum number of mixture
    N_MIX = 15;
    IterA = 250;
    IterB = 50;

    hx  = 1;    hy  = 1;
    hvx = 1;    hvy = 1;
    steps = [hx;hy;hvx;hvy];

    %% Define State Transition Process
    % here, we focus on 2-D model. 
    Mark.dT     = TestPath.T(2)-TestPath.T(1);
    TransProc.A = [1    0    Mark.dT     0;
                   0    1    0           Mark.dT;
                   0    0    1           0;
                   0    0    0           1];
    TransProc.S = [ 6.0     0    0      0;
                    0      6.0   0      0;
                    0       0    12     0;
                    0       0    0      12];

    %% Build initial mixture model - posterior model
    % ---- Model A
    % model A is based on shrinkage model
    INIT_VAR   = 3;
    INIT_V_VAR = 6;
    M = [TestPath.X(ind_a) TestPath.Y(ind_a)    TestPath.VX(ind_a) TestPath.VY(ind_a)];
    Q = [INIT_VAR    0          0           0;
         0           INIT_VAR   0           0;
         0           0         INIT_V_VAR   0;
         0           0         0    INIT_V_VAR];

    % ----- Model A
    % model A is based on mixture addition
    PostMixedModelA.n_mix       = 1;
    PostMixedModelA.Model{1}.M  = M;
    PostMixedModelA.Model{1}.S  = Q;       
    PostMixedModelA.Model{1}.W  = 1;

    % % ------- Complete Model
    EXTND = 10;
    mx = min(TestPath.X)-EXTND:2:max(TestPath.X)+EXTND;
    my = min(TestPath.Y)-EXTND:2:max(TestPath.Y)+EXTND;
    vx = min(TestPath.VX)-EXTND:2:max(TestPath.VX)+EXTND;
    vy = min(TestPath.VY)-EXTND:2:max(TestPath.VY)+EXTND;

    PV  = zeros(length(mx)*length(my)*length(vx)*length(vy),4);
    PXY = zeros(length(mx)*length(my),2);
    ind = 1;
    ind_xy = 1;
    for i=1:length(mx)
        for j=1:length(my)
            PXY(ind_xy,:)=[mx(i) my(j)];
            ind_xy =  ind_xy + 1;
            for vi=1:length(vx)
                for vj=1:length(vy)
                    PV(ind,:)=[mx(i) my(j) vx(vi) vy(vj)];
                    ind = ind + 1;
                end
            end
        end
    end
    dE  = [mx(2)-mx(1);my(2)-my(1);mx(2)-mx(1);my(2)-my(1)];
    HPD = 0.95;

end

%% Exact Model
if which_method ==1
    ind_a = 1;
    EXTND = 10;
    mx    = min(TestPath.X)-EXTND:2:max(TestPath.X)+EXTND;
    my    = min(TestPath.Y)-EXTND:2:max(TestPath.Y)+EXTND;

    PV  = zeros(length(mx)*length(my),2);
    ind = 1;
    for i=1:length(mx)
        for j=1:length(my)
            PV(ind,:)=[mx(i) my(j)];
            ind = ind + 1;
        end
    end
    dE  = [mx(2)-mx(1);my(2)-my(1)];
    HPD = 0.95;

    PostComplete = zeros(length(mx),length(my));
    [~,ind_x]=min(abs(TestPath.X(ind_a)-mx));
    [~,ind_y]=min(abs(TestPath.Y(ind_a)-my));
    PostComplete(ind_x,ind_y)=1;
    h  = fspecial('gaussian',10,1); 
    PostComplete = imfilter(PostComplete,h,'replicate');
    
    %% build L for non-spike time
    Mo =[0*MS(ind_a,2:end-4) 0];
    num_time_step = 1;
    Mark.dT = TestPath.T(2)-TestPath.T(1);
    L00  = ay_likelihood(Mo,Mark,num_time_step,mx,my);
end
            

%% Particle Filter Solution, build initial set of particles
if which_method ==2
    ind_a      = 1;
    INIT_VAR   = 3;
    M   =    [TestPath.X(ind_a) TestPath.Y(ind_a)];
    Q   = [INIT_VAR    0        ;
            0           INIT_VAR];
    Ps  = mvnrnd(M,Q,PARTCILE_NO);
    
    %% Graph Section, keep mixture model
    file_name = [base_file_to_save '_particle2_' num2str(ind_a)];
    save(file_name,'Ps','M');

end

%% Main Loop which Runs Different Methods
% Run the decoder-encoder loop
ind_b = length(TestPath.X);
for t= 2:1:ind_b
    t
    %% Exact Solution (2D)
    if which_method == 1
        %% start timer
        tic
        %% Decoding step
        Mo      = [MS(t,2:end-4) 0];
        Mark.dT = TestPath.T(t)-TestPath.T(t-1);
        %% this is memory intensive method, where we build L00 ahead of time
        %if sum(Mo)== 0
        %    L = L00;
        %else
        %    % calculate likelihood
        %    L  = ay_likelihood(Mo,Mark,num_time_step,mx,my);
        %end
        %% build the likelihood on each time point
        L  = ay_likelihood(Mo,Mark,num_time_step,mx,my);
        % run filter
        num_time_step = 1;
        PostComplete = ay_filter(L,PostComplete,mx,my,TransProc,num_time_step);
        %% stop timer, calculate processing time
        proc_time = toc; 
        
        %% Performance Section, save result
        xPostComplete = zeros(length(mx)*length(my),1);
        ind  = 1;
        for xx=1:length(mx)
             for yy=1:length(my)
                 xPostComplete(ind)=PostComplete(xx,yy);
                ind = ind +1;
             end
        end
        Real_XY = [MS(t,end-3) MS(t,end-2)];
        [inHPD_Complete,SErr_Complete,MSErr_Complete,AreaInd,MX] = ay_performance_estimate_2d(1,Real_XY,xPostComplete,PV,dE,HPD);
        Result =[t sum(Mo) MS(t,1) MS(t,end-3) MS(t,end-2)   inHPD_Complete   SErr_Complete   MSErr_Complete   AreaInd MX' proc_time];
        fid  = fopen([res_file_to_save '_exact2.txt'],'at');
        for z = 1:length(Result)
             fprintf(fid,'%f  ',Result(z));
        end
        fprintf(fid,'\r\n');
        fclose(fid);

        %% Save mixture model per each time step
        file_name = [base_file_to_save '_exact2_' num2str(t)];
        save(file_name,'PostComplete','mx','my');
    end
    %% Particle Filter Solution
    if which_method == 2   
        %% Start timer
        tic
        %% Decoding Step
        Mo            = [MS(t,2:end-4) 0];
        Mark.dT       = TestPath.T(t)-TestPath.T(t-1);
        %% we need to draw n sampler per each time point
        Ps  = ay_gen_particle_samples(Ps,TransProc,Mo,Mark);
        %% Caculate processing time, stop timer
        proc_time = toc;
        
        %% Perfromance Result
        % mean of the estimate
        Mx = sum(Ps)/size(Ps,1);
        % keep mean here
        Err = sqrt((Mx(1)-MS(t,end-3)).^2 + (Mx(2)-MS(t,end-2)).^2); 
        Result =[t sum(Mo) MS(t,1 ) MS(t,end-3)   MS(t,end-2) Err  Mx proc_time];
        fid = fopen([res_file_to_save '_partcile2.txt'],'at');
        for z=1:length(Result)
             fprintf(fid,'%f  ',Result(z));
        end
        fprintf(fid,'\r\n');
        fclose(fid);

        %% Graph Section, keep mixture model
        file_name = [base_file_to_save '_particle_' num2str(t)];
        save(file_name,'Ps','Mx');
    end
    %% Single Gaussian
    if which_method == 3
        tic
        %% Decoding Step
        num_time_step = 1;
        Mo = [MS(t,2:end-4) 0];
        Mark.dT = TestPath.T(t)-TestPath.T(t-1);
    
        PreMixedModelA = ay_prior(PostMixedModelA,TransProc,num_time_step);
        if sum(Mo)==0   % no spike time
            % a more accurate solution is Gaussian Approximation Here
            % mode=0, is just copying Prior to Posterior
            % mode=1, is using the idea proposed in paper
            [tempPostMixedModel,update] = ay_gmm_update(PreMixedModelA,Mark,steps);
            if update
                Ps  = ay_gen_particle(PARTCILE_NO,Mo,PreMixedModelA,Mark);
                PostMixedModelA= ay_gmm_select(1,1,IterA,IterB,Ps,choose_criteria);
            else
               PostMixedModelA = tempPostMixedModel;
            end
        end
        if sum(Mo)>0     % on spike time
            % we draw samples, and then try to re-estimate Guassian mixture
            Ps = ay_gen_particle(PARTCILE_NO,Mo,PreMixedModelA,Mark);
            % run gaussian mixture process
            PostMixedModelA= ay_gmm_select(1,1,IterA,IterB,Ps,choose_criteria);
        end
        proc_time = toc;
        %% Performance Section, save result
        Real_XY = [MS(t,end-3) MS(t,end-2)];
        [inHPD_Mixed,SErr_Mixed,MSErr_Mixed,AreaInd,MX]    = ay_performance_estimate_2d(2,Real_XY,PostMixedModelA,PV,dE,HPD);
        Result =[t sum(Mo) MS(t,1) MS(t,end-3)   MS(t,end-2)   inHPD_Mixed  SErr_Mixed  MSErr_Mixed  PostMixedModelA.n_mix  AreaInd MX' proc_time];
        fid = fopen([res_file_to_save '_gauss2.txt'],'at');
        for z=1:length(Result)
             fprintf(fid,'%f  ',Result(z));
        end
        fprintf(fid,'\r\n');
        fclose(fid);

        %% Graph Section, keep mixture model
        file_name = [base_file_to_save '_gauss2_' num2str(t)];
        save(file_name,'PostMixedModelA');
    end
    %% GMM (2D)
    if which_method == 4
        
        tic 
        %% Decoding loop
        num_time_step = 1;
        Mo = [MS(t,2:end-4) 0];
        Mark.dT = TestPath.T(t)-TestPath.T(t-1);
        PreMixedModelA = ay_prior(PostMixedModelA,TransProc,num_time_step);
        if sum(Mo)==0   % no spike time
            % a more accurate solution is Gaussian Approximation Here
            % mode=0, is just copying Prior to Posterior
            % mode=1, is using the idea proposed in paper
            [tempPostMixedModel,update] = ay_gmm_update(PreMixedModelA,Mark,steps);
            if update
                Ps  = ay_gen_particle(PARTCILE_NO,Mo,PreMixedModelA,Mark);
                PostMixedModelA= ay_gmm_select(1,N_MIX,IterA,IterB,Ps,choose_criteria);
            else
               PostMixedModelA = tempPostMixedModel;
            end
        end
        if sum(Mo)>0     % on spike time
            % we draw samples, and then try to re-estimate Guassian mixture
            Ps = ay_gen_particle(PARTCILE_NO,Mo,PreMixedModelA,Mark);
            % run gaussian mixture process
            PostMixedModelA= ay_gmm_select(1,N_MIX,IterA,IterB,Ps,choose_criteria);
        end
        proc_time = toc;
        %% Performance Section, save result
        Real_XY = [MS(t,end-3) MS(t,end-2)];
        [inHPD_Mixed,SErr_Mixed,MSErr_Mixed,AreaInd,MX]    = ay_performance_estimate_2d(2,Real_XY,PostMixedModelA,PV,dE,HPD);
        Result =[t sum(Mo) MS(t,1 ) MS(t,end-3)   MS(t,end-2)   inHPD_Mixed  SErr_Mixed  MSErr_Mixed  PostMixedModelA.n_mix  AreaInd MX' proc_time];
        fid = fopen([res_file_to_save '_gmm2.txt'],'at');
        for z=1:length(Result)
             fprintf(fid,'%f  ',Result(z));
        end
        fprintf(fid,'\r\n');
        fclose(fid);

        %% Graph Section, keep mixture model
        file_name = [base_file_to_save '_gmm2_' num2str(t)];
        save(file_name,'PostMixedModelA');
   
    end
    if which_method == 5
        tic
        %% Decodign Loop
        num_time_step      = 1;
        Mo = [MS(t,2:end-4) 0];
        Mark.dT = TestPath.T(t)-TestPath.T(t-1);

        PreMixedModelA = ay_prior(PostMixedModelA,TransProc,num_time_step);
        if sum(Mo)==0   % no spike time
            % a more accurate solution is Gaussian Approximation Here
            % mode=0, is just copying Prior to Posterior
            % mode=1, is using the idea proposed in paper
            [tempPostMixedModelA,update] = ay_gmm_update(PreMixedModelA,Mark,steps);
            if update
               Ps  = ay_gen_particle(PARTCILE_NO,Mo,PreMixedModelA,Mark);
               PostMixedModelA= ay_gmm_select(1,N_MIX,IterA,IterB,Ps,choose_criteria);
            else
               PostMixedModelA = tempPostMixedModelA;
            end
        end
        if sum(Mo)>0     % on spike time
            % we draw samples, and then try to re-estimate Guassian mixture
            Ps = ay_gen_particle(PARTCILE_NO,Mo,PreMixedModelA,Mark);
            % run gaussian mixture process
            PostMixedModelA= ay_gmm_select(1,N_MIX,IterA,IterB,Ps,choose_criteria);
        end
        proc_time = toc;

        %% Performance Section, save result
        Real_XY_VXVY = [MS(t,end-3) MS(t,end-2) MS(t,end-1) MS(t,end)];
        [inHPD_Mixed,SErr_Mixed,MSErr_Mixed,PErr,AreaInd,MX] = ay_performance_estimate_xx(Real_XY_VXVY,PostMixedModelA,PV,dE,HPD,PXY);
        Result =[t sum(Mo) MS(t,end-3 ) MS(t,end-2)   MS(t,end-1)   MS(t,end) inHPD_Mixed   SErr_Mixed     MSErr_Mixed     PostMixedModelA.n_mix   PErr  AreaInd MX' proc_time];
        fid = fopen([res_file_to_save '_gmm4.txt'],'at');
        for z=1:length(Result)
             fprintf(fid,'%f  ',Result(z));
        end
        fprintf(fid,'\r\n');
        fclose(fid);

        %% Graph Section, keep mixture model
        file_name = [base_file_to_save '_gmm4_' num2str(t)];
        save(file_name,'PostMixedModelA');
    end
    
    
%% To generate pictures, movies, etc
%     subplot(1,4,2)
%     Lt=ay_gmm_plot(PostMixedModelB,mx,my);
%     %imagesc(mx,my,sqrt(Lt'/sum(sum(Lt))));
%     imagesc(mx,my,sqrt(Lt'));
%     hold on
%     plot(MS(:,end-1),MS(:,end),'.','LineWidth',0.1,'Color',[0.9 0.9 0.9]);
%     imagesc(mx,my,sqrt(Lt'/sum(sum(Lt))));
%     alpha 0.8
%     plot(MS(t,end-1),MS(t,end),'ro','Markersize',8,'LineWidth',2);
%     for zz=1:PostMixedModelB.n_mix
%           plot(PostMixedModelB.Model{zz}.M(1),PostMixedModelB.Model{zz}.M(2),'k+','Markersize',8,'LineWidth',1);
%     end
%     hold off
%     xlabel('X')
%     ylabel('Y')
% 
%     title(['Single'])
%     axis square
%     set(gca,'fontsize',12);
%     colormap parula
%     
%     subplot(1,4,3)
%     Lt=ay_gmm_plot(PostMixedModelA,mx,my);
%     %imagesc(mx,my,sqrt(Lt'/sum(sum(Lt))));
%     imagesc(mx,my,sqrt(Lt'));
%     hold on
%     plot(MS(:,end-1),MS(:,end),'.','LineWidth',0.1,'Color',[0.9 0.9 0.9]);
%     imagesc(mx,my,sqrt(Lt'/sum(sum(Lt))));
%     alpha 0.8
%     plot(MS(t,end-1),MS(t,end),'ro','Markersize',8,'LineWidth',2);
%     for zz=1:PostMixedModelA.n_mix
%           plot(PostMixedModelA.Model{zz}.M(1),PostMixedModelA.Model{zz}.M(2),'k+','Markersize',8,'LineWidth',1);
%     end
%     hold off
%     xlabel('X')
%     ylabel('Y')
% 
%     title(['Hybrid #' num2str(PostMixedModelA.n_mix)])
%     axis square
%     set(gca,'fontsize',12);
%     colormap parula
%     
%     subplot(1,4,4)
%     %imagesc(mx,my,sqrt(Lt'/sum(sum(Lt))));
%     imagesc(mx,my,sqrt(PostComplete'));
%     hold on
%     plot(MS(:,end-1),MS(:,end),'.','LineWidth',0.1,'Color',[0.9 0.9 0.9]);
%     alpha 0.8
%     plot(MS(t,end-1),MS(t,end),'ro','Markersize',8,'LineWidth',2);
%     imagesc(mx,my,sqrt(PostComplete'/sum(sum(PostComplete))));
%     alpha 0.8
%     hold off
%     xlabel('X')
%     ylabel('Y')
% 
%     title(['Exact'])
%     axis square
%     set(gca,'fontsize',12);
%     colormap parula
%     
% 
%     subplot(1,4,1)
%     ind = find(Mo);
%     plot(Mark.Path.X,Mark.Path.Y,'.','MarkerSize',3,'COLOR',[0.8 0.8 0.8]);
%     hold on
%     text_title = [];
%     if ~isempty(ind)
%         text_title= 'A Neuron Fired ';
%         for ll=1:length(ind)
%             plot(Mark.Cell{ind(ll)}.X,Mark.Cell{ind(ll)}.Y,'*','MarkerSize',3);
%             text_title = strcat(text_title,[' ' num2str(ind(ll))]);
%         end
%     end    
%     xlabel('X')
%     ylabel('Y')
%     title_str =['Time:' num2str(t)];title_str=strcat(title_str,text_title);
%     title(title_str);
%     set(gca,'YDir','Reverse');
%     axis square
%     set(gca,'fontsize',12);
%     hold off
% 
%     set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 8])
% 
%     print('-dpng',['images\result_' num2str(fig_ind) '.png']);
    

end
