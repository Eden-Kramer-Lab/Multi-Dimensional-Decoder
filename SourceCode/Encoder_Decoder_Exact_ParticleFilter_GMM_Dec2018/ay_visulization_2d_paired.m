session            = 10;
sub_session_a      = 4;
sub_session_b      = 4;
scale              = 1;

FileName = ['BonWMaze\bonpos' num2str(session) '.mat'];
load(FileName);
PathData     = pos{length(pos)};

%% Mouse Path
P   = PathData{sub_session_a}.data(:,[1 2 3 4 5]);
blocked_ind = find(P(:,2)==0);
correct_ind = find(P(:,2)~=0);
for i=1:length(blocked_ind)
    ind_after  = find(correct_ind >= blocked_ind(i));
    ind_after  = ind_after(1);
    
    ind_before = find(correct_ind <= blocked_ind(i));
    ind_before = ind_before(end);

    if ~isempty(ind_after) && ~isempty(ind_before)
        xp = P(correct_ind(ind_after),2)+ P(correct_ind(ind_before),2);
        yp = P(correct_ind(ind_after),3)+ P(correct_ind(ind_before),3);
        P(blocked_ind(i),2) = xp*0.5;
        P(blocked_ind(i),3) = yp*0.5;
    end
    if isempty(ind_after) && ~isempty(ind_before)
        P(blocked_ind(i),2) = P(correct_ind(ind_before),2);
        P(blocked_ind(i),3) = P(correct_ind(ind_before),3);
    end
    if ~isempty(ind_after) && isempty(ind_before)
        P(blocked_ind(i),2) = P(correct_ind(ind_after),2);
        P(blocked_ind(i),3) = P(correct_ind(ind_after),3);
    end
end
lp = size(P,1);
% Training Set
Pb = P(1:round(lp/2),:);
% Test Set
Pa = P(round(lp/2)+1:end,:);

%% Place Cell Activity
%FileName = ['BonWMaze\emilespikes' num2str(session) '.mat'];
FileName = ['BonWMaze\bonspikes' num2str(session) '.mat'];
load(FileName);
CellData     = spikes{1,session};

%% session A
ms_a   = []; 
ex_ind = 0;
for i=1:length(CellData{1})
    temp = CellData{sub_session_a}{i};
    if ~isempty(temp)
        for j=1:length(temp)
            ex_ind   = ex_ind + 1;
            if ~isempty(temp{j})
                if ~isempty(temp{j}.data)
                    ex_temp  = temp{j}.data(:,1); 
                    ms_a     = [ms_a; ex_temp   ex_ind*ones(size(ex_temp,1),1)]; 
                end
            end
        end
    end
end


%% Create Path Data
TestPath.T = Pa(:,1);
TestPath.X = Pa(:,2);
TestPath.Y  = Pa(:,3);
TestPath.VX = Pa(:,5).*cos(Pa(:,4));
TestPath.VY = Pa(:,5).*sin(Pa(:,4));

%% Keep Path Label Result
MS   = [TestPath.T   zeros(length(TestPath.T),1)   TestPath.X  TestPath.Y   TestPath.VX  TestPath.VY]; 
for h=1:size(ms_a,1)
    [~,ind]= min(abs(TestPath.T-ms_a(h,1)));
    MS(ind,2) = ms_a(h,2);
end
MSA = MS;
TestPathA = TestPath;

%% session B
ms_a   = []; 
ex_ind = 0;
for i=1:length(CellData{1})
    temp = CellData{sub_session_b}{i};
    if ~isempty(temp)
        for j=1:length(temp)
            ex_ind   = ex_ind + 1;
            if ~isempty(temp{j})
                if ~isempty(temp{j}.data)
                    ex_temp  = temp{j}.data(:,1); 
                    ms_a     = [ms_a; ex_temp   ex_ind*ones(size(ex_temp,1),1)]; 
                end
            end
        end
    end
end
num_of_cell = ex_ind;

%% Create Path Data
TestPath.T = Pb(:,1);
TestPath.X = Pb(:,2);
TestPath.Y = Pb(:,3);
TestPath.VX = Pb(:,5).*cos(Pb(:,4));
TestPath.VY = Pb(:,5).*sin(Pb(:,4));

%% Keep Path Label Result
MS   = [TestPath.T   zeros(length(TestPath.T),1)   TestPath.X  TestPath.Y   TestPath.VX  TestPath.VY]; 
for h=1:size(ms_a,1)
    [~,ind]= min(abs(TestPath.T-ms_a(h,1)));
    MS(ind,2) = ms_a(h,2);
end
MSB = MS;
TestPathB =TestPath;

%% Build plot
fig_ind = 1;
for i=1:num_of_cell
        
    subplot(1,3,1)
    ind = find(MSA(:,2)==i);
    plot(TestPathA.X,TestPathA.Y,'COLOR',[0.8 0.8 0.8]);
    hold on
    plot(MSA(ind,3),MSA(ind,4),'r*','Markersize',12,'LineWidth',2);
    hold off
    xlabel('X');ylabel('Y');axis tight
    grid on;grid minor
    set(gca,'fontsize',18);
    inda=ind;
  %  xlim([150 250])
  %  ylim([50 180])
    axis square

    subplot(1,3,2)
    ind = find(MSB(:,2)==i);
    
    plot(TestPathB.X,TestPathB.Y,'COLOR',[0.6 0.6 0.6]);
    hold on
    plot(MSB(ind,3),MSB(ind,4),'bo','Markersize',12,'LineWidth',2);
    hold off
    xlabel('X');ylabel('Y');axis tight
    grid on;grid minor
    set(gca,'fontsize',18);
    indb = ind;
 %   xlim([150 250])
 %   ylim([50 180])
    axis square
    title(['Cell ' num2str(i)])
    
    
    subplot(1,3,3)
    plot(TestPathA.X,TestPathA.Y,'COLOR',[0.8 0.8 0.8]);
    hold on
    plot(TestPathB.X,TestPathB.Y,'COLOR',[0.6 0.6 0.6]);
    plot(MSA(inda,3),MSA(inda,4),'r*','Markersize',12,'LineWidth',2);
    hold on
    plot(MSB(indb,3),MSB(indb,4),'bo','Markersize',12,'LineWidth',2);
    hold off
    xlabel('X');ylabel('Y');axis tight
    grid on;grid minor
    set(gca,'fontsize',18);
  %  xlim([150 250])
  %  ylim([50 180])
    axis square
    
    
    
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 60 30])
    
    print('-dpng',['Visualization\cell_compare_' num2str(fig_ind)   '.png']);
    fig_ind = fig_ind + 1;
    pause(1)
end


ind_a=1;
ind_b=62;
%%------------------
outputVideo = VideoWriter(fullfile('Visualization','cell_compare.avi'));
outputVideo.FrameRate = 1;
open(outputVideo);
for ii = ind_a:ind_b
   img = imread(['Visualization\cell_compare_' num2str(ii) '.png']);
   writeVideo(outputVideo,img)
end
close(outputVideo);



