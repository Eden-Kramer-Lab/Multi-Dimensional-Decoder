session            = 20;
sub_session_a      = 5;
sub_session_b      = 5;
scale              = 10;

FileName = ['BonWMaze\emilepos' num2str(session) '.mat'];
load(FileName);
PathData     = pos{length(pos)};

%% Mouse Path
P   = PathData{sub_session_a}.data(:,[1 6 7 8 9]);
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
Pa  = P;

%% Place Cell Activity
%FileName = ['BonWMaze\emilespikes' num2str(session) '.mat'];
FileName = ['BonWMaze\spikesnoreplay' num2str(session) '.mat'];
load(FileName);
%CellData     = spikes{1,session};
CellData     = spikesnoreplay{1,session};

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
num_of_cell = ex_ind;

%% Create Path Data
TestPath.T = Pa(:,1);
TestPath.X = Pa(:,2);
TestPath.Y = Pa(:,3);
TestPath.VX = Pa(:,5).*cos(Pa(:,4));
TestPath.VY = Pa(:,5).*sin(Pa(:,4));

%% Keep Path Label Result
MS   = [TestPath.T   zeros(length(TestPath.T),1)   TestPath.X  TestPath.Y   TestPath.VX  TestPath.VY]; 
for h=1:size(ms_a,1)
    [~,ind]= min(abs(TestPath.T-ms_a(h,1)));
    MS(ind,2) = ms_a(h,2);
end

%% Build plot
fig_ind = 1;
for i=1:num_of_cell
    ind = find(MS(:,2)==i);
    
    subplot(2,3,1)
    plot(TestPath.X,TestPath.Y,'.','COLOR',[0.8 0.8 0.8]);
    hold on
    plot(MS(ind,3),MS(ind,4),'*','Markersize',12);
    hold off
    xlabel('X');ylabel('Y');axis tight
    grid on;grid minor
    set(gca,'fontsize',18);

    
    subplot(2,3,4)
    plot(TestPath.VX,TestPath.VY,'.','COLOR',[0.8 0.8 0.8]);
    hold on
    plot(MS(ind,5),MS(ind,6),'*','Markersize',12);
    hold off
    xlabel('VX');ylabel('VY');axis tight
    grid on;grid minor
    ylim([-42 45])
    xlim([-55 55])
    set(gca,'fontsize',18);
    
    subplot(2,3,2)
    plot(TestPath.X,TestPath.VX,'.','COLOR',[0.8 0.8 0.8]);
    hold on
    plot(MS(ind,3),MS(ind,5),'*','Markersize',12);
    hold off
    xlabel('X');ylabel('VX');axis tight
    grid on;grid minor
    ylim([-55 55])
    set(gca,'fontsize',18);
    title(['CELL  ' num2str(i)]);
    
    subplot(2,3,3)
    plot(TestPath.X,TestPath.VY,'.','COLOR',[0.8 0.8 0.8]);
    hold on
    plot(MS(ind,3),MS(ind,6),'*','Markersize',12);
    hold off
    xlabel('X');ylabel('VY');axis tight
    grid on;grid minor
    ylim([-42 45])
    set(gca,'fontsize',18);
    
    subplot(2,3,5)
    plot(TestPath.Y,TestPath.VX,'.','COLOR',[0.8 0.8 0.8]);
    hold on
    plot(MS(ind,4),MS(ind,5),'*','Markersize',12);
    hold off
    xlabel('Y');ylabel('VX');axis tight
    grid on;grid minor
    ylim([-55 55])
    set(gca,'fontsize',18);
    
    subplot(2,3,6)
    plot(TestPath.Y,TestPath.VY,'.','COLOR',[0.8 0.8 0.8]);
    hold on
    plot(MS(ind,4),MS(ind,6),'*','Markersize',12);
    hold off
    xlabel('Y');ylabel('VY');axis tight
    grid on;grid minor
    ylim([-42 45])
    set(gca,'fontsize',18);
    
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 60 30])
    
    print('-dpng',['Visualization\cell_compare_' num2str(fig_ind)   '.png']);
    fig_ind = fig_ind + 1;
    
end


% ind_a=1;
% ind_b=98;
% %%------------------
% outputVideo = VideoWriter(fullfile('images','cell_visulaization.avi'));
% outputVideo.FrameRate = 1;
% open(outputVideo);
% for ii = ind_a:ind_b
%    img = imread(['Visualization\cell_' num2str(ii) '.png']);
%    writeVideo(outputVideo,img)
% end
% close(outputVideo);



