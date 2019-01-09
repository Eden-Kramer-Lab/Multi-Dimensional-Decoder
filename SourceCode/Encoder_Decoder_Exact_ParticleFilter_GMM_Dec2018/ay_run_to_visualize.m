%% Which Session to Call
session            = 10;
sub_session_a      = 2;
sub_session_b      = 4;
scale              = 1;66.8000000005122/10;

%% Mouse Path
FileName = ['BonWMaze\bonpos' num2str(session) '.mat'];
load(FileName);
PathData     = pos{length(pos)};
% fill zero points
P   = PathData{sub_session_a}.data(:,[1 2 3]);
blocked_ind = find(P(:,2)==0);
correct_ind = find(P(:,2)~=0);
for i=1:length(blocked_ind)
    ind_after  = find(correct_ind >= blocked_ind(i));
    ind_after  = ind_after(1);
    
    ind_before = find(correct_ind <= blocked_ind(i));
    ind_before = ind_before(end);
    
%     if ~isempty(ind_after) && ~isempty(ind_before)
%         dT = P(correct_ind(ind_after),1)- P(correct_ind(ind_before),1); 
%         dX = P(correct_ind(ind_after),2)- P(correct_ind(ind_before),2);
%         dY = P(correct_ind(ind_after),3)- P(correct_ind(ind_before),3);
%         dt = P(blocked_ind(i),1) - P(correct_ind(ind_before),1);
%         xp = P(correct_ind(ind_after),2)+ (dX * dt)/dT;
%         yp = P(correct_ind(ind_after),3)+ (dY * dt)/dT;
%         P(blocked_ind(i),2) = xp;
%         P(blocked_ind(i),3) = yp;
%     end
    if ~isempty(ind_after) && ~isempty(ind_before)
       % dT = P(correct_ind(ind_after),1)- P(correct_ind(ind_before),1); 
       % dX = P(correct_ind(ind_after),2)- P(correct_ind(ind_before),2);
       % dY = P(correct_ind(ind_after),3)- P(correct_ind(ind_before),3);
       % dt = P(blocked_ind(i),1) - P(correct_ind(ind_before),1);
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

P   = PathData{sub_session_b}.data(:,[1 2 3]);
blocked_ind = find(P(:,2)==0);
correct_ind = find(P(:,2)~=0);
for i=1:length(blocked_ind)
    ind_after  = find(correct_ind >= blocked_ind(i));
    ind_after  = ind_after(1);
    
    ind_before = find(correct_ind <= blocked_ind(i));
    ind_before = ind_before(end);
    
%     if ~isempty(ind_after) && ~isempty(ind_before)
%         dT = P(correct_ind(ind_after),1)- P(correct_ind(ind_before),1); 
%         dX = P(correct_ind(ind_after),2)- P(correct_ind(ind_before),2);
%         dY = P(correct_ind(ind_after),3)- P(correct_ind(ind_before),3);
%         dt = P(blocked_ind(i),1) - P(correct_ind(ind_before),1);
%         xp = P(correct_ind(ind_after),2)+ (dX * dt)/dT;
%         yp = P(correct_ind(ind_after),3)+ (dY * dt)/dT;
%         P(blocked_ind(i),2) = xp;
%         P(blocked_ind(i),3) = yp;
%     end
    if ~isempty(ind_after) && ~isempty(ind_before)
       % dT = P(correct_ind(ind_after),1)- P(correct_ind(ind_before),1); 
       % dX = P(correct_ind(ind_after),2)- P(correct_ind(ind_before),2);
       % dY = P(correct_ind(ind_after),3)- P(correct_ind(ind_before),3);
       % dt = P(blocked_ind(i),1) - P(correct_ind(ind_before),1);
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
Pb = P;

%% We have two figure




%% Place Cell Activity
FileName = ['BonWMaze\bonspikes' num2str(session) '.mat'];
load(FileName);
CellData     = spikes{1,session};

ms_a   = []; 
ex_ind = 0;
for i=1:length(CellData{1})
    temp = CellData{sub_session_a}{i};
    if ~isempty(temp)
        for j=1:length(temp)
            ex_ind = ex_ind+1;
            if ~isempty(temp{j})
                if ~isempty(temp{j}.data)
                    ex_temp  = temp{j}.data(:,1:3); 
                    ms_a     = [ms_a; ex_temp   ex_ind*ones(size(ex_temp,1),1)]; 
                end
            end
        end
    end
end
num_of_cell_a = ex_ind;

ms_b   = [];
ex_ind = 0;
for i=1:length(CellData{1})
    temp = CellData{sub_session_b}{i};
    if ~isempty(temp)
        for j=1:length(temp)
            ex_ind   = ex_ind + 1;
            if ~isempty(temp{j})
                if ~isempty(temp{j}.data)
                    ex_temp  = temp{j}.data(:,1:3); 
                    ms_b     = [ms_b; ex_temp   ex_ind*ones(size(ex_temp,1),1)]; 
                end
            end
        end
    end
end
num_of_cell_b = ex_ind;

for h=1: 1;%num_of_cell_a
  figure(1)
  subplot(10,10,h)
  plot(Pa(:,2),Pa(:,3),'.','LineWidth',0.2,'Color',[0.8 0.8 0.8]);hold on;
  % plot(Pb(:,2),Pb(:,3),'o','Color',[0.1 0.1 0.1]);
  axis tight
  ind= find(ms_a(:,4)==h);
  plot(ms_a(ind,2),ms_a(ind,3),'*','LineWidth',0.2);

  figure(2)
  subplot(10,10,h)
  % plot(Pa(:,2),Pa(:,3),'.','Color',[0.8 0.8 0.8]);hold on;
  plot(Pb(:,2),Pb(:,3),'.','LineWidth',0.2,'Color',[0.8 0.8 0.8]);hold on
  axis tight
  ind= find(ms_b(:,4)==h);
  plot(ms_b(ind,2),ms_b(ind,3),'*','LineWidth',0.2);
end

min_pt_pa = Pa(1,1);   max_pt_pa = Pa(end,1); d_pt_pa   = Pb(2,1)-Pb(1,1);%Pa(2,1)-Pa(1,1);
min_pt_pb = Pb(1,1);   max_pt_pb = Pb(end,1); d_pt_pb   = Pb(2,1)-Pb(1,1);

fig_ind=1;
for h=1: num_of_cell_a
  figure(h)
  subplot(2,2,1)   
  plot(Pa(:,2),'r.');hold on;
  title(['Pa - cell:'  num2str(h)])
  axis tight
  ylabel('X')
  t_ind = find(ms_a(:,end)==h);
  t_ind = round(1+(ms_a(t_ind,1)-min_pt_pa)/d_pt_pa);
  plot(t_ind,Pa(t_ind,2),'b+'); hold off;
  
  subplot(2,2,3)   
  plot(Pa(:,3),'r.');hold on;
  ylabel('Y')
  plot(t_ind,Pa(t_ind,3),'b+'); hold off;
  axis tight
  
  %-----------------------------------
  
  
  subplot(2,2,2)   
  plot(Pb(:,2),'r.');hold on;
  title(['Pb - cell:'  num2str(h)])
  ylabel('X')
  t_ind = find(ms_b(:,end)==h);
  t_ind = round(1+(ms_b(t_ind,1)-min_pt_pb)/d_pt_pb);
  plot(t_ind,Pb(t_ind,2),'b+'); hold off;
  axis tight
  
  
  subplot(2,2,4)   
  plot(Pb(:,3),'r.');hold on;
  ylabel('Y')
  plot(t_ind,Pb(t_ind,3),'b+'); hold off;
  axis tight
  
  %----------------------------------
  
  set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 60 10])
  print('-dpng',['Visualization\path_' num2str(fig_ind) '.png']);
  fig_ind = fig_ind + 1;
  axis tight  
end

ind_a = 14550;
ind_b = 14850;
for h=1: num_of_cell_a
  figure(h)
  %-----------------------------------
  subplot(2,1,1)   
  plot(Pb(:,2),'r.');hold on;
  title(['Pb - cell:'  num2str(h)])
  ylabel('X')
  t_ind = find(ms_b(:,end)==h);
  t_ind = round(1+(ms_b(t_ind,1)-min_pt_pb)/d_pt_pb);
  plot(t_ind,Pb(t_ind,2),'b+'); hold off;
  axis tight
  xlim([ind_a ind_b])
  
  
  subplot(2,1,2)   
  plot(Pb(:,3),'r.');hold on;
  ylabel('Y')
  plot(t_ind,Pb(t_ind,3),'b+'); hold off;
  axis tight
  xlim([ind_a ind_b])
  
  %----------------------------------
  pause() 
  
end



