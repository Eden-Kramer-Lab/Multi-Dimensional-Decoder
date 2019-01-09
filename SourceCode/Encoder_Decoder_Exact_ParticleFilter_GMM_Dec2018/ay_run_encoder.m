function ay_run_encoder(which_file_to_call,spike_file_name,session,sub_session,training_percentage,encoder_file_name)
%% For example
% session     = 20;
% sub_session = 2;
% which_file_to_call = 'BonWMaze\emilepos';
% spike_file_name = 'BonWMaze\spikesnoreplay'


%% Which Session to Call
scale       = 11.2;
d_pt        = 0.0336199999999735;

%% Mouse Path
FileName = [which_file_to_call num2str(session) '.mat'];
load(FileName);
PathData = pos{length(pos)};
% fill zero points
P       = PathData{sub_session}.data(:,[1 2 3 4 5]);
ind     = find( ~isnan(P(:,2)) & ~isnan(P(:,3)) & ~isnan(P(:,4)) & ~isnan(P(:,5)) );
P       = P(ind,:);

ang     = P(:,4);
P(:,4)  = cos(ang).* P(:,5);
P(:,5)  = sin(ang).* P(:,5);
%%------------------
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
        
        vx = P(correct_ind(ind_after),4)+ P(correct_ind(ind_before),4);
        vy = P(correct_ind(ind_after),5)+ P(correct_ind(ind_before),5);
        P(blocked_ind(i),4) = vx*0.5;
        P(blocked_ind(i),5) = vy*0.5;
    end
    if isempty(ind_after) && ~isempty(ind_before)
        P(blocked_ind(i),2) = P(correct_ind(ind_before),2);
        P(blocked_ind(i),3) = P(correct_ind(ind_before),3);
        
        P(blocked_ind(i),4) = P(correct_ind(ind_before),4);
        P(blocked_ind(i),5) = P(correct_ind(ind_before),5);
    end
    if ~isempty(ind_after) && isempty(ind_before)
        P(blocked_ind(i),2) = P(correct_ind(ind_after),2);
        P(blocked_ind(i),3) = P(correct_ind(ind_after),3);
        
        P(blocked_ind(i),4) = P(correct_ind(ind_after),4);
        P(blocked_ind(i),5) = P(correct_ind(ind_after),5);
    end
end
lp = size(P,1);
% Training Set
ind = round(lp*training_percentage/100);
Pa = P(1:ind,:);
% Test Set
Pb = P(ind+1:end,:);

%% Place Cell Activity
FileName = [spike_file_name num2str(session) '.mat'];
load(FileName);
CellData    = spikes{1,session};
ms_t   = [];
ex_ind = 0;
for i=1:length(CellData{1})
    temp = CellData{sub_session}{i};
    if ~isempty(temp)
        for j=1:length(temp)
            ex_ind   = ex_ind + 1;
            if ~isempty(temp{j})
                if ~isempty(temp{j}.data)
                    ex_temp  = temp{j}.data(:,1); 
                    ms_t     = [ms_t;
                                ex_temp   ex_ind*ones(size(ex_temp,1),1)]; 
                end
            end
        end
    end
end
num_of_cell = ex_ind;
[~,ind] = sort(ms_t(:,1),'ascend');
ms_t    = ms_t(ind,:);
[~,ind] = min(abs(Pa(end,1)-ms_t(:,1)));
ms_a    = ms_t(1:ind,:);
ms_b    = ms_t(ind:end,:);

%% Kernel
Kernel.Sm = 6*[1 0;
                 0 1];
Kernel.St = 9*[1 0;
                 0 1];
Kernel.Vm = 6.0*[1 0;
                 0 1];

%% Path Trajectory for Test
min_pt = Pb(1,1);   
max_pt = Pb(end,1);
%d_pt   = Pb(2,1)-Pb(1,1);
TestPath.T  = interp1(Pb(:,1),Pb(:,1),min_pt:d_pt/scale:max_pt);
TestPath.X  = interp1(Pb(:,1),Pb(:,2),min_pt:d_pt/scale:max_pt);
TestPath.Y  = interp1(Pb(:,1),Pb(:,3),min_pt:d_pt/scale:max_pt);
TestPath.VX = interp1(Pb(:,1),Pb(:,4),min_pt:d_pt/scale:max_pt);
TestPath.VY = interp1(Pb(:,1),Pb(:,5),min_pt:d_pt/scale:max_pt);

%% Keep Path Label Result
MS   = [TestPath.T'   zeros(length(TestPath.T'),num_of_cell)   TestPath.X'  TestPath.Y' TestPath.VX' TestPath.VY']; 
for h=1:size(ms_b,1)
    [~,ind] = min(abs(TestPath.T-ms_b(h,1)));
    MS(ind,1+ms_b(h,2)) = 1;
end

%% Cells are based on the first half of data
CellPath.T  = Pa(:,1);
CellPath.X  = Pa(:,2);
CellPath.Y  = Pa(:,3);
CellPath.VX = Pa(:,4);
CellPath.VY = Pa(:,5);

for i=1:num_of_cell
    i
    ind = find( ms_a(:,2)== i);
    if ~isempty(ind)
        for j=1:length(ind)
            [~,t_ind]= min(abs(CellPath.T-ms_a(ind(j),1)));
            PCell{i}.X(j)  = CellPath.X(t_ind);
            PCell{i}.Y(j)  = CellPath.Y(t_ind);
            PCell{i}.VX(j) = CellPath.VX(t_ind);
            PCell{i}.VY(j) = CellPath.VY(t_ind);
            PCell{i}.T(j)  = ms_a(ind(j),1);
        end
    else
        PCell{i}.X = [];
        PCell{i}.Y = [];
        PCell{i}.VX = [];
        PCell{i}.VY = [];
        PCell{i}.T = [];
    end
end

save(encoder_file_name,'TestPath','PCell','CellPath','Kernel','MS');

