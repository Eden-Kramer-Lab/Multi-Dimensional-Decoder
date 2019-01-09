METHOD = 2;
%% Decoder-Encoder Loop
ind_a = 2;%423000; 
ind_b = 27000;

%% Load Data For Intensity Calculation 
% load Intensity Calculation Components - this bring path as well
load('MarkIntensityNew');
% ------- Complete Model
EXTND = 30;
if METHOD > 1
    mx = min(TestPath.X)-EXTND:1:max(TestPath.X)+EXTND;
    my = min(TestPath.Y)-EXTND:2:max(TestPath.Y)+EXTND;
end

IMAGE_X = zeros(ind_b-ind_a+1,length(mx));
IMAGE_Y = zeros(ind_b-ind_a+1,length(my));
for t=ind_a:1:ind_b
    % run filter
    if METHOD==1
        file_name = ['images\exact_' num2str(t)];
    end
    if METHOD==2
        file_name = ['images\gmm_' num2str(t)];
    end
    if METHOD==3
        file_name = ['images\gauss_' num2str(t)];
    end
    load(file_name);
    if METHOD==1
        %% exact solution (mx, my are loaded)
        Lx = sum(PostComplete,1);
        Ly = sum(PostComplete,2);
    else
        %% other methods
        Lx =ay_gmm_1d(PostMixedModelA,mx,1);
        Ly =ay_gmm_1d(PostMixedModelA,my,2);
    end
    IMAGE_X(t-ind_a+1,:)=Lx;
    IMAGE_Y(t-ind_a+1,:)=Ly;
end
%% subplot 1
subplot(2,1,1)
imagesc(ind_a:ind_b,mx,sqrt(IMAGE_X'));
hold on;
plot(ind_a:ind_b,TestPath.X(ind_a:ind_b),'w','LineWidth',2)
hold off;
xlabel('time step')
ylabel('X (cm)')
%% subplot 2
subplot(2,1,2)
imagesc(ind_a:ind_b,my,sqrt(IMAGE_Y'));
hold on;
plot(ind_a:ind_b,TestPath.Y(ind_a:ind_b),'w','LineWidth',2)
hold off;
xlabel('time step')
ylabel('Y (cm)')
%% call magnify
magnify



