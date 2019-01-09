%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Encoder Function %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data files are in BonWMaze folder
%% There are two files per each session and a rat
%% For example, for the rat called "emile" in session "20" there is a "positon file"  called "emilepos20"
%% For this session, there is also a "spike" files which is called "emilespikes20"
%% Each session consists of multiple subsessions, here we call sub_session 4 from session 20 using "emile" data
run_encoder          = 1;
which_file_to_call   = 'rat_cell_data\bonpos';
spike_file_name      = 'rat_cell_data\bonspikes';
session              = 10;
sub_session          = 4;
training_percentage  = 85; % training data is the first 85% of the data
encoder_file_name    = 'encoder_files\emile';  % this will be the data file being generated for the decoder model 
if run_encoder ==1
    ay_run_encoder(which_file_to_call,spike_file_name,session,sub_session,training_percentage,encoder_file_name);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Decoder Function %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
which_method     = 5; % 1= exact, 2= particle filtering(2d) 3= single gaussian(2d) 4= gmm (2d) 5=gmm(4d)
performance_file = 'decoder_performance\emile_perf';
detail_base_file = 'decoder_per_timestep\emile_files';
ay_run_decoder(encoder_file_name,which_method,performance_file,detail_base_file);


