clc;
clear;
%% public
hand = 'NineGrid';
data_dir_root = 'E:\zju\bci analysis\monkey\';
%data_dir_root = 'F:\';
traject_source = 'NS3'; %% traject_source = 'CSV'; %% traject_source = 'NS3';
csv_screensize = [1920,1080];
savedir =  [data_dir_root 'Data.Session\'];

%% seperate config
date = ['2024' '-' '02' '-' '26'];
datafile_id = '1';
session = '1';

area = ['M1L';'PMd';'M1R'];
area_ch= [96,64,96];
notmember_list=[];
%data_dir = 'E:\zju\bci analysis\monkey\nsp\';
%notmember_list=[333,561,579,639,939,1071]; % 20221110
%% seperate config
date = ['2024' '-' '03' '-' '04'];

%% seperate config
date = ['2024' '-' '03' '-' '11'];

%% seperate config
date = ['2024' '-' '03' '-' '25'];

%% seperate config
date = ['2024' '-' '03' '-' '26'];

%%
data_dir = [data_dir_root 'nsp\16318ninegrid\' date(1:4) date(6:7) date(9:10) '\'];
data_dir_joystick = [data_dir_root 'joystick\ui\' date '\Session-' session '\'];

%% Blackrock .nev 
%filename_nev = [data_dir 'datafile00' datafile_id '.nev'];
filename_nev = [data_dir 'datafile00' datafile_id '-sorted.nev'];
%filename_nev2 = [data_dir 'datafile00' num2str(datafile_id+1) '.nev'];

%filename_nev = [data_dir 'datafile001-sorted.nev'];
%fn = 'E:\zju\1.18\datafile001-202201228.nev';
%openNEV(filename_nev2, 'nomat', 'nosave');
%NEV2 = NEV;
openNEV(filename_nev, 'nomat', 'nosave');

%% Blackrock .ns3 .ns6
filename_ns6 = [data_dir 'datafile00' datafile_id '.ns6'];
filename_ns3 = [data_dir 'datafile00' datafile_id '.ns3'];
if strcmp(traject_source,'NS3')
    raw_data = openNSx('report','read',filename_ns3, 'p:double', 's:20');
    f=2000;
else
    raw_data = openNSx('report','read',filename_ns6, 'p:double', 's:300');
    f=30000;
end
f_joystick =100; % =f/s
% openNSx('report','read',filename_ns6,'t:0:10','sec', 'p:double', 's:600');
% mode:The user can specify the mode of duration in [duration],such as 'sec', 'min', 'hour', or 'sample'
% precision
% skipfactor

%notmember_list = [187,259,313,379,385,433,445,446,547,571,583,637,643,751,757,758,763,811,812,813,817,937,985,1027,1028,1029,1045,1051,1052,1053,1069,1070];

%notmember_list =[561,573,1047,1071]; 
%notmember_list =[735]; % 20230313
%floor(1+double(readNeurons)/6)
%% Read Data
save_path = [savedir hand '.' date '.DatafileID' datafile_id '.Session.' session '.' traject_source];
[readNeurons,nsp_ch] = InitializationReadNeuron(NEV,notmember_list)
session_csv = csvread([data_dir_joystick date '.Session-' session '.csv.csv'],1,1);
Data = InitializationAnalogInputMain_9Grid(NEV,readNeurons,nsp_ch,raw_data,date,data_dir,data_dir_joystick,session,session_csv,savedir,f,f_joystick,save_path);

%% CutTrial 1 NS6
% load([savedir hand '.' date '.DatafileID' datafile_id '.Session.' session '.NS6' '.mat']);
% save_trial_dir =  [savedir date '.DatafileID' datafile_id '.S' session '.' 'NS6' '\'];
% InitializationCutTrial(Data,date,datafile_id,save_trial_dir,session,traject_source,csv_screensize);
%% CutTrial 2 CSV
% load([savedir hand '.' date '.DatafileID' datafile_id '.Session.' session '.CSV' '.mat']);
% save_trial_dir =  [savedir date '.DatafileID' datafile_id '.S' session '.CSV' '\'];
% InitializationCutTrial(Data,date,datafile_id,save_trial_dir,session,'CSV',csv_screensize);


%% CutTrial 3 NS3 & Plot Trajectory
%%%%%%%%%%%%20240122 判断除了center success,距离多少以内的也可以算对（不一定握在了棒上
save_path = [savedir hand '.' date '.DatafileID' datafile_id '.Session.' session '.' traject_source];
load([save_path '.mat']);
save_trial_dir =  [save_path '\'];
InitializationCutTrial_9Grid(Data,date,datafile_id,save_trial_dir,session,'NS3',csv_screensize); 

%% Selected Neurons & Trials 1 -- Regression
% threshFR = 0.3;
% traject_source='NS3'
% savedir_session_dir =  [savedir date '.DatafileID' datafile_id '.S' session '.' traject_source '\'];
% 
% save_selectedtrial_dir = 'E:\qaas\bci analysis\monkey\Data.SelectedTrials.NoJudge\';
% save_autoselectedtrial_dir = [save_selectedtrial_dir date '.DatafileID' datafile_id '.S' session '.' traject_source '.Auto\'];
% save_manualselectedtrial_dir = [save_selectedtrial_dir date '.DatafileID' datafile_id '.S' session '.' traject_source '.Manual\'];
% if ~exist(save_autoselectedtrial_dir, 'dir')
%     mkdir(save_autoselectedtrial_dir);
% end
% if ~exist(save_manualselectedtrial_dir, 'dir')
%     mkdir(save_manualselectedtrial_dir);
% end
% path_judge = 'N';
% InitializationSelectedTrialAuto(date,session,savedir_session_dir,savedir,hand,datafile_id,save_autoselectedtrial_dir,save_manualselectedtrial_dir,threshFR,traject_source, path_judge);

%% Selected Neurons & Trials 1 -- Classification NS3
%%%%%%%%%%%%20240122 判断除了center success,距离多少以内的也可以算对（不一定握在了棒上
threshFR = 0.3;
traject_source='NS3'
savedir_session_dir =  [savedir hand '.' date '.DatafileID' datafile_id '.Session.' session '.' traject_source '\'];
save_path = [savedir hand '.' date '.DatafileID' datafile_id '.Session.' session '.' traject_source];

save_selectedtrial_dir = [data_dir_root 'Data.SelectedTrials\'];
save_autoselectedtrial_dir = [save_selectedtrial_dir date '.DatafileID' datafile_id '.S' session '.' traject_source '.Auto\'];
save_manualselectedtrial_dir = [save_selectedtrial_dir date '.DatafileID' datafile_id '.S' session '.' traject_source '.Manual\'];

path_judge = 'Y';
InitializationSelectedTrialAuto_9Grid(date,session,savedir_session_dir,save_path,datafile_id,save_autoselectedtrial_dir,save_manualselectedtrial_dir,threshFR,path_judge);

%% Selected Neurons & Trials 2 -- Classification NS6
% threshFR = 0.3;
% traject_source='NS6';
% savedir_session_dir =  [savedir date '.DatafileID' datafile_id '.S' session '.' traject_source '\'];
% 
% save_selectedtrial_dir = 'E:\qaas\bci analysis\monkey\Data.SelectedTrials\';
% save_autoselectedtrial_dir = [save_selectedtrial_dir date '.DatafileID' datafile_id '.S' session '.' traject_source '.Auto\'];
% save_manualselectedtrial_dir = [save_selectedtrial_dir date '.DatafileID' datafile_id '.S' session '.' traject_source '.Manual\'];
% if ~exist(save_autoselectedtrial_dir, 'dir')
%     mkdir(save_autoselectedtrial_dir);
% end
% if ~exist(save_manualselectedtrial_dir, 'dir')
%     mkdir(save_manualselectedtrial_dir);
% end
% path_judge = 'Y';
% InitializationSelectedTrialAuto(date,session,savedir_session_dir,savedir,hand,datafile_id,save_autoselectedtrial_dir,save_manualselectedtrial_dir,threshFR,traject_source,path_judge);


%% Manual 
