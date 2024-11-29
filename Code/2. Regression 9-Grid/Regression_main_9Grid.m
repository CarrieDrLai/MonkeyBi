clc;
clear;

%% public
hand = 'NineGrid';
data_dir_root = 'E:\zju\bci analysis\monkey\';
save_selectedtrial_dir = [data_dir_root 'Data.SelectedTrials\'];
traject_source = 'NS3'; %% traject_source = 'CSV'; %% traject_source = 'NS3';

savedir =  [data_dir_root 'Process.Decode\'];
%% seperate config

date_list = [
    ['2024' '-' '02' '-' '26'];
    ['2024' '-' '03' '-' '04'];
    ['2024' '-' '03' '-' '11'];
    ['2024' '-' '03' '-' '25'];
    ['2024' '-' '03' '-' '26'];]

%% path



%% Decode
% Neural data: Single Units, Area, Ensemble + Ensemble (selected neurons);
% Context: Within\Across Hand
threshFR = 0.3;

NineGridRegression(savedir,save_selectedtrial_dir,date_list,threshFR)


%% Further Analysis for Decode

date_list = [
    ['2024' '-' '02' '-' '26'];
    ['2024' '-' '03' '-' '04'];
    ['2024' '-' '03' '-' '11'];
    ['2024' '-' '03' '-' '25'];
    ['2024' '-' '03' '-' '26'];]
area_ch= [0,96,64,96];
area = ['M1L';'PMd';'M1R'];

% ensemble
[left_xy,right_xy,bi_LxyRxy,left_xy_shuffle,right_xy_shuffle,bi_LxyRxy_shuffle,p_list,p_hat_list,h_list]=AnalysisEnsembleWithinHandCC(savedir,date_list);
[Uni2Bi_Lxy,Uni2Bi_Rxy,Bi2Uni_LxyRxy,Uni2Bi_Lxy_shuffle,Uni2Bi_Rxy_shuffle,Bi2Uni_LxyRxy_shuffle,p_list_Across,p_hat_list_Across,h_list_Across]=AnalysisEnsembleAcrossHandCC(savedir,date_list);

graphpad_table_L = [left_xy',left_xy_shuffle'];
graphpad_table_R = [right_xy',right_xy_shuffle'];
graphpad_table_B = [bi_LxyRxy',bi_LxyRxy_shuffle'];
p_hat_list;

graphpad_table_Uni2Bi_L = [Uni2Bi_Lxy',Uni2Bi_Lxy_shuffle'];
graphpad_table_Uni2Bi_R = [Uni2Bi_Rxy',Uni2Bi_Rxy_shuffle'];
graphpad_table_Bi2Uni = [Bi2Uni_LxyRxy',Bi2Uni_LxyRxy_shuffle'];
p_hat_list_Across

% units
savedir =  [data_dir_root 'Process.Decode\'];
area = ['M1L';'PMd';'M1R'];
area_ch= [0,96,64,96];

[Within,Across_Uni_2_Bi,Across_Bi_2_Uni,p_list,p_hat_list]=AnalysisUnitsCC(savedir,date_list,area,area_ch);

%[temp_within,temp_across,cc_temp_within,cc_temp_across] = AnalysisEnsembleWithinHandCC(savedir)

cc_thresh =0.1;
type_list = AnalysisClassifyNeurons(savedir,area,cc_thresh);



