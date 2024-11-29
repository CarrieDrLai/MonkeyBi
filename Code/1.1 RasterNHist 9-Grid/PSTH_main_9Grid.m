clc;
clear;

%% public
hand = 'NineGrid';
data_dir_root = 'E:\zju\bci analysis\monkey\';
%data_dir_root = 'F:\';
traject_source = 'NS3'; %% traject_source = 'CSV'; %% traject_source = 'NS3';
savedir =  [data_dir_root 'Data.Session\'];

area = ['M1L';'PMd';'M1R'];
area_ch= [96,64,96];
datafile_id = '1';
session = '1';

%% seperate config
date_list = [
    ['2024' '-' '02' '-' '26'];
    ['2024' '-' '03' '-' '04'];
    ['2024' '-' '03' '-' '11'];
    ['2024' '-' '03' '-' '25'];
    ['2024' '-' '03' '-' '26'];]

%% path
for date_i = 5:size(date_list,1)
    date = date_list(date_i,:);
    
    save_selectedtrial_dir = [data_dir_root 'Data.SelectedTrials\'];
    filename = [date '.DatafileID' datafile_id  '.S' session  '.' traject_source];
    data_dir = [save_selectedtrial_dir filename '.Manual\'];
    %data_dir_session = 'E:\qaas\bci analysis\monkey\Data.Session\';
    save_autoselectedtrial_dir = [save_selectedtrial_dir filename '.Auto\'];

    savedir_raster = [data_dir_root 'Data.Raster\' filename '\'];
    if ~exist(savedir_raster, 'dir')
        mkdir(savedir_raster);
    end

    %%
    threshFR = 0.3;
    data = load([savedir '\' hand '.' date '.DatafileID' datafile_id '.Session.' session '.' traject_source '.mat']);
    PSTHRasterGo_interp(data,date,session,area,area_ch,data_dir,save_autoselectedtrial_dir,savedir_raster,threshFR);
    %PSTHRasterGo(date,session,area,area_ch,data_dir,save_autoselectedtrial_dir,savedir_raster,threshFR);
    PSTHRasterReturn(date,session,area,area_ch,data_dir,save_autoselectedtrial_dir,savedir_raster,threshFR);
    %%
    HistGo_interp(date,area,area_ch,savedir_raster);
    HistReturn(date,area,area_ch,savedir_raster);
    MeanMap(date,area,area_ch,savedir_raster);
    
end
