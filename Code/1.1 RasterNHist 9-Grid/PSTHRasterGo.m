function [] = PSTHRasterGo(date,session,area,area_ch,data_dir,save_autoselectedtrial_dir,savedir_raster,threshFR)

    %%% Behavioral Plot

    %%% Mean FR Plot
    %load([data_dir_session date '.Session.' session '.mat']);
    %load([data_dir date '.S' num2str(session) '\' date '.Session.' session '.mat']);

    %% 0. Load all files of ManualSeletedTrials
    trial_dir = dir(data_dir);
    
    file_mat_id = [];     % id in size(trial_dir,1)
    file_trialno_id = []; % trial no
    selected_trial= [];  
    
    for t = 1:size(trial_dir,1)
        trial_no_curr_split = split(trial_dir(t).name,'.');
        if size(trial_no_curr_split,1)==6
            file_mat_id = [file_mat_id;t];
            file_trialno_id = [file_trialno_id;str2num(cell2mat(trial_no_curr_split(3)))];
        end
        trial_no_curr_split = split(trial_dir(t).name,';');
        if size(trial_no_curr_split,1)==2
            file_trialname = cell2mat(trial_no_curr_split(1));
            file_targetname = cell2mat(trial_no_curr_split(2));
            selected_trial = [selected_trial;str2num(file_trialname(9:end))];
        end
    end

    %load([trial_dir(t).folder '\' date '.' file_trialname(9:end) '.TargetNo.' file_targetname(9:end-4) '.mat']);
        
    %% 1. load RasterGo of AutoSeletedTrials, Remap RasterGo -> RasterGo of ManualSeletedTrials
    
    load([save_autoselectedtrial_dir date '.DataRasterGo.mat'])
    
    [selected_trial, ia ,ib] = intersect(selected_trial, DataRasterGo.TargetGo.trial_no);
    
    selected_raster = DataRasterGo.RasterGo.raster_go(:,ib,:);
    selected_trajectory_map = DataRasterGo.TrajectoryGo.trajectory_map(ib,:,:);
    selected_timemap = DataRasterGo.RasterGo.time_map_go(ib,:,:);
    selected_target_list = DataRasterGo.TargetGo.target_list(ib);
 
    
    target_num = [];
    target_id = [];
    for t = 1:size(DataRasterGo.TargetGo.target_unique,1)
        count = sum(target_num);
        target_id_curr = find(selected_target_list == DataRasterGo.TargetGo.target_unique(t));
        target_num = [target_num;size(target_id_curr,1)];
        [move_len_sort,newid] = sort(sum(selected_timemap(target_id_curr,:,2)')); 
        target_id = [target_id;target_id_curr(newid)];
    end
    nid = find(DataRasterGo.FiringRateGo.mean_fr_session>threshFR);
    selected_raster_sort = selected_raster(nid,target_id,:);
    selected_trajectory_map_sort = selected_trajectory_map(target_id,:,:);
    selected_timemap_sort = selected_timemap(target_id,:,:);
    selected_target_list_sort = selected_target_list(target_id);
    selected_trial_sort = selected_trial(target_id);    
    selected_trial_mat_sort = [];
    for i = 1:size(selected_trial,1)
         selected_trial_mat_sort = [selected_trial_mat_sort;file_mat_id(find(file_trialno_id == selected_trial_sort(i)))];
    end       
    
    RasterGo = struct('raster_go',selected_raster_sort,'before_go_binlen',DataRasterGo.RasterGo.before_go_binlen,'psth_binlen',DataRasterGo.RasterGo.psth_binlen,'time_map_go',selected_timemap_sort,'NeuronNo',DataRasterGo.RasterGo.NeuronNo(nid),'area_ch',area_ch,'area',area);
    TargetGo = struct('target_num',target_num,'target_list',selected_target_list_sort,'target_unique',DataRasterGo.TargetGo.target_unique,'trial_no',selected_trial_sort);
    FiringRateGo = struct('mean_fr_session',DataRasterGo.FiringRateGo.mean_fr_session(nid));
    TrajectoryGo = struct('trajectory_map',selected_trajectory_map_sort);
    DataRasterGo_ManualSelected = struct('RasterGo',RasterGo,'TargetGo',TargetGo,'FiringRateGo',FiringRateGo,'TrajectoryGo',TrajectoryGo);
    
    save([data_dir date '.DataRasterGo_ManualSelected.mat'],'DataRasterGo_ManualSelected','selected_trial_mat_sort');
    save([savedir_raster date '.DataRasterGo_ManualSelected.mat'],'DataRasterGo_ManualSelected','selected_trial_mat_sort');
    
    %% 2. draw new raster
    [' ======= Start draw Raster -- Go Phase ======= ']
    load([savedir_raster date '.DataRasterGo_ManualSelected.mat'],'DataRasterGo_ManualSelected');
    before_go_binlen = DataRasterGo_ManualSelected.RasterGo.before_go_binlen;
    psth_binlen = DataRasterGo_ManualSelected.RasterGo.psth_binlen;
    target_num = DataRasterGo_ManualSelected.TargetGo.target_num;
    target_unique = DataRasterGo_ManualSelected.TargetGo.target_unique;
    mean_fr = DataRasterGo_ManualSelected.FiringRateGo.mean_fr_session;
     
    for n = 1 : size(DataRasterGo_ManualSelected.RasterGo.raster_go,1)
        clf;    
        fig = figure('Visible', 'off');
        set(fig, 'Position', [1,1,960,960]);
        %
        %fig = figure(1);
        raster_curr = DataRasterGo_ManualSelected.RasterGo.raster_go(n,:,:);
        raster_curr = reshape(raster_curr,[size(DataRasterGo_ManualSelected.RasterGo.raster_go,2),size(DataRasterGo_ManualSelected.RasterGo.raster_go,3)]);
        imagesc(raster_curr>0);
        hold on;

        plot([before_go_binlen+0.5,before_go_binlen+0.5],[0,psth_binlen],'-r','LineWidth',3);
        hold on;

        ytick = [];
        for t = 1:size(target_unique,1)-1
            plot([0,psth_binlen],[sum(target_num(1:t))+0.5,sum(target_num(1:t))+0.5],'-g','LineWidth',3);
            if t >1
                ytick = [ytick;(sum(target_num(1:t))+sum(target_num(1:t-1)))*0.5];
            else
                ytick=[ytick;target_num(1)*0.5];
            end
            hold on;
        end

        ytick = [ytick;(sum(target_num(1:end-1))+sum(target_num))*0.5];
        for i = 1:sum(target_num)
            plot(min(find(selected_timemap_sort(i,:,1)>0)),i,'go','Markersize',3,'LineWidth',2);
            hold on;
            plot(max(find(selected_timemap_sort(i,:,2)>0)),i,'ro','Markersize',5,'LineWidth',2);
            hold on;
            plot(max(find(selected_timemap_sort(i,:,3)>0)),i,'go','Markersize',5,'LineWidth',2);
            hold on;
            if size(max(find(selected_timemap_sort(i,:,4)>0)),2)==0
                plot(max(find(selected_timemap_sort(i,:,3)>0))+1,i,'ro','Markersize',3,'LineWidth',2);
            else
                plot(max(find(selected_timemap_sort(i,:,4)>0)),i,'ro','Markersize',3,'LineWidth',2);
            end
            hold on;
            if size(max(find(selected_timemap_sort(i,:,5)>0)),2)==0
                if size(max(find(selected_timemap_sort(i,:,4)>0)),2)>0
                    plot(max(find(selected_timemap_sort(i,:,4)>0))+1,i,'go','Markersize',3,'LineWidth',2);
                else
                    plot(max(find(selected_timemap_sort(i,:,3)>0))+2,i,'go','Markersize',3,'LineWidth',2);
                end
            else
                plot(max(find(selected_timemap_sort(i,:,5)>0)),i,'go','Markersize',3,'LineWidth',2);
            end
        end

        title([date ' CH' num2str(ceil(double(DataRasterGo_ManualSelected.RasterGo.NeuronNo(n))/6)) ' NeuronNo.' num2str(DataRasterGo_ManualSelected.RasterGo.NeuronNo(n)) ' MeanFR' num2str(mean_fr(n))]);
        xlabel = [before_go_binlen];
        xtick = 'on';
        set(gca,'xtick',before_go_binlen);
        set(gca,'xticklabel','start');
        set(gca,'ytick',ytick);
        yticklabel_list = [' T';'TR';' R';'BR';' B';'BL';' L';'TL'];
		if target_unique(1)>400
			yticklabel_list = ['LT  ';'L TR';'LR  ';'L BR';'LB  ';'  RT';'R TL';'  RL';'R BL';'  RB';'T  T';'TRTR';'LRRL';'BRBR';'B  B';];

        yticklabel=[];
        for i = 1:size(target_unique,1)
            yticklabel = [yticklabel;yticklabel_list(i,:)];
        end
        set(gca,'yticklabel',yticklabel);
        %set(gca,'yticklabel',[' T';'TR';' R';'BR';' B';'BL';' L';'TL']);
        hold off;
        fig_name = [date ' CH' num2str(ceil(double(DataRasterGo_ManualSelected.RasterGo.NeuronNo(n))/6)) ' NeuronNo.' num2str(DataRasterGo_ManualSelected.RasterGo.NeuronNo(n)) '.tif']
        saveas(fig, fullfile(savedir_raster, fig_name) );
    end
    [' ======= End draw Raster -- Go Phase ======= ']

end