function [] = InitializationSelectedTrialAuto_9Grid(date,session,savedir_session_dir,save_path,datafile_id,save_autoselectedtrial_dir,save_manualselectedtrial_dir,threshFR,path_judge)
    %%
    load([savedir_session_dir '\' date '.DatafileID' datafile_id '.Session.' session '.mat']);
    load([save_path '.mat']);

    %mean_fr_session = mean(mean_fr_list);
    
    if ~exist(save_autoselectedtrial_dir, 'dir')
        mkdir(save_autoselectedtrial_dir);
    end
    if ~exist(save_manualselectedtrial_dir, 'dir')
        mkdir(save_manualselectedtrial_dir);
    end
    
    %%
    %%%%%%%% filename
    trial_dir = dir(savedir_session_dir);
    trial_no_filename=[];
    filename_id = [];
    for t = 1:size(trial_dir,1)
        trial_no_curr_split = split(trial_dir(t).name,'.');
        if size(trial_no_curr_split,1)==6
            if sum(cell2mat(trial_no_curr_split(6))=='mat')==3
                trial_no_curr = cell2mat(trial_no_curr_split(3));
                trial_no_filename = [trial_no_filename;str2num(trial_no_curr)];
                filename_id = [filename_id;t];
            end
        end
    end


    %%  Raster Each Neuron: Trial * Time
    %%%%%%%%%%%% raster %%%%%%%%%%%%
    % id In all trials;
    
    trial_id_success = trial_no_list(find(flag.success_hold_center_flag ==1 & flag.auto_flag<1 & sum(movement_onset_mask')'>4)); % 20ms/bin
    id_success = find(flag.success_hold_center_flag ==1 & flag.auto_flag<1 & sum(movement_onset_mask')'>4 ); % 20ms/bin
    
    if path_judge=='Y' %% for classification
        judge_move_maxlen=200;
        trial_id_judge = trial_no_list(find(path_eff<3.5 & path_eff>0.5 & sum(movement_onset_mask')'<judge_move_maxlen));
        id_judge = find(path_eff<3.5 & path_eff>0.5 & sum(movement_onset_mask')'<judge_move_maxlen);
        before_go_binlen = 50;% before go: 50 bins; Go:100bins; Return:75bins; Reward & Intertrial: 125
        psth_binlen = 350;
    else               %% for regression
        judge_move_maxlen=250;
        trial_id_judge = trial_no_list(find(sum(movement_onset_mask')'<judge_move_maxlen));
        id_judge = find(sum(movement_onset_mask')'<judge_move_maxlen);;
        before_go_binlen = 100;% before go: 50 bins; Go:100bins; Return:75bins; Reward & Intertrial: 125
        psth_binlen = 500;
    end
    
    
    trial_id_selected = intersect(trial_id_success,trial_id_judge);
    id_selected = intersect(id_success,id_judge);    
    % id In normally recorded trials;
    % success trial; Not passive obs trial;
    % hold judge in the setup system is set as 100ms (5 bins), if the move len of a trial is
    % less than 100ms, it may be a failed trial

    target_list_new = target_list(id_selected);
    target_unique = unique(target_list_new);
	if target_list_new(1)>400
		target_unique = [410;420;430;440;450;401;408;407;406;405;411;428;437;446;455];

    trial_no = [];
    raster = zeros([size(Data.spike.readNeurons,2),length(id_selected),psth_binlen]);
    trajectory_map = zeros([length(id_selected),psth_binlen,4]);
    time_map = zeros([length(id_selected),psth_binlen,5]);
    del_id = [];
    for i = 1:size(id_selected,1)
        load([savedir_session_dir trial_dir(filename_id(find(trial_no_filename == trial_id_selected(i)))).name]);
        %[savedir_session_dir trial_dir(filename_id(find(trial_no_filename == trial_no_list(id(i))))).name]
        return_len = size(find(Data_triali.trial_info_i.phase_no_i==4),2)+size(find(Data_triali.trial_info_i.phase_no_i==6),2);
        reward_len = size(find(Data_triali.trial_info_i.phase_no_i==7),2);
        intertrial_len = size(find(Data_triali.trial_info_i.phase_no_i==8),2);
        if return_len > (psth_binlen - before_go_binlen - reward_len - intertrial_len)
            del_id = [del_id;i];
        else
            trial_no = [trial_no;trial_id_selected(i)];
            move_id = find(movement_onset_mask(id_selected(i),:)>0);
            
            spike_i = Data_triali.spike_i.spikes_50Hz_i;
            trajectory_i = Data_triali.trajectory_i.traject_trial;
            %before Move#-1
            %before_go0 = min(find(Data_triali.trial_info_i.phase_no_i==3)) - 1 + Data_triali.spike_i.onset_bin - 1; 
            before_go0 = min(move_id) - 1;
            start_point = max(before_go_binlen+1-before_go0,1); % Move#1
            end_point = min(size(Data_triali.trial_info_i.phase_no_i,2)-before_go0 + before_go_binlen,psth_binlen);
            start_point_spike = max(before_go0 - before_go_binlen + 1,1);
            end_point_spike = start_point_spike + end_point-start_point;
            %[start_point,end_point,start_point_spike,end_point_spike]
            %end_point_spike = min(end_point-start_point+1,size(spike_i,2));        

            raster(:,i,start_point:end_point) = spike_i(:,start_point_spike:end_point_spike);
            trajectory_map(i,start_point:end_point,:) = trajectory_i(start_point_spike:end_point_spike,:);
            time_map(i,start_point:before_go_binlen,1) = 1; 
            time_map(i,(before_go_binlen+1):(before_go_binlen+size(move_id,2)),2) = 1; 
            time_map(i,(before_go_binlen+size(move_id,2)+1):(before_go_binlen+size(move_id,2)+return_len),3) = 1; 
            time_map(i,(before_go_binlen+size(move_id,2)+return_len+1):(before_go_binlen+size(move_id,2)+return_len+reward_len),4) = 1; 
            %time_map(i,(before_go_binlen+size(move_id,2)+return_len+reward_len+1):(before_go_binlen+size(move_id,2)+return_len+reward_len+intertrial_len),5) = 1; 
            time_map(i,(before_go_binlen+size(move_id,2)+return_len+reward_len+1):end_point,5) = 1; 
        end
    end
    %time_map = [before_go,move_mark,return_len,reward_len,intertrial_len];
    mask1 = zeros([size(raster,2),1]);
    mask1(find(sum(sum(raster,3),1)>10))=1;
    mask2 = ones([size(raster,2),1]);
    mask2(del_id)=0;
    cut_zero_id = find(mask1&mask2);
    
%     for i = 1:size(del_id,1)
%         id_selected(del_id(i)) = [];
%         raster(:,del_id(i),:) = [];
%         trajectory_map(del_id(i),:,:) = [];
%         time_map(del_id(i),:,:) = [];
%         target_list_new(del_id(i)) = [];
%     end

    cut_zero_id = find(sum(sum(raster,3),1)>10);
    raster = raster(:,cut_zero_id,:);
    trajectory_map = trajectory_map(cut_zero_id,:,:);
    time_map = time_map(cut_zero_id,:,:);
    target_list_new = target_list_new(cut_zero_id,:);
    id_selected = id_selected(cut_zero_id);
    %trial_no = trial_no(cut_zero_id);
	
	mean_fr_session_new=[]; %% mean fr, go phase
	time_trial = sum(sum(time_map,3),2);
	for i = 1:size(raster,1)
		mean_fr_session_new = [mean_fr_session_new;50*sum(sum(reshape(raster(i,:,:),[size(raster(i,:,:),2) psth_binlen]).* sum(time_map,3)))./sum(time_trial)];
	end
    
    nid = find(mean_fr_session_new>threshFR);
    new_readNeurons = Data.spike.readNeurons(nid);
    %mean_fr_session_new = mean_fr_session(nid);
    raster = raster(nid,:,:);
    
    for i = 1:size(trial_no,1)
        copyfile([savedir_session_dir trial_dir(filename_id(find(trial_no_filename == trial_no(i)))).name],save_autoselectedtrial_dir);
        copyfile([savedir_session_dir trial_dir(filename_id(find(trial_no_filename == trial_no(i)))).name],save_manualselectedtrial_dir);
        copyfile([savedir_session_dir 'TrialNo.' num2str(trial_no(i)) '; Target.' num2str(target_list_new(i)) '.tif'],save_autoselectedtrial_dir);
        copyfile([savedir_session_dir 'TrialNo.' num2str(trial_no(i)) '; Target.' num2str(target_list_new(i)) '.tif'],save_manualselectedtrial_dir);
    end


    %% reordering the spike_map, according to the target_id 
    target_num = [];
    target_id = [];
    target_map = zeros(size(target_unique,1),length(target_list_new));
    for t = 1:size(target_unique,1)
        count = sum(target_num);
        target_id_curr = find(target_list_new == target_unique(t));
        target_num = [target_num;size(target_id_curr,1)];
        [move_len_sort,newid] = sort(sum(time_map(target_id_curr,:,2)')); 
        target_id = [target_id;target_id_curr(newid)];
        target_map(t,(sum(target_num)-size(target_id_curr,1)+1):sum(target_num))=1;
    end
    new_raster = raster(:,target_id,:);
    new_trajectory_map = trajectory_map(target_id,:,:);
    new_time_map = time_map(target_id,:,:);
    target_list_new = target_list_new(target_id);
    id_selected = id_selected(target_id);
    trial_no = trial_no(target_id);

    RasterGo = struct('raster_go',new_raster,'before_go_binlen',before_go_binlen,'psth_binlen',psth_binlen,'time_map_go',new_time_map,'NeuronNo',new_readNeurons);
    TargetGo = struct('target_num',target_num,'target_list',target_list_new,'target_unique',target_unique,'target_map',target_map,'trial_no',trial_no);
    FiringRateGo = struct('mean_fr_session',mean_fr_session_new(nid));
    TrajectoryGo = struct('trajectory_map',new_trajectory_map);
    DataRasterGo = struct('RasterGo',RasterGo,'TargetGo',TargetGo,'FiringRateGo',FiringRateGo,'TrajectoryGo',TrajectoryGo);

    save([save_autoselectedtrial_dir date '.DataRasterGo.mat'],'DataRasterGo');

end