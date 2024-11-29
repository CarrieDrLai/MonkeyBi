function [] = PSTHRasterGo_interp(data,date,session,area,area_ch,data_dir,save_autoselectedtrial_dir,savedir_raster,threshFR)

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
    
    %% plot example raster (after interp)
    interp_num0 = 20*[0 0.25,1.25,1.25,0.25,0.5]; %prepareºÍrewardËãÒ»°ë
    interp_raster = zeros([size(nid,1) size(selected_trial_sort,1),sum(interp_num0)]);
    interp_trajectory = zeros([size(selected_trial_sort,1),sum(interp_num0),size(selected_trajectory_map_sort,3)]);
    judge_len = 6;
    phase_list = [3,3,4,7,8];
    
    for i = 1:size(selected_trial_sort,1)
        for phase = 1:5
            id_curr = find(data.Data.trial_info.trial_no == (selected_trial_sort(i)) & data.Data.trial_info.phase_no == phase_list(phase));
            if phase == 3
                id_curr = find(data.Data.trial_info.trial_no == (selected_trial_sort(i)) & data.Data.trial_info.phase_no == phase_list(phase));
                id_curr = [id_curr,find(data.Data.trial_info.trial_no == (selected_trial_sort(i)) & data.Data.trial_info.phase_no == 6)];
            end
            if size(id_curr,2)>0
                if phase == 1
                    id_curr = id_curr(1:sum(selected_timemap_sort(i,:,1)));
                end
                if phase == 2
                    id_curr = id_curr(sum(selected_timemap_sort(i,:,1))+1:end);
                end
            end
            for n = 1:size(nid,1)        
                curr_nid = find(data.Data.spike.readNeurons == DataRasterGo.RasterGo.NeuronNo(n));
                if size(id_curr,2)<judge_len
                    if phase==1 || phase==3
                        id_next1 = find(data.Data.trial_info.trial_no == (selected_trial_sort(i)) & data.Data.trial_info.phase_no == (phase_list(phase+1)));
                        id_next2 = find(data.Data.trial_info.trial_no == (selected_trial_sort(i)) & data.Data.trial_info.phase_no == (phase_list(phase+2)));
                        id_next = [id_next1,id_next2];
                        id_curr= [id_curr,id_next(1:judge_len-size(id_curr,2))];  
                    end
                    
                    if phase==2 
                        id_next1 = find(data.Data.trial_info.trial_no == (selected_trial_sort(i)) & data.Data.trial_info.phase_no == 4);
                        id_next2 = find(data.Data.trial_info.trial_no == (selected_trial_sort(i)) & data.Data.trial_info.phase_no == 6);
                        id_next3 = find(data.Data.trial_info.trial_no == (selected_trial_sort(i)) & data.Data.trial_info.phase_no == 7);
                        id_next = [id_next1,id_next2,id_next3];
                        id_curr= [id_next(end-(judge_len-1-size(id_curr,2)):end),id_curr];
                    end
                    
                    if phase==4 
                        id_last1 = find(data.Data.trial_info.trial_no == (selected_trial_sort(i)) & data.Data.trial_info.phase_no == 6);
                        id_last2 = find(data.Data.trial_info.trial_no == (selected_trial_sort(i)) & data.Data.trial_info.phase_no == 4);
                        id_last3 = find(data.Data.trial_info.trial_no == (selected_trial_sort(i)) & data.Data.trial_info.phase_no == 3);
                        id_last = [id_last1,id_last2,id_last3];
                        id_curr= [id_last(end-(judge_len-1-size(id_curr,2)):end),id_curr];
                    end        
                    
                    if phase==5 
                        id_last1 = find(data.Data.trial_info.trial_no == (selected_trial_sort(i)) & data.Data.trial_info.phase_no == 7);
                        id_last2 = find(data.Data.trial_info.trial_no == (selected_trial_sort(i)) & data.Data.trial_info.phase_no == 6);
                        id_last3 = find(data.Data.trial_info.trial_no == (selected_trial_sort(i)) & data.Data.trial_info.phase_no == 4);
                        id_last = [id_last1,id_last2,id_last3];
                        id_curr= [id_last(end-(judge_len-1-size(id_curr,2)):end),id_curr];
                    end
                end
                
                raster0 = data.Data.spike.spikes_50Hz(curr_nid,id_curr);
                
                x = (1:size(raster0,2))'; 
                xi = (1:interp_num0(phase+1))';   
                data_interp = interp1(x,raster0,xi,'nearest','extrap')';
                interp_raster(n,i,(1+sum(interp_num0(1:phase))):sum(interp_num0(1:(1+phase)))) =data_interp;
            end
            [a,intersect_id] = intersect(find(data.Data.trial_info.trial_no == (selected_trial_sort(i))),id_curr);
            
%             target_pos0 = 0.24*[0,0;0,1;1,1;1,0;1,-1;0,-1;-1,-1;-1,0;-1,1;];
%             ratio = 0.04;
%             dist = max(abs(target_pos0'))+ratio.*[1,1,1,1,1,1,1,1,1];
% 			B = rem((selected_target_list_sort(i)-rem(selected_target_list_sort(i),10))/10,10);
% 			C = rem(selected_target_list_sort(i),10);
%             max_dist = max(dist(B+1),dist(C+1));
            v_num = size(data.Data.Analog.joystick.joystick_trials,1)/size(data.Data.Analog.joystick.go_last_time,1);
            curr_trialno = selected_trial_sort(i);
            %curr_trial_len = size(find(data.Data.trial_info.trial_no == (selected_trial_sort(i))),2);            
            load([save_autoselectedtrial_dir date '.TrialNo.' num2str(curr_trialno) '.TargetNo.' num2str(selected_target_list_sort(i)) '.mat'])
            curr_joystick = Data_triali.trajectory_i.traject_trial';
            
            %exchange_id = [3,4,1,2];
            
            for vi = 1 :v_num
                trajectory0 = curr_joystick(vi,intersect_id);    
%                 if selected_target_list_sort(i)>400
%                     if B*C==0
%                         trajectory_trials_voltage_x = curr_joystick(exchange_id(vi*2-1),intersect_id)';
%                         trajectory_trials_voltage_y = curr_joystick(exchange_id(vi*2),intersect_id)';
%                         trajectory0_x = 0.5-max_dist*(trajectory_trials_voltage_x-0.5)./0.5;
%                         trajectory0_y = 0.5+max_dist*(trajectory_trials_voltage_y-0.5)./0.5
%                     else
%                         trajectory0 = curr_joystick(vi,intersect_id);    
%                     end
%                 end
                
%             x_trials_voltage = joystick_trials_i(1,:);
%             y_trials_voltage = joystick_trials_i(2,:)';
%             x_trials_voltage_r = joystick_trials_i(3,:)'; 
%             y_trials_voltage_r = joystick_trials_i(4,:)';
%             x_trials = 
%             y_trials = 0.5+max_dist*(y_trials_voltage-0.5)./0.5;
%             x_trials_r = 0.5-max_dist*(x_trials_voltage_r-0.5)./0.5;   % x_trials_r = 0.5 +- (-1,1) max_dist_r
%             y_trials_r = 0.5+max_dist*(y_trials_voltage_r-0.5)./0.5;
%                 
%             if B*C==0  % Unimanual    
%                 curr_x_trials = x_trials_r - 0.24;
%                 curr_y_trials = y_trials_r;
%                 x_trials_r = x_trials + 0.24;
%                 y_trials_r = y_trials;
%                 x_trials = curr_x_trials;
%                 y_trials = curr_y_trials;
%             else
%                 x_center = x_center - 0.24;    
%                 x_return = x_return - 0.24;
%                 x_trials = x_trials - 0.24;
%                 x_center_r = x_center_r + 0.24;
%                 x_return_r = x_return_r + 0.24;
%                 x_trials_r = x_trials_r + 0.24;
%             end
            
                trajectory0 = trajectory0';
                x = (1:size(trajectory0,1))'; 
                xi = (1:interp_num0(phase+1))';   
                data_interp = interp1(x,trajectory0,xi,'nearest','extrap')';
                interp_trajectory(i,(1+sum(interp_num0(1:phase))):sum(interp_num0(1:(1+phase))),vi) = data_interp;
            end
        end
    end

    RasterGo = struct('raster_go',selected_raster_sort,'before_go_binlen',DataRasterGo.RasterGo.before_go_binlen,'psth_binlen',DataRasterGo.RasterGo.psth_binlen,'time_map_go',selected_timemap_sort,'NeuronNo',DataRasterGo.RasterGo.NeuronNo(nid),'area_ch',area_ch,'area',area);
    TargetGo = struct('target_num',target_num,'target_list',selected_target_list_sort,'target_unique',DataRasterGo.TargetGo.target_unique,'trial_no',selected_trial_sort);
    FiringRateGo = struct('mean_fr_session',DataRasterGo.FiringRateGo.mean_fr_session(nid));
    TrajectoryGo = struct('trajectory_map',selected_trajectory_map_sort);
    InterpData = struct('interp_raster',interp_raster,'interp_num0',interp_num0,'interp_trajectory',interp_trajectory);
    
    DataRasterGo_ManualSelected = struct('InterpData',InterpData,'RasterGo',RasterGo,'TargetGo',TargetGo,'FiringRateGo',FiringRateGo,'TrajectoryGo',TrajectoryGo);
    
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
			yticklabel_list = ['LT  ';'L TR';'LR  ';'L BR';'LB  ';'  RT';'R TL';'  RL';'R BL';'  RB';'T  T';'TRTL';'LRRL';'BRBL';'B  B';];
        end
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
    
    %% 3. draw new raster - interp
    [' ======= Start draw Raster -- Interp ======= ']
    load([savedir_raster date '.DataRasterGo_ManualSelected.mat'],'DataRasterGo_ManualSelected');
    target_num = DataRasterGo_ManualSelected.TargetGo.target_num;
    target_unique = DataRasterGo_ManualSelected.TargetGo.target_unique;
    mean_fr = DataRasterGo_ManualSelected.FiringRateGo.mean_fr_session;
     
    for n = 1 : size(mean_fr,1)
        clf;    
        fig = figure('Visible', 'off');
        set(fig, 'Position', [1,1,960,960]);
        %
        %fig = figure(1);
        raster_curr = DataRasterGo_ManualSelected.InterpData.interp_raster(n,:,:);
        raster_curr = reshape(raster_curr,[size(DataRasterGo_ManualSelected.InterpData.interp_raster,2),size(DataRasterGo_ManualSelected.InterpData.interp_raster,3)]);
        imagesc(raster_curr>0);
        hold on;
        for i = 1:size(DataRasterGo_ManualSelected.InterpData.interp_num0,2)-2
            plot([sum(DataRasterGo_ManualSelected.InterpData.interp_num0(1:i+1))+0.5,sum(DataRasterGo_ManualSelected.InterpData.interp_num0(1:i+1))+0.5],[0,sum(target_num)],'-r','LineWidth',3);
            hold on;
        end
        
        ytick = [];
        for t = 1:size(target_unique,1)-1
            plot([0,sum(DataRasterGo_ManualSelected.InterpData.interp_num0)],[sum(target_num(1:t))+0.5,sum(target_num(1:t))+0.5],'-g','LineWidth',3);
            if t >1
                ytick = [ytick;(sum(target_num(1:t))+sum(target_num(1:t-1)))*0.5];
            else
                ytick=[ytick;target_num(1)*0.5];
            end
            hold on;
        end
        ytick = [ytick;(sum(target_num(1:end-1))+sum(target_num))*0.5];

        title([date ' CH' num2str(ceil(double(DataRasterGo_ManualSelected.RasterGo.NeuronNo(n))/6)) ' NeuronNo.' num2str(DataRasterGo_ManualSelected.RasterGo.NeuronNo(n)) ' MeanFR' num2str(mean_fr(n))]);
        xtick = [];
        for i = 1:size(DataRasterGo_ManualSelected.InterpData.interp_num0,2)-2
            xtick = [xtick;sum(DataRasterGo_ManualSelected.InterpData.interp_num0(1:i+1))];
        end
        set(gca,'xtick',xtick);
        xticklabel = [' MoveStart ';'ReturnStart';'RewardStart';' RestStart '];
        set(gca,'xticklabel',xticklabel);
        
        set(gca,'ytick',ytick);
        yticklabel_list = [' T';'TR';' R';'BR';' B';'BL';' L';'TL'];
		if target_unique(1)>400
			yticklabel_list = ['LT  ';'L TR';'LR  ';'L BR';'LB  ';'  RT';'R TL';'  RL';'R BL';'  RB';'T  T';'TRTL';'LRRL';'BRBL';'B  B';];
        end
        yticklabel=[];
        for i = 1:size(target_unique,1)
            yticklabel = [yticklabel;yticklabel_list(i,:)];
        end
        set(gca,'yticklabel',yticklabel);
        %set(gca,'yticklabel',[' T';'TR';' R';'BR';' B';'BL';' L';'TL']);
        hold off;
        fig_name = [date ' CH' num2str(ceil(double(DataRasterGo_ManualSelected.RasterGo.NeuronNo(n))/6)) ' NeuronNo.' num2str(DataRasterGo_ManualSelected.RasterGo.NeuronNo(n)) '.Interp.tif']
        saveas(fig, fullfile(savedir_raster, fig_name) );
    end
    [' ======= End draw Raster -- Interp ======= ']
    
    %% 4. plot trajectory
    savedir_trajectory = [savedir_raster 'Trajectory\'];
    if ~exist(savedir_trajectory, 'dir')
        mkdir(savedir_trajectory);
    end

    for t = 1:size(DataRasterGo_ManualSelected.TargetGo.trial_no,1)
        clf;    
        fig = figure('Visible', 'off');
        set(fig, 'Position', [1,1,480,480]);
        plot_data = reshape(DataRasterGo_ManualSelected.InterpData.interp_trajectory(t,:,:),[sum(DataRasterGo_ManualSelected.InterpData.interp_num0),size(DataRasterGo_ManualSelected.InterpData.interp_trajectory,3)]);
        plot(plot_data);
        hold on;
        for i = 1:size(DataRasterGo_ManualSelected.InterpData.interp_num0,2)-2
            plot([sum(DataRasterGo_ManualSelected.InterpData.interp_num0(1:i+1))+0.5,sum(DataRasterGo_ManualSelected.InterpData.interp_num0(1:i+1))+0.5],[0,1],'-r','LineWidth',3);
            hold on;
        end
        xtick = [];
        for i = 1:size(DataRasterGo_ManualSelected.InterpData.interp_num0,2)-2
            xtick = [xtick;sum(DataRasterGo_ManualSelected.InterpData.interp_num0(1:i+1))];
        end
        set(gca,'xtick',xtick);
        xticklabel = [' Move ';'Return';'Reward';' Rest '];
        set(gca,'xticklabel',xticklabel);
        legend(['vxL';'vyL';'vxR';'vyR']);
        fig_name = [date '.TrialNo' num2str(DataRasterGo_ManualSelected.TargetGo.trial_no(t)) '.Target' num2str(DataRasterGo_ManualSelected.TargetGo.target_list(t)) '.Interp.tif']
        saveas(fig, fullfile(savedir_trajectory, fig_name));
        
        clf;    
        fig = figure('Visible', 'off');
        set(fig, 'Position', [1,1,480,480]);
        plot_data = DataRasterGo_ManualSelected.TrajectoryGo.trajectory_map(t,:,:);
        time_curr = sum(DataRasterGo_ManualSelected.RasterGo.time_map_go(t,:,:),3);
        plot_data = plot_data(:,find(time_curr>0),:);
        plot_data = reshape(plot_data,[size(plot_data,2),size(plot_data,3)]);
        plot(plot_data);
        hold on;
        xtick = [];
        for i = 1:size(DataRasterGo_ManualSelected.InterpData.interp_num0,2)-2
            t_curr = sum(sum(DataRasterGo_ManualSelected.RasterGo.time_map_go(t,:,1:i)));
            if i>2
                if t_curr==xtick(i-1)
                    t_curr = t_curr+1;
                end
            end
            xtick = [xtick;t_curr];
        end
        
        for i = 1:size(DataRasterGo_ManualSelected.InterpData.interp_num0,2)-2
            plot([xtick(i)+0.5,xtick(i)+0.5],[0,1],'-r','LineWidth',3);
            hold on;
        end
        set(gca,'xtick',xtick);
        set(gca,'xticklabel',xticklabel);
        legend(['vxL';'vyL';'vxR';'vyR']);
        fig_name = [date '.TrialNo' num2str(DataRasterGo_ManualSelected.TargetGo.trial_no(t)) '.Target' num2str(DataRasterGo_ManualSelected.TargetGo.target_list(t)) '.tif']
        saveas(fig, fullfile(savedir_trajectory, fig_name) );        
    end
    
end