function [] = InitializationCutTrial_9Grid(Data,date,datafile_id,save_trial_dir,session,traject_source,csv_screensize)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%   7 Save Trial Data   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~exist(save_trial_dir, 'dir')
        mkdir(save_trial_dir);
    end

    %trial_id = find((Data.trial_info.flag.success_hold_center_flag==1).*(Data.trial_info.flag.success_hold_flag==1).*(~Data.trial_info.flag.auto_flag));
    trial_list = Data.trial_info.trial_list';
    success_id = (~Data.trial_info.flag.auto_flag);

    target_pos_voltage = [0,0;0,1;0.71,0.71;1,0;0.71,0.71;0,-1;-0.71,-0.71;-1,0;-0.71,0.71];
    target_pos_voltage = target_pos_voltage./2+0.5;
    %target_pos = 0.32*[0,1;0.71,0.71;1,0;0.71,-0.71;0,-1;-0.71,-0.71;-1,0;-0.71,0.71;];
	target_pos0 = 0.24*[0,0;0,1;1,1;1,0;1,-1;0,-1;-1,-1;-1,0;-1,1;];
    ratio = 0.04;
    dist = max(abs(target_pos0'))+ratio.*[1,1,1,1,1,1,1,1,1];
    %dist = max(abs(target_pos'))+ratio.*[1,1.41,1,1.41,1,1.41,1,1.41];
    target_pos = target_pos0+0.5;
	target_pos(:,1) = target_pos(:,1)-0.24;
	target_pos_r = target_pos;
	target_pos_r(:,1) = target_pos_r(:,1)+2*0.24;

    target_list = [];
    mean_fr_list = [];
    mean_center_fr_list = [];
    mean_peripheral_fr_list = [];
    baseline_list = [];
    baseline_sd_list = [];
    time_center_list = [];
    time_peripheral_list = [];
    movement_onset_mask = zeros((size(trial_list,1)),20000);
    trial_no_list = [];
    count = 0;
    spikes_50Hz_raster = zeros((size(trial_list,1)),20000);

    success_hold_center_flag = [];
    success_hold_flag = [];
    repeat_timeout_flag = [];
    repeat_flag = [];
    auto_flag = [];
    
    dist_x_square_list = [];
    dist_y_square_list = [];
    min_sqr_dist_list = [];
    path_eff=[];
    for i = 1: size(trial_list,1)
    % 1. success; 2. phase的时间节点要对； 3. 到达target位置不能偏离过多
        trial_no = trial_list(i); % curr_trial_no
        trial_id = find(Data.trial_info.trial_list==trial_no);

        phase_no_i = Data.trial_info.phase_no(find(Data.trial_info.trial_no == trial_no)); %mark the time of phases
        %if size(unique(phase_no_i),2)<5
        if size(find(phase_no_i==3),2)==0|size(find(phase_no_i==6),2)==0
            continue
        else
            meta_trial_i = Data.meta;
            target_i = Data.trial_info.target_list(trial_id);
            sensor_i = Data.Analog.sensor.sensor_trials(find(Data.trial_info.trial_no == trial_id));
            joystick_hold_i = Data.Analog.joystick.joystick_hold(4*trial_id-3:4*trial_id,:);
            joystick_hold_i = joystick_hold_i(:,(1:max(find(joystick_hold_i(1,:)~=0))));
            joystick_center_hold_i = Data.Analog.joystick.joystick_center_hold(4*trial_id-3:4*trial_id,:);
            joystick_center_hold_i = joystick_center_hold_i(:,(1:max(find(joystick_center_hold_i(1,:)~=0))));
            % joystick_hold_csv_i = Data.Analog.joystick.joystick_hold_csv(2*trial_id-1:2*trial_id,:);
            % joystick_hold_csv_i = joystick_hold_csv_i(:,(1:max(find(joystick_hold_csv_i(1,:)~=0))));
            % joystick_center_hold_csv_i = Data.Analog.joystick.joystick_center_hold_csv(2*trial_id-1:2*trial_id,:);
            % joystick_center_hold_csv_i = joystick_center_hold_csv_i(:,(1:max(find(joystick_center_hold_csv_i(1,:)~=0))));
            joystick_trials_i = Data.Analog.joystick.joystick_trials(4*trial_id-3:4*trial_id,:);   
            joystick_trials_i = joystick_trials_i(:,(1:max(find(joystick_trials_i(1,:)~=0))));
            
			C = fix(rem(rem(target_i,100),10));
			B = fix((rem(target_i,100) - C)/10);
            max_dist = max(dist(B+1),dist(C+1));
            
            x_trials_voltage = joystick_trials_i(1,:)';
            y_trials_voltage = joystick_trials_i(2,:)';
            x_trials_voltage_r = joystick_trials_i(3,:)'; 
            y_trials_voltage_r = joystick_trials_i(4,:)';
            x_trials = 0.5-max_dist*(x_trials_voltage-0.5)./0.5;
            y_trials = 0.5+max_dist*(y_trials_voltage-0.5)./0.5;
            x_trials_r = 0.5-max_dist*(x_trials_voltage_r-0.5)./0.5;   % x_trials_r = 0.5 +- (-1,1) max_dist_r
            y_trials_r = 0.5+max_dist*(y_trials_voltage_r-0.5)./0.5;
			
            x_center_voltage = joystick_center_hold_i(1,:)';
            y_center_voltage = joystick_center_hold_i(2,:)';
            x_center_voltage_r = joystick_center_hold_i(3,:)';
            y_center_voltage_r = joystick_center_hold_i(4,:)';
            x_center = 0.5-max_dist*(x_center_voltage-0.5)./0.5;
            y_center = 0.5+max_dist*(y_center_voltage-0.5)./0.5;
            x_center_r = 0.5-max_dist*(x_center_voltage_r-0.5)./0.5;
            y_center_r = 0.5+max_dist*(y_center_voltage_r-0.5)./0.5;

            x_return_voltage = joystick_hold_i(1,:)';
            y_return_voltage = joystick_hold_i(2,:)';
            x_return_voltage_r = joystick_hold_i(3,:)';
            y_return_voltage_r = joystick_hold_i(4,:)';
            x_return = 0.5-max_dist*(x_return_voltage-0.5)./0.5;
            y_return = 0.5+max_dist*(y_return_voltage-0.5)./0.5;
            x_return_r = 0.5-max_dist*(x_return_voltage_r-0.5)./0.5;
            y_return_r = 0.5+max_dist*(y_return_voltage_r-0.5)./0.5;

            % x_center_csv0 = joystick_center_hold_csv_i(1,:)./(csv_screensize(1));
            % y_center_csv0 = joystick_center_hold_csv_i(2,:)./(csv_screensize(2));
            % x_return_csv0 = joystick_hold_csv_i(1,:)./(csv_screensize(1));
            % y_return_csv0 = joystick_hold_csv_i(2,:)./(csv_screensize(2));
            
            % x_center_csv = TrajectoryInterp(size(x_center,1),x_center_csv0);
            % y_center_csv = 1 - TrajectoryInterp(size(y_center,1),y_center_csv0);
            % x_return_csv = TrajectoryInterp(size(x_return,1),x_return_csv0);
            % y_return_csv = 1 - TrajectoryInterp(size(y_return,1),y_return_csv0);            

            if B*C==0  % Unimanual
                curr_x_center = x_center_r - 0.24; % x_center_r LeftX
                curr_y_center = y_center_r;        % y_center_r LeftY   
                curr_x_return = x_return_r - 0.24; 
                curr_y_return = y_return_r;        
                curr_x_trials = x_trials_r - 0.24;
                curr_y_trials = y_trials_r;
                x_center_r = x_center + 0.24;      % x_center RightX
                y_center_r = y_center;             % y_center RightY
                x_return_r = x_return + 0.24;     
                y_return_r = y_return;
                x_trials_r = x_trials + 0.24;
                y_trials_r = y_trials;
                x_center = curr_x_center;    
                y_center = curr_y_center;
                x_return = curr_x_return;
                y_return = curr_y_return;
                x_trials = curr_x_trials;
                y_trials = curr_y_trials;
            else
                x_center = x_center - 0.24;    
                x_return = x_return - 0.24;
                x_trials = x_trials - 0.24;
                x_center_r = x_center_r + 0.24;
                x_return_r = x_return_r + 0.24;
                x_trials_r = x_trials_r + 0.24;
            end
            dist_x_square = (x_center-target_pos(B+1,1)).^2;
            dist_y_square = (y_center-target_pos(B+1,2)).^2;           
            dist_x_square_r = (x_center_r-target_pos_r(C+1,1)).^2;
            dist_y_square_r = (y_center_r-target_pos_r(C+1,2)).^2;              
			
			if B>0 & C==0     % Left
				judge_onset = [diff(x_center),diff(y_center)];
				[min_dist,min_idx] = min(dist_x_square+dist_y_square);
            else
                if B==0 & C>0 % Right
                    judge_onset = [diff(x_center_r),diff(y_center_r)];
                    [min_dist,min_idx] = min(dist_x_square_r+dist_y_square_r);
                else         % Bimanual
                    judge_onset = [diff(x_center),diff(y_center),diff(y_center_r),diff(y_center_r)];
                    [min_dist,min_idx] = min(dist_x_square+dist_y_square+dist_x_square_r+dist_y_square_r);
                end
            end
            
            if min_dist>100 | size(x_return,1)==0 | sum(phase_no_i==3)<2
            %if abs([x_center(end),y_center(end)] - target_pos(target_i-100,:))>0.25
                continue
            else
                trial_no_list = [trial_no_list;trial_no];
                count = count + 1;

                success_hold_center_flag = [success_hold_center_flag;Data.trial_info.flag.success_hold_center_flag(trial_id)];
                success_hold_flag = [success_hold_flag;Data.trial_info.flag.success_hold_flag(trial_id)];
                repeat_timeout_flag = [repeat_timeout_flag;Data.trial_info.flag.repeat_timeout_flag(trial_id)];
                repeat_flag = [repeat_flag;Data.trial_info.flag.repeat_flag(trial_id)];
                auto_flag = [auto_flag;Data.trial_info.flag.auto_flag(trial_id)];
                flag_i = struct('success_hold_center_flag_i',Data.trial_info.flag.success_hold_center_flag(trial_id),'success_hold_flag_i',Data.trial_info.flag.success_hold_flag(trial_id),'repeat_timeout_flag_i',Data.trial_info.flag.repeat_timeout_flag(trial_id),'repeat_flag_i',Data.trial_info.flag.repeat_flag(trial_id),'auto_flag_i',Data.trial_info.flag.auto_flag(trial_id));
                
                onset_bin=0;
                move_index = find(sum(abs(judge_onset'))>0.002);
                dot_list =[];
				if B*C==0
					for ind = 1 : size(move_index,2)
						dot_list = [dot_list;dot(judge_onset(move_index(ind),:),target_pos0(B+C+1,:))];
                    end 
                else
					for ind = 1 : size(move_index,2)
						dot_list = [dot_list;dot(judge_onset(move_index(ind),1:2),target_pos0(B+1,:))+dot(judge_onset(move_index(ind),3:4),target_pos0(C+1,:))];
                    end
                end
                
                for ind = 1 : size(move_index,2)-5
                    if move_index(ind+3) - move_index(ind) == 3 & abs(x_center(move_index(ind))-0.26)<0.12 & abs(y_center(move_index(ind))-0.5)<0.12 & abs(x_center_r(move_index(ind))-0.74)<0.12 & abs(y_center_r(move_index(ind))-0.5)<0.12 & dot_list(ind)>0
                        onset_bin = move_index(ind);
                        onset_bin = onset_bin+1;
                        break
                    end
                end
                

                if onset_bin==0
                    onset_bin = size(x_center,1);
                end
                %move_index_diff = find(diff(move_index)==1);
                %onset_bin = move_index(move_index_diff)+1;
                %onset_bin = min(find(sum(abs(judge_onset'))>0.002))+1; % 'onset_bin' for center phase
                %onset_bin = onset_bin + min(find(phase_no_i==3)) - 1;  % 'onset_bin' for a trial
                %Data.spike.spikes_50Hz = ones([3,200000]);
                movement_onset_mask(count,onset_bin:max(find(phase_no_i==3))) = 1;
                move_id = movement_onset_mask(count,:)>0; % 1: go phase & Move
                time_center = sum(movement_onset_mask(count,:)>0)/meta_trial_i.f0;
                if Data.trial_info.flag.success_hold_center_flag(trial_id) ==1 % path efficiency应该为倒数
                    if C==0
						path_eff0 = sum(sqrt(sum(judge_onset(onset_bin:end,:)'.^2)))./sqrt((x_center(min_idx)-0.26)^2+(y_center(min_idx)-0.5)^2);
                    elseif B==0
						path_eff0 = sum(sqrt(sum(judge_onset(onset_bin:end,:)'.^2)))./sqrt((x_center_r(min_idx)-0.74)^2+(y_center_r(min_idx)-0.5)^2);
                    else
                        path_eff0 = (sum(sqrt(sum(judge_onset(onset_bin:end,1:2)'.^2)))+sum(sqrt(sum(judge_onset(onset_bin:end,3:4)'.^2))))./(sqrt((x_center(min_idx)-0.26)^2+(y_center(min_idx)-0.5)^2)+sqrt((x_center_r(min_idx)-0.74)^2+(y_center_r(min_idx)-0.5)^2));
                    end
                else
					if C==0
						path_eff0 = sum(sqrt(sum(judge_onset'.^2)))./ sqrt((0.26-target_pos(B+1,1)).^2 + (0.5-target_pos(B+1,2)).^2);
                    elseif B==0
						path_eff0 = sum(sqrt(sum(judge_onset'.^2)))./ sqrt((0.74-target_pos_r(C+1,1)).^2 + (0.5-target_pos_r(C+1,2)).^2);
                    else
						path_eff0 = (sum(sqrt(sum(judge_onset(:,1:2)'.^2)))+sum(sqrt(sum(judge_onset(:,3:4)'.^2))))./ (sqrt((0.26-target_pos(B+1,1)).^2 + (0.5-target_pos(B+1,2)).^2)+sqrt((0.74-target_pos_r(C+1,1)).^2 + (0.5-target_pos_r(C+1,2)).^2));
                    end
                end
                %path_eff = sqrt((0.5-target_pos(target_i-100,1)).^2 + (0.5-target_pos(target_i-100,2)).^2)./sum(abs(judge_onset));
                
                path_eff = [path_eff;path_eff0];
                % path_eff larger, trajectory worse
                spikes_50Hz_i = Data.spike.spikes_50Hz(:,find(Data.trial_info.trial_no == trial_no)); 
                %waveforms_50Hz_i = Data.spike.waveforms_50Hz(:,find(Data.trial_info.trial_no == trial_no)); 
				
                if success_id(i) == 1
                    fig = figure('Visible', 'off');
                    % fig = figure(1);
                    set(fig, 'Position', [1,1,480,480]);
                    %xlabel = [0:30:360];
                    %xtick = pi*xlabel/180;
					if B>0
						target_x = [target_pos(B+1,1)-0.08;target_pos(B+1,1)-0.08;target_pos(B+1,1)+0.08;target_pos(B+1,1)+0.08;target_pos(B+1,1)-0.085];
						target_y = [target_pos(B+1,2)-0.08;target_pos(B+1,2)+0.08;target_pos(B+1,2)+0.08;target_pos(B+1,2)-0.08;target_pos(B+1,2)-0.085];
						plot(target_x,target_y,'b-','Markersize',20,'LineWidth',8);
						hold on;
						plot(0.26,0.5,'go','Markersize',30,'LineWidth',10);
						hold on;			
						plot(x_center(1:2:onset_bin),y_center(1:2:onset_bin),'om','Markersize',10,'LineWidth',3);
						hold on;
						plot(x_center(onset_bin:end),y_center(onset_bin:end),'ok','Markersize',15,'LineWidth',3);
						hold on;
						min_dist_id = max(find(dist_x_square + dist_y_square == min(dist_x_square + dist_y_square)));% 计算最近的位置，x方+y方，求最小值
						plot(x_center(min_dist_id),y_center(min_dist_id),'bo','Markersize',30,'LineWidth',5);
						hold on;		
						plot(x_return(1:2:end),y_return(1:2:end),'--r','Markersize',25,'LineWidth',5);
						hold on;
						plot(x_return(end),y_return(end),'ro','Markersize',30,'LineWidth',5);
						hold on;
						plot(x_trials(1:2:end),y_trials(1:2:end),'g--','Markersize',20,'LineWidth',3);
                    end
                    if C>0
						target_xr = [target_pos_r(C+1,1)-0.08;target_pos_r(C+1,1)-0.08;target_pos_r(C+1,1)+0.08;target_pos_r(C+1,1)+0.08;target_pos_r(C+1,1)-0.085];
						target_yr = [target_pos_r(C+1,2)-0.08;target_pos_r(C+1,2)+0.08;target_pos_r(C+1,2)+0.08;target_pos_r(C+1,2)-0.08;target_pos_r(C+1,2)-0.085];
						plot(target_xr,target_yr,'b-','Markersize',20,'LineWidth',8);
						hold on;
						plot(0.74,0.5,'go','Markersize',30,'LineWidth',10);
						hold on;			
						plot(x_center_r(1:2:onset_bin),y_center_r(1:2:onset_bin),'om','Markersize',10,'LineWidth',3);
						hold on;	
						plot(x_center_r(onset_bin:end),y_center_r(onset_bin:end),'ok','Markersize',15,'LineWidth',3);
						hold on;
						min_dist_id = max(find(dist_x_square_r + dist_y_square_r == min(dist_x_square_r + dist_y_square_r)));% 计算最近的位置，x方+y方，求最小值
						plot(x_center_r(min_dist_id),y_center_r(min_dist_id),'bo','Markersize',30,'LineWidth',5);
						hold on;		
						plot(x_return_r(1:2:end),y_return_r(1:2:end),'--r','Markersize',25,'LineWidth',5);
						hold on;
						plot(x_return_r(end),y_return_r(end),'ro','Markersize',30,'LineWidth',5);
						hold on;
						plot(x_trials_r(1:2:end),y_trials_r(1:2:end),'g--','Markersize',20,'LineWidth',3);
                    end
                    hold on;

                    set(gca,'XLim',[0 1]);
                    set(gca,'YLim',[0 1]);    
                    fig_name = ['TrialNo.' num2str(trial_no) '; Target.' num2str(Data.trial_info.target_list(trial_id)) '.tif'];
                    title([ ' Trial.' num2str(trial_no) '; Tar.' num2str(Data.trial_info.target_list(trial_id)) ' ;Go' num2str(time_center) ' ;Time' num2str(size(x_trials,1)/meta_trial_i.f0) ' ;Cen' num2str(Data.trial_info.flag.success_hold_center_flag(trial_id)) ' ;Succ' num2str(Data.trial_info.flag.success_hold_flag(trial_id)) ';PEff' num2str(path_eff0)]);
                    saveas(fig, fullfile(save_trial_dir, fig_name) );
                    hold off;
                    clf
                end
                

                target_list = [target_list;target_i];
                mean_fr_list = [mean_fr_list;mean(spikes_50Hz_i')*meta_trial_i.f0];

                %mean_center_fr = Data.spike.spikes_50Hz(:,find(Data.trial_info.trial_no == trial_id & Data.trial_info.phase_no==3))';
                %mean_center_fr_list = [mean_center_fr_list;mean(mean_center_fr)*meta_trial_i.f0];
                if sum(move_id)>0
                    mean_center_fr = spikes_50Hz_i(:,move_id)';
                    if size([mean_center_fr],1) ==1
                        mean_center_fr_list = [mean_center_fr_list;[mean_center_fr]*meta_trial_i.f0];
                    else
                        mean_center_fr_list = [mean_center_fr_list;mean([mean_center_fr])*meta_trial_i.f0];
                    end
                else
                    mean_center_fr=0;
                    mean_center_fr_list = [mean_center_fr_list;zeros([1,size(spikes_50Hz_i,1)])];
                end
                time_center_list = [time_center_list;time_center];
                mean_peripheral_fr = Data.spike.spikes_50Hz(:,[find(Data.trial_info.trial_no == trial_no & Data.trial_info.phase_no==4),find(Data.trial_info.trial_no == trial_no & Data.trial_info.phase_no==6)])';
                mean_peripheral_fr_list = [mean_peripheral_fr_list;mean(mean_peripheral_fr)*meta_trial_i.f0];
                time_peripheral = size(mean_peripheral_fr,1)./meta_trial_i.f0;
                time_peripheral_list = [time_peripheral_list;time_peripheral];
                if size(find(Data.trial_info.trial_no == trial_no & Data.trial_info.phase_no==8),2)<9
                    baseline_id = find(Data.trial_info.trial_no == trial_no);
                    baseline_id = baseline_id(end-9:end);
                else
                    baseline_id = find(Data.trial_info.trial_no == trial_no & Data.trial_info.phase_no==8);
                end
                baseline = mean(Data.spike.spikes_50Hz(:,baseline_id)');
                baseline_sd = std(Data.spike.spikes_50Hz(:,baseline_id)');
                baseline_list = [baseline_list;baseline*meta_trial_i.f0];
                baseline_sd_list = [baseline_sd_list;baseline_sd*meta_trial_i.f0];

                %dist_x_square_list = [dist_x_square_list;min(dist_x_square)];
                %dist_y_square_list = [dist_y_square_list;min(dist_y_square)];
                %min_sqr_dist_list = [min_sqr_dist_list;min(dist_x_square+dist_y_square)];

                spike_i = struct('onset_bin',onset_bin,'readNeurons',Data.spike.readNeurons,'spikes_50Hz_i',spikes_50Hz_i);
                joystick_i = struct('sensor_i',sensor_i,'joystick_trials_i',joystick_trials_i,'joystick_hold_i',joystick_hold_i,'joystick_center_hold_i',joystick_center_hold_i,'go_last_time_i',Data.Analog.joystick.go_last_time(trial_id));
                trajectory_i = struct('traject_center',[x_center,y_center,x_center_r,y_center_r],'traject_peripheral',[x_return,y_return,x_return_r,y_return_r],'traject_trial',[x_trials,y_trials,x_trials_r,y_trials_r],'path_eff',path_eff0) % 'path_eff',path_eff
                %trajectory_i = struct('traject_center',[x_center,y_center],'traject_peripheral',[x_return,y_return],'traject_center_csv',[x_center_csv,y_center_csv],'traject_peripheral_csv',[x_return_csv,y_return_csv],'traject_trial',[x_trials,y_trials],'dist_x_square',dist_x_square,'dist_y_square',dist_y_square,'min_sqr_dist',min(dist_x_square+dist_y_square),'path_eff',path_eff0) % 'path_eff',path_eff

                FR_i = struct('mean_fr',mean(spikes_50Hz_i')*meta_trial_i.f0,'mean_center_fr',mean([mean_center_fr])*meta_trial_i.f0,'mean_peripheral_fr',mean(mean_peripheral_fr)*meta_trial_i.f0,'mean_baseline_fr',baseline*meta_trial_i.f0,'mean_baseline_sd_fr',baseline_sd*meta_trial_i.f0);
                time_i = struct('time_center',time_center,'time_peripheral',time_peripheral);
                
                trial_info_i = struct('trial_list',trial_no,'target_list_i',target_i,'phase_no_i',phase_no_i,'flag_i',flag_i);
                Data_triali = struct('meta_trial_i',meta_trial_i,'trial_info_i',trial_info_i,'joystick_i',joystick_i,'spike_i',spike_i,'FR_i',FR_i,'time_i',time_i,'trajectory_i',trajectory_i);
                save([save_trial_dir date '.TrialNo.' num2str(trial_no) '.TargetNo.' num2str(Data.trial_info.target_list(trial_id)) '.mat'],'Data_triali');

                ['Save TrialNo' num2str(trial_no)]
                
            end
        end
    end
    movement_onset_mask = movement_onset_mask(1:count,:);
    flag = struct('success_hold_center_flag',success_hold_center_flag,'success_hold_flag',success_hold_flag,'repeat_timeout_flag',repeat_timeout_flag,'repeat_flag',repeat_flag,'auto_flag',auto_flag);
    save([save_trial_dir date '.DatafileID' datafile_id '.Session.' num2str(session) '.mat'],'trial_no_list','movement_onset_mask','target_list','time_center_list','time_peripheral_list','mean_fr_list','mean_center_fr_list','mean_peripheral_fr_list','baseline_list','baseline_sd_list','min_sqr_dist_list','path_eff','flag');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%   7 Save Trial Data   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
end