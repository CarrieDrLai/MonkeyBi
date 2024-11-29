function [Data] = InitializationAnalogInputMain_9Grid(NEV,readNeurons,nsp_ch,raw_data,date,data_dir,data_dir_joystick,session,session_csv,savedir,f,f_joystick,save_path)

    %%% 20220803- %%%%
    %%% 128CH left M1 +PMd
    %%% 256CH left M1 + right + PMd

    %%% AnalogInput 完成信号alignment代码 - sorting - 数据保存 - 行为学筛选

    %%% joystick with sensor
    %%% 1. joystick prepare & prepare: 先观察，考虑去除 
    %%% 2. 判断条件：go阶段(a)起始位置 (b) 完成时长 (c)轨迹


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Config %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % nsp_ch = 128
    % trajectory = NS6.Data(nsp_ch+10:nsp_ch+11,:)./32768
    % trajectory = trajectory';
    %fig = figure(1);
    %plot(trajectory);

    % savedir = data_dir;
    if ~exist(savedir, 'dir')
        mkdir(savedir);
    end

    %f = 30000;
    %f_joystick = 100;
    f0 = 50;
    start_prepare_id = 1;%一般=1 %可能有时候开了nsprecording，关闭ui,重新打开ui,会造成重新计数的情况

    [trial_list,trial_no,phase_no,spikes_50Hz,start_point,start_point_30k] = InitializationcompileNEVNS6file(NEV,readNeurons,nsp_ch,start_prepare_id,f,f0,f_joystick);
    ['compile NEV file 1 -- done']
    %[trial_list2,trial_no2,phase_no2,spikes_50Hz2,readNeurons2,nsp_ch2,start_point2] = compileNEVfile(NEV2,start_prepare_id,f0,f_joystick);
    %['compile NEV file 2 -- done']


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Behavior Task
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%  4  Joystick   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [sensor1,trajectory,stable,trajectory2] = InitializationAnalogInputLabel(raw_data);
    clear raw_data;
    trajectory = DownSampleRawData(trajectory,trial_no,start_point,f,f0,f_joystick);
    trajectory2 = DownSampleRawData(trajectory2,trial_no,start_point,f,f0,f_joystick);
    sensor_compress = [sensor1;stable];
    sensor_compress = DownSampleRawData(sensor_compress,trial_no,start_point,f,f0,f_joystick);
%     traject_start_point = uint32(f_joystick*(double(start_point)/double(f)));
%     temp = trajectory(:,traject_start_point:traject_start_point+uint32(double(size(trial_no,2)-1)/f0*f_joystick-1));
%     temp = reshape(temp,[2,uint32(f_joystick)/uint32(f0),size(temp,2)/(uint32(f_joystick)/uint32(f0))]);
%     trajectory = mean(temp,2);
%     trajectory = reshape(trajectory,[2,size(trial_no,2)-1]);
%     temp = trajectory2(:,traject_start_point:traject_start_point+uint32(double(size(trial_no,2)-1)/f0*f_joystick-1));
%     temp = reshape(temp,[2,uint32(f_joystick)/uint32(f0),size(temp,2)/(uint32(f_joystick)/uint32(f0))]);
%     trajectory2 = mean(temp,2);
%     trajectory2 = reshape(trajectory2,[2,size(trial_no,2)-1]);
    
    joystick_center_hold = ones(4*(size(trial_list,1)-1),12000)*(-1);
    joystick_hold = ones(4*(size(trial_list,1)-1),5000)*(-1);
    joystick_trials = ones(4*(size(trial_list,1)-1),20000)*(-1);
    joystick_center_hold_csv = ones(4*(size(trial_list,1)-1),12000)*(-1);
    joystick_hold_csv = ones(4*(size(trial_list,1)-1),5000)*(-1);
    sensor_center_hold = ones(2*(size(trial_list,1)-1),12000)*(-1);
    sensor_hold = ones(2*(size(trial_list,1)-1),5000)*(-1);
    sensor_trials = ones(2*(size(trial_list,1)-1),20000)*(-1);
    % phase_hold_center,3; phase_hold, 6
    for i = 1:(size(trial_list,2)-1)
        trial_i = trial_list(i);
        trajectory_trial = find(trial_no == trial_i);
        peripheral = [find(trial_no == trial_i & phase_no ==4),find(trial_no == trial_i & phase_no ==6)];
        center = find(trial_no == trial_i & phase_no>0);
        if size(peripheral,2)>0
            center = center(center<peripheral(1));
            trajectory_temp = [trajectory(:,peripheral(1):peripheral(end))];
            trajectory2_temp = [trajectory2(:,peripheral(1):peripheral(end))];
            joystick_hold((4*i-3):(4*i-2),1:size(trajectory_temp,2))=trajectory_temp;      
            joystick_hold((4*i-1):(4*i),1:size(trajectory_temp,2))=trajectory2_temp;      
            sensor_hold((2*i-1):(2*i),1:size(trajectory_temp,2))=sensor_compress(:,peripheral(1):peripheral(end));  
            fn=readtable([data_dir_joystick date '.Session-' session '.Trial-' num2str(trial_i) '.Trajectory.csv']);
            if size(fn,2)<3
                trajectory_temp = [[960;540],[960;540]];
                trajectory2_temp = [[960;540],[960;540]];
            else
                trajectory_temp = csvread([data_dir_joystick date '.Session-' session '.Trial-' num2str(trial_i) '.Trajectory.csv'],1,1);
                trajectory2_temp = trajectory_temp(3:4,:);
                trajectory_temp = trajectory_temp(1:2,:);
            end
            joystick_hold_csv((4*i-3):(4*i-2),1:size(trajectory_temp,2))=trajectory_temp; 
            joystick_hold_csv((4*i-1):(4*i),1:size(trajectory_temp,2))=trajectory2_temp; 
        end
        if size(center,2)>0
            trajectory_center = [trajectory(:,center(1):center(end))];
            trajectory2_center = [trajectory2(:,center(1):center(end))];
            joystick_center_hold((4*i-3):(4*i-2),1:size(trajectory_center,2))=trajectory_center;
            joystick_center_hold((4*i-1):(4*i),1:size(trajectory_center,2))=trajectory2_center;
            sensor_center_hold((2*i-1):(2*i),1:size(trajectory_center,2))=sensor_compress(:,center(1):center(end));  
            trajectory_center_temp = csvread([data_dir_joystick date '.Session-' session '.Trial-' num2str(trial_i) '.TrajectoryHOLD.csv'],1,1);
            joystick_center_hold_csv((4*i-3):(4*i),1:size(trajectory_center_temp,2))=trajectory_center_temp;
        end
        if size(trajectory_trial,2)>0
            trajectory_temp = [trajectory(:,trajectory_trial(1):min(size(trajectory,2),trajectory_trial(end)))];
            trajectory2_temp = [trajectory2(:,trajectory_trial(1):min(size(trajectory,2),trajectory_trial(end)))];
            joystick_trials((4*i-3):(4*i-2),1:size(trajectory_temp,2))=trajectory_temp;
            joystick_trials((4*i-1):(4*i),1:size(trajectory_temp,2))=trajectory2_temp;
            sensor_trials((2*i-1):(2*i),1:size(trajectory_temp,2))=sensor_compress(:,trajectory_trial(1):min(size(trajectory,2),trajectory_trial(end)));  
        end               
    end

    ['compile Joystick file -- done']

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%  4  Joystick   %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%   5 Target & Flag   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % session_csv = csvread([data_dir_joystick date '.Session-' session '.csv.csv'],1,1);

    target_list = session_csv(:,3);
    target_list = target_list(trial_list(1:end-1));
    success_hold_center_flag = session_csv(:,5);
    success_hold_center_flag = success_hold_center_flag(trial_list(1:end-1));
    success_hold_flag = session_csv(:,6);
    success_hold_flag = success_hold_flag(trial_list(1:end-1));
    repeat_timeout_flag = session_csv(:,7);
    repeat_timeout_flag = repeat_timeout_flag(trial_list(1:end-1));
    auto_flag = session_csv(:,8);
    auto_flag = auto_flag(trial_list(1:end-1));
    repeat_flag = session_csv(:,9);
    repeat_flag = repeat_flag(trial_list(1:end-1));

    go_last_time = session_csv(:,11);
    go_last_time = go_last_time(trial_list(1:end-1));
    ['compile CSV file -- done']
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%   5 Target & Flag   %%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%   6 Save Data   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    meta = struct('f',f,'f0',f0,'f_joystick',f_joystick,'date',date,'session',session);
    flag = struct('success_hold_center_flag',success_hold_center_flag,'success_hold_flag',success_hold_flag,'repeat_timeout_flag',repeat_timeout_flag,'auto_flag',auto_flag,'repeat_flag',repeat_flag);
    trial_info = struct('trial_list',trial_list(1:end-1),'target_list',target_list,'trial_no',trial_no,'phase_no',phase_no,'flag',flag)
    %trial_info = struct('trial_list',trial_list(1:end-1),'target_list',target_list,'trial_no',trial_no,'phase_no',phase_no,'trial_no2',trial_no2,'phase_no2',phase_no2,'flag',flag)
    spike = struct('readNeurons',readNeurons,'spikes_50Hz',spikes_50Hz);
    %spike2 = struct('readNeurons2',readNeurons2,'spikes_50Hz2',spikes_50Hz2);
    joystick = struct('joystick_trials',joystick_trials,'joystick_hold',joystick_hold,'joystick_center_hold',joystick_center_hold,'joystick_hold_csv',joystick_hold_csv,'joystick_center_hold_csv',joystick_center_hold_csv,'go_last_time',go_last_time);
    sensor = struct('sensor_center_hold',sensor_center_hold,'sensor_hold',sensor_hold,'sensor_trials',sensor_trials);
    Analog = struct('sensor',sensor,'joystick',joystick);
    Data = struct('meta',meta,'trial_info',trial_info,'Analog',Analog,'spike',spike);
    %Data = struct('meta',meta,'trial_info',trial_info,'Analog',Analog,'spike',spike,'spike2',spike2);
    save([save_path '.mat'],'Data');
    ['Save Session file -- done']

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%   6 Save Data   %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%    MUA    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %MUA = zeros(numOfTimestamps, 96);
    %for idx = 1:length(readElectrodes)
    %    MUA(timeStamps(NEV.Data.Spikes.Electrode == readElectrodes(idx)),readElectrodes(idx)) = 1;
    %end

    %MUA_50Hz = [];
    %t = 1;
    %while t*f/50 < size(MUA,1)
    %    MUA_50Hz = [MUA_50Hz;sum(MUA(((t-1)*(f/50)+1):t*f/50,:))];
    %    t=t+1;
    %end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%    MUA    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function timeStampsSec = timeStampsToTimestampSeconds(timeStamps, samplingRate)

        timeStampsSec = double(timeStamps) / double(samplingRate);

    end

end