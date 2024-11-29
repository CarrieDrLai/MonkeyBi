function [] = PSTHRasterGo(date,session,area,area_ch,data_dir,save_autoselectedtrial_dir,savedir_raster,threshFR)

    %%% Behavioral Plot

    %%% Mean FR Plot
    %load([data_dir_session date '.Session.' session '.mat']);
    %load([data_dir date '.S' num2str(session) '\' date '.Session.' session '.mat']);
    trial_dir = dir(data_dir);
    
    %% 3.1 Raster Each Neuron: Trial * Time; Align on Return Phase
    load([savedir_raster date '.DataRasterGo_ManualSelected.mat']);
    target_num = DataRasterGo_ManualSelected.TargetGo.target_num;
    target_unique = DataRasterGo_ManualSelected.TargetGo.target_unique;
    psth_binlen = DataRasterGo_ManualSelected.RasterGo.psth_binlen;
    before_return_binlen = 200;
    raster_return = zeros([size(DataRasterGo_ManualSelected.RasterGo.raster_go,1),size(DataRasterGo_ManualSelected.RasterGo.raster_go,2),psth_binlen]);
    trajectory_map_return = zeros([size(DataRasterGo_ManualSelected.RasterGo.raster_go,2),psth_binlen,4]);
    time_map_return = zeros([size(DataRasterGo_ManualSelected.RasterGo.raster_go,2),psth_binlen,5]);
    
    load([data_dir trial_dir(selected_trial_mat_sort(1)).name]);
    [n,nid,nidb]= intersect(DataRasterGo_ManualSelected.RasterGo.NeuronNo,Data_triali.spike_i.readNeurons);
    
    
    for i = 1:sum(target_num) 
        load([data_dir date '.TrialNo.' num2str(DataRasterGo_ManualSelected.TargetGo.trial_no(i)) '.TargetNo.' num2str(DataRasterGo_ManualSelected.TargetGo.target_list(i)) '.mat']);
        spike_i = Data_triali.spike_i.spikes_50Hz_i(nidb,:);
        trajectory_i = Data_triali.trajectory_i.traject_trial;
        before_return0 = min([find(Data_triali.trial_info_i.phase_no_i==4),find(Data_triali.trial_info_i.phase_no_i==6)]) - 1; % before Return#-1
        
        start_point = max(before_return_binlen + 1 - before_return0,1);    
        end_point = min(size(Data_triali.trial_info_i.phase_no_i,2)-before_return0 + before_return_binlen,psth_binlen);
        start_point_spike = max(before_return0 - before_return_binlen + 1,1);
        end_point_spike = start_point_spike + end_point-start_point;
        %[before_return0,start_point,end_point,start_point_spike,end_point_spike]
        raster_return(:,i,start_point:end_point) = spike_i(:,start_point_spike:end_point_spike);
        trajectory_map_return(i,start_point:end_point,:) = trajectory_i(start_point_spike:end_point_spike,:);
        return_len = min(sum(DataRasterGo_ManualSelected.RasterGo.time_map_go(i,:,3)),100);
        go_len = min(sum(DataRasterGo_ManualSelected.RasterGo.time_map_go(i,:,2)),before_return_binlen-1);
        
        time_map_return(i,start_point:(before_return_binlen-go_len),1) = 1; %beforego
        time_map_return(i,(before_return_binlen-go_len+1):before_return_binlen,2) = 1; %go
        time_map_return(i,(before_return_binlen+1):(before_return_binlen+return_len),3) = 1; % return
        time_map_return(i,(before_return_binlen+return_len+1):(before_return_binlen+return_len+sum(sum(DataRasterGo_ManualSelected.RasterGo.time_map_go(i,:,4)))),4) = 1; 
        time_map_return(i,(before_return_binlen+return_len+sum(sum(DataRasterGo_ManualSelected.RasterGo.time_map_go(i,:,4)))+1):end_point,5) = 1; 
    end
    %move_mark = [move_mark;size(find(Data_triali.trial_info_i.phase_no_i==3),2)-Data_triali.spike_i.onset_bin+1];
    %time_map_return = [before_return,time_map(:,2:5)];
    go_len=[];
    for i = 1:sum(target_num) 
        load([data_dir trial_dir(selected_trial_mat_sort(i)).name]);
        go_len = [go_len;min(sum(DataRasterGo_ManualSelected.RasterGo.time_map_go(i,:,2)),before_return_binlen-1)];
    end
    
    %new_raster_return = raster_return(:,target_id,:);
    %new_time_map_return = time_map_return(target_id,:);
    RasterReturn = struct('raster_return',raster_return,'before_return_binlen',before_return_binlen,'psth_binlen',psth_binlen,'time_map_return',time_map_return,'NeuronNo',DataRasterGo_ManualSelected.RasterGo.NeuronNo);
    TargetReturn = struct('target_num',target_num,'target_list',DataRasterGo_ManualSelected.TargetGo.target_list,'target_unique',target_unique,'trial_no',DataRasterGo_ManualSelected.TargetGo.trial_no);
    FiringRateReturn = struct('mean_fr_session',DataRasterGo_ManualSelected.FiringRateGo.mean_fr_session);
    TrajectoryReturn = struct('trajectory_map',trajectory_map_return);
    DataRasterReturn = struct('RasterReturn',RasterReturn,'TargetReturn',TargetReturn,'FiringRateReturn',FiringRateReturn,'TrajectoryReturn',TrajectoryReturn);

    save([savedir_raster date '.DataRasterReturn.mat'],'DataRasterReturn');

     %% 3.2 draw new raster Return
    load([savedir_raster date '.DataRasterReturn.mat'],'DataRasterReturn');
    psth_binlen = DataRasterReturn.RasterReturn.psth_binlen;
    before_return_binlen = DataRasterReturn.RasterReturn.before_return_binlen;
    target_num = DataRasterReturn.TargetReturn.target_num;
    target_unique = DataRasterReturn.TargetReturn.target_unique;    
    time_map_return = DataRasterReturn.RasterReturn.time_map_return;
    raster_return = DataRasterReturn.RasterReturn.raster_return;
    
    [' ======= Start draw Raster -- Return Phase ======= ']
    for n = 1 : size(raster_return,1)
        clf;    
        fig = figure('Visible', 'off');
        set(fig, 'Position', [1,1,960,960]);
        %
        %fig = figure(1);
        raster_curr = raster_return(n,:,:);
        raster_curr = reshape(raster_curr,[size(raster_return,2),size(raster_return,3)]);
        imagesc(raster_curr>0);
        hold on;

        plot([before_return_binlen+0.5,before_return_binlen+0.5],[0,psth_binlen],'-r','LineWidth',3);
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

        a=[];
        for i = 1:sum(target_num) 
            a = [a,size(find(time_map_return(i,:,1)>0),2)];
        end
        
        ytick = [ytick;(sum(target_num(1:end-1))+sum(target_num))*0.5];
        for i = 1:sum(target_num)
            if size(find(time_map_return(i,:,1)>0),2)==0
                plot(min(find(time_map_return(i,:,2)>0)),i,'go','Markersize',3,'LineWidth',2);
                hold on;
            else
                plot(min(find(time_map_return(i,:,1)>0)),i,'go','Markersize',3,'LineWidth',2);
                hold on;                
            end

            plot(min(find(time_map_return(i,:,2)>0)),i,'ro','Markersize',5,'LineWidth',2);
            hold on;
            plot(max(find(time_map_return(i,:,3)>0)),i,'go','Markersize',5,'LineWidth',2);
            hold on;
            if size(max(find(time_map_return(i,:,4)>0)),2)==0
                plot(max(find(time_map_return(i,:,3)>0))+1,i,'ro','Markersize',3,'LineWidth',2);
            else
                plot(max(find(time_map_return(i,:,4)>0)),i,'ro','Markersize',3,'LineWidth',2);
            end
            hold on;
            if size(max(find(time_map_return(i,:,5)>0)),2)==0
                if size(max(find(time_map_return(i,:,4)>0)),2)>0
                    plot(max(find(time_map_return(i,:,4)>0))+1,i,'go','Markersize',3,'LineWidth',2);
                else
                    plot(max(find(time_map_return(i,:,3)>0))+2,i,'go','Markersize',3,'LineWidth',2);
                end
            else
                plot(max(find(time_map_return(i,:,5)>0)),i,'go','Markersize',3,'LineWidth',2);
            end
        end

        title([date ' CH' num2str(ceil(double(DataRasterReturn.RasterReturn.NeuronNo(n))/6)) ' NeuronNo.' num2str(DataRasterReturn.RasterReturn.NeuronNo(n)) '.Return']);
        xlabel = [before_return_binlen];
        xtick = 'on';
        set(gca,'xtick',before_return_binlen);
        set(gca,'xticklabel','Return');
        set(gca,'ytick',ytick);
        yticklabel_list = [' T';'TR';' R';'BR';' B';'BL';' L';'TL'];
		if target_unique(1)>400
			yticklabel_list = ['LT  ';'L TR';'LR  ';'L BR';'LB  ';'  RT';'R TL';'  RL';'R BL';'  RB';'T  T';'TRTL';'LRRL';'BRBL';'B  B';];
        end
        yticklabel=[];
        for i = 1:size(target_unique,1)
            yticklabel = [yticklabel;yticklabel_list((i),:)];
        end
        set(gca,'yticklabel',yticklabel);
        hold off;
        fig_name2 = [date ' CH' num2str(ceil(double(DataRasterReturn.RasterReturn.NeuronNo(n))/6)) ' NeuronNo.' num2str(DataRasterReturn.RasterReturn.NeuronNo(n)) '.Return.tif']
        saveas(fig, fullfile(savedir_raster, fig_name2) );
    end
    [' ======= End draw Raster -- Return Phase ======= ']
end