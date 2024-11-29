function [] = HistGo_interp(date,area,area_ch,savedir_raster)
    %load([data_dir_session date '.Session.' session '.mat']);
    %load([data_dir date '.S' num2str(session) '\' date '.Session.' session '.mat']);
    
    load([savedir_raster date '.DataRasterGo_ManualSelected.mat']);
    
    target_num = DataRasterGo_ManualSelected.TargetGo.target_num;
    target_num_curr = [0;target_num];
    target_unique = DataRasterGo_ManualSelected.TargetGo.target_unique;
    target_list = DataRasterGo_ManualSelected.TargetGo.target_list;

    mean_fr_session = DataRasterGo_ManualSelected.FiringRateGo.mean_fr_session;
    readNeurons = DataRasterGo_ManualSelected.RasterGo.NeuronNo;
    
    raster = DataRasterGo_ManualSelected.RasterGo.raster_go;
    time_map = DataRasterGo_ManualSelected.RasterGo.time_map_go;

    psth_binlen = DataRasterGo_ManualSelected.RasterGo.psth_binlen;
    before_go_binlen = DataRasterGo_ManualSelected.RasterGo.before_go_binlen;

    %% Hist Go Alignment

    %idlen = [0;];
    %id = [];
    %for p = 1:5
    %    time_map_curr = time_map(:,:,p);
    %    idlen = [idlen;max(find(sum(time_map_curr)>0))-min(find(sum(time_map_curr)>0))+1];
    %    id = [id;min(find(sum(time_map_curr)>0));max(find(sum(time_map_curr)>0))];
    %end

    

    %for n = 1 : size(new_raster,1)
    %    for t = 1 : size(target_unique,1)
    %        spike_curr = new_raster(n,(sum(target_num_curr(1:t))+1):sum(target_num_curr(1:t+1)),:);
    %        for p = 1:5
    %            time_map_curr = new_time_map((sum(target_num_curr(1:t))+1):sum(target_num_curr(1:t+1)),:,p);
    %            spike_curr0 = reshape(spike_curr,size(spike_curr,[2,3])).* time_map_curr;
    %            spike_curr0 = sum(spike_curr0)';
    %            hist_go(n,t,(1+sum(idlen(1:p))):sum(idlen(1:(p+1)))) =  spike_curr0(id(p*2-1):id(p*2))./ (sum(time_map_curr(:,id(p*2-1):id(p*2)))'+0.00001);
    %        end
    %    end
    %end
    hist_go = zeros([size(raster,1),size(target_unique,1),psth_binlen]); %NeuronNum*8*binlen
    for n = 1 : size(raster,1)
        for t = 1 : size(target_unique,1)
            spike_curr = raster(n,(sum(target_num_curr(1:t))+1):sum(target_num_curr(1:t+1)),:);
            time_map_curr = sum(time_map((sum(target_num_curr(1:t))+1):sum(target_num_curr(1:t+1)),:,:),3);
            spike_curr0 = reshape(spike_curr,[size(spike_curr,2),size(spike_curr,3)]).* time_map_curr;
            hist_go(n,t,:) =  sum(spike_curr0,1)./ (sum(time_map_curr,1)+0.00001);
        end
    end

    [' ======= Start draw Hist -- Go Phase ======= ']

    for n = 1 : size(raster,1)
        clf;
        fig = figure('Visible', 'off');
        set(fig, 'Position', [1,1,960,960]);
        value=hist_go(n,:,:);
        value = reshape(value,[size(hist_go,2),size(hist_go,3)]);

        %plot(value','DisplayName','hist_go','LineWidth',2)
        %smoothdata(value,2,'gaussian',15)
        value = smoothdata(value,2,'gaussian',15);
        ylim = max([max(max(value)),0.6]);
        plot(value','DisplayName','hist_go','LineWidth',2)
        hold on;

        plot([before_go_binlen,before_go_binlen],[0,ylim],'--r','Markersize',3,'LineWidth',2);
        set(gca,'xtick',[before_go_binlen]);
        set(gca,'xticklabel',['Move']);
        %set(gca,'xticklabel',['Prep';' GO ';'Back';'Rewa';'Rest']);
        set(gca,'YLim',[0 ylim]);    
%         legend_list = ['101';'102';'103';'104';'105';'106';'107';'108'];
%         legend_curr=[];
%         for i = 1:size(target_unique,1)
%             legend_curr = [legend_curr;legend_list((target_unique(i)-100),:)];
%         end
        legend_curr=[];
        for i = 1:size(target_unique,1)
            legend_curr = [legend_curr;num2str(target_unique(i))];
        end
            
        legend(legend_curr);
        title([date  ' CH' num2str(ceil(double(readNeurons(n))/6))  ' NeuronNo.' num2str(readNeurons(n))]);
        fig_name4 = [date  ' CH' num2str(ceil(double(readNeurons(n))/6)) ' NeuronNo.' num2str(readNeurons(n)) ' HistGo' '.tif']
        saveas(fig, fullfile(savedir_raster, fig_name4) );

    end
    [' ======= End draw Hist -- Go Phase ======= ']
    %% Hist Interp
    
    hist_interp = zeros([size(raster,1),size(target_unique,1),sum(DataRasterGo_ManualSelected.InterpData.interp_num0)]); %NeuronNum*8*binlen
    for n = 1 : size(raster,1)
        for t = 1 : size(target_unique,1)
            spike_curr = DataRasterGo_ManualSelected.InterpData.interp_raster(n,(sum(target_num_curr(1:t))+1):sum(target_num_curr(1:t+1)),:);
            spike_curr0 = reshape(spike_curr,[size(spike_curr,2),size(spike_curr,3)]);
            hist_interp(n,t,:) =  mean(spike_curr0,1);
        end
    end
    save([savedir_raster date '.HistInterp.mat'],'hist_interp');
    
    [' ======= Start draw Hist -- Interp ======= ']

    for n = 1 : size(raster,1)
        clf;
        fig = figure('Visible', 'off');
        set(fig, 'Position', [1,1,960,960]);
        value=hist_interp(n,:,:);
        value = reshape(value,[size(hist_interp,2),size(hist_interp,3)]);

        %plot(value','DisplayName','hist_go','LineWidth',2)
        %smoothdata(value,2,'gaussian',15)
        value = smoothdata(value,2,'gaussian',15);
        ylim = max([max(max(value)),0.6]);
        plot(value','DisplayName','hist_interp','LineWidth',2)
        hold on;
        
        for i = 1:size(DataRasterGo_ManualSelected.InterpData.interp_num0,2)-2
            plot([sum(DataRasterGo_ManualSelected.InterpData.interp_num0(1:i+1)),sum(DataRasterGo_ManualSelected.InterpData.interp_num0(1:i+1))],[0,ylim],'--r','Markersize',3,'LineWidth',2);
        end
        xtick = [];
        for i = 1:size(DataRasterGo_ManualSelected.InterpData.interp_num0,2)-2
            xtick = [xtick;sum(DataRasterGo_ManualSelected.InterpData.interp_num0(1:i+1))];
        end
        set(gca,'xtick',xtick);
        xticklabel = [' MoveStart ';'ReturnStart';'RewardStart';' RestStart '];
        set(gca,'xticklabel',xticklabel);
        set(gca,'XLim',[0 sum(DataRasterGo_ManualSelected.InterpData.interp_num0)]);  
        set(gca,'YLim',[0 ylim]);    
%         legend_list = ['101';'102';'103';'104';'105';'106';'107';'108'];
%         legend_curr=[];
%         for i = 1:size(target_unique,1)
%             legend_curr = [legend_curr;legend_list((target_unique(i)-100),:)];
%         end
        legend_curr=[];
        for i = 1:size(target_unique,1)
            legend_curr = [legend_curr;num2str(target_unique(i))];
        end
            
        legend(legend_curr);
        title([date  ' CH' num2str(ceil(double(readNeurons(n))/6))  ' NeuronNo.' num2str(readNeurons(n))]);
        fig_name4 = [date  ' CH' num2str(ceil(double(readNeurons(n))/6)) ' NeuronNo.' num2str(readNeurons(n)) ' Hist' '.Interp.tif']
        saveas(fig, fullfile(savedir_raster, fig_name4) );

    end
    [' ======= End draw Hist -- Interp ======= ']
end
