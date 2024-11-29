function [] = MeanMap(date,area,area_ch,savedir_raster)
    %load([data_dir_session date '.Session.' session '.mat']);
    %load([data_dir date '.S' num2str(session) '\' date '.Session.' session '.mat']);
    
    load([savedir_raster date '.DataRasterReturn.mat']);
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

    %% 4.1 Sum & Mean for each phase
    % Sum map : NeuronNum * TrialNum * 5Phases
    % Mean map : NeuronNum * TargetNum * 5Phases
    sum_map = zeros([size(raster,1),size(raster,2),5]); 
    for n = 1 : size(raster,1)
        for p = 1:5
            sum_map(n,:,p) = sum((reshape(raster(n,:,:),[size(raster,2),size(raster,3)]).*time_map(:,:,p))')./(sum(time_map(:,:,p),2)+0.00001)';
        end
    end

    target_num_curr = [0;target_num];
    mean_map = zeros([size(raster,1),size(target_unique,1),5]); 
    for t = 1:size(target_unique,1)
        mean_map(:,t,:) = nanmean(sum_map(:,(sum(target_num_curr(1:t))+1):sum(target_num_curr(1:(t+1))),:),2);
    end
    %for n = 1 : size(new_raster,1)
    %    for t = 1:size(target_unique,1)
    %        mean_map(n,t,:) = nanmean(sum_map(n,(sum(target_num_curr(1:t))+1):sum(target_num_curr(1:(t+1))),:),2);
    %    end
    %end
    save([savedir_raster date '.RasterMean.mat'],'sum_map','mean_map','target_num','target_unique','target_list','mean_fr_session','readNeurons');

    %% 4.2 draw Sum & Mean Map for each neuron
    load([savedir_raster date '.RasterMean.mat']);
    
    [' ======= Start draw SumMap  ======= ']

    for n = 1 : size(raster,1) 
        clf;
        fig = figure('Visible', 'off');
        set(fig, 'Position', [1,1,960,960]);
        mean_value=mean_map(n,:,:);
        mean_value = reshape(mean_value,[size(target_unique,1),5]);
        ylim = max([max(max(mean_value)),0.4]);
        plot(mean_value','DisplayName','mean_list_curr_list','LineWidth',4)
        hold on;

        for i = 1:3
            plot([i+1,i+1],[0,ylim],'--r','Markersize',3,'LineWidth',2);
        end
        set(gca,'xtick',1:1:5);
        set(gca,'xticklabel',['Prep';' GO ';'Back';'Rewa';'Rest']);
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
        title([date ' CH' num2str(ceil(double(readNeurons(n))/6))  ' NeuronNo.' num2str(readNeurons(n))]);
        fig_name3 = [date ' CH' num2str(ceil(double(readNeurons(n))/6)) ' NeuronNo.' num2str(readNeurons(n)) 'MeanFR' '.tif']
        saveas(fig, fullfile(savedir_raster, fig_name3) );

    end

    [' ======= End draw SumMap  ======= ']

end
