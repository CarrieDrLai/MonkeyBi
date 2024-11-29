function [] = HistReturn(date,area,area_ch,savedir_raster)
    %load([data_dir_session date '.Session.' session '.mat']);
    %load([data_dir date '.S' num2str(session) '\' date '.Session.' session '.mat']);
    
    load([savedir_raster date '.DataRasterReturn.mat']);
    
    target_num = DataRasterReturn.TargetReturn.target_num;
    target_num_curr = [0;target_num];
    target_unique = DataRasterReturn.TargetReturn.target_unique;
    target_list = DataRasterReturn.TargetReturn.target_list;

    mean_fr_session = DataRasterReturn.FiringRateReturn.mean_fr_session;
    readNeurons = DataRasterReturn.RasterReturn.NeuronNo;
    
    raster = DataRasterReturn.RasterReturn.raster_return;
    time_map = DataRasterReturn.RasterReturn.time_map_return;

    psth_binlen = DataRasterReturn.RasterReturn.psth_binlen;
    before_return_binlen = DataRasterReturn.RasterReturn.before_return_binlen;


    %% Hist Return Alignment

    [' ======= Start draw Hist -- Return Phase ======= ']
    raster_return = DataRasterReturn.RasterReturn.raster_return;
    time_map_return = DataRasterReturn.RasterReturn.time_map_return;

    hist_return = zeros([size(raster_return,1),size(target_unique,1),psth_binlen]); %NeuronNum*8*binlen

    for n = 1 : size(raster_return,1)
        for t = 1 : size(target_unique,1)
            spike_curr = raster_return(n,(sum(target_num_curr(1:t))+1):sum(target_num_curr(1:t+1)),:);
            time_map_curr = sum(time_map_return((sum(target_num_curr(1:t))+1):sum(target_num_curr(1:t+1)),:,:),3);
            spike_curr0 = reshape(spike_curr,[size(spike_curr,2),size(spike_curr,3)]).* time_map_curr;
            hist_return(n,t,:) =  sum(spike_curr0,1)./ (sum(time_map_curr,1)+0.00001);
        end
    end

    for n = 1 : size(raster_return,1)
        clf;
        fig = figure('Visible', 'off');
        set(fig, 'Position', [1,1,960,960]);
        value=hist_return(n,:,:);
        value = reshape(value,[size(hist_return,2),size(hist_return,3)]);

        %plot(value','DisplayName','hist_go','LineWidth',2)
        %smoothdata(value,2,'gaussian',15)
        value = smoothdata(value,2,'gaussian',15);
        ylim = max([max(max(value)),0.6]);
        plot(value','DisplayName','hist_return','LineWidth',2)
        hold on;

        plot([before_return_binlen,before_return_binlen],[0,ylim],'--r','Markersize',3,'LineWidth',2);
        set(gca,'xtick',[before_return_binlen]);
        set(gca,'xticklabel',['Return']);
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
        fig_name5 = [date  ' CH' num2str(ceil(double(readNeurons(n))/6)) ' NeuronNo.' num2str(readNeurons(n)) ' HistReturn' '.tif']
        saveas(fig, fullfile(savedir_raster, fig_name5) );

    end
    [' ======= End draw Hist -- Return Phase ======= ']
end
