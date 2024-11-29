function [] = NineGridRegression(savedir,save_selectedtrial_dir,date_list,threshFR)

    %% config & path
    
    datafile_id = '1';
    session = '1';
    area = ['M1L';'PMd';'M1R'];
    area_ch= [0,96,64,96];
    
    %% 1. Prepare Dataset
    for date_i = 1:size(date_list,1)
       % 1.1 config & path
        date = date_list(date_i,:);
        save_pls_data_dir = [savedir date '\pls_data\'];
        if ~exist(save_pls_data_dir, 'dir')
            mkdir(save_pls_data_dir);
        end

        filename = [date '.DatafileID' datafile_id  '.S' session  '.' traject_source];
        save_autoselectedtrial_dir = [save_selectedtrial_dir filename '.Auto\'];    
        load([save_autoselectedtrial_dir date '.DataRasterGo.mat'])
        
        n_id_thresh = DataRasterGo.FiringRateGo.mean_fr_session>threshFR;
        raster = DataRasterGo.RasterGo.raster_go(n_id_thresh,:,:);
        trajectory_map = DataRasterGo.TrajectoryGo.trajectory_map;
        time_map = DataRasterGo.RasterGo.time_map_go;
        target_list = DataRasterGo.TargetGo.target_list;
        %target_unique = DataRasterGo.TargetGo.target_unique;
        %trial_no = DataRasterGo.TargetGo.trial_no;
        NeuronNo = DataRasterGo.RasterGo.NeuronNo(n_id_thresh);
        %before_go_binlen = DataRasterGo.RasterGo.before_go_binlen;
        %psth_binlen = DataRasterGo.RasterGo.psth_binlen;
        %target_num = DataRasterGo.TargetGo.target_num;
        %target_num_curr = [0;target_num];
        %target_map = DataRasterGo.TargetGo.target_map;

        % 1.2 context & neuron masks
%         NeuronMask = [];
%         for a = 1:size(area,1)
%             %NeuronMask = [NeuronMask;(NeuronNo_SDJudge <= 6*sum(area_ch(1:a+1)) & NeuronNo_SDJudge >= 1+6*sum(area_ch(1:a)))];
%             NeuronMask = [NeuronMask;(NeuronNo <= 6*sum(area_ch(1:a+1)) & NeuronNo >= 1+6*sum(area_ch(1:a)))];
%         end
    %     mask = find(NeuronNo_SDJudge <= 6*sum(area_ch(1:a+1)) & NeuronNo_SDJudge >= 1+6*sum(area_ch(1:a))) ;
    %     left_Nid_mask = ceil(NeuronNo/6)<97;
    %     right_Nid_mask = ceil(NeuronNo/6)>160;
    %     PMd_Nid_mask = (ceil(NeuronNo/6)<160 & ceil(NeuronNo/6)>96);
    %     NeuronMask = [left_Nid_mask;PMd_Nid_mask;right_Nid_mask]';

        left_hand_mask = mod(target_list,10)==0;
        right_hand_mask = (target_list-400)<10;
        bi_hand_id = find(~(left_hand_mask+right_hand_mask));
        left_hand_id = find(left_hand_mask);
        right_hand_id = find(right_hand_mask);

        left_hand_id_shuffle = left_hand_id(randperm(size(left_hand_id,1)));
        right_hand_id_shuffle = right_hand_id(randperm(size(right_hand_id,1)));
        bi_hand_id_shuffle = bi_hand_id(randperm(size(bi_hand_id,1)));

        % 1.3 dataset
        frequency = 50;
        data0 = frequency * smoothdata(raster,2,'gaussian',15); % smoothed spike count
        window_bins = 5; % 100ms window
        %var_num = 8;
        smooth_sigma = 1;
        smooth_win = 20;
        save_pls_data_dir = [savedir date '\pls_data\'];
        [dataset_L] = RegressionDataSpecificWin(save_pls_data_dir,'L',data0,target_list,window_bins,left_hand_id_shuffle,time_map,trajectory_map,smooth_sigma,smooth_win);
        [dataset_R] = RegressionDataSpecificWin(save_pls_data_dir,'R',data0,target_list,window_bins,right_hand_id_shuffle,time_map,trajectory_map,smooth_sigma,smooth_win);
        [dataset_Bi] = RegressionDataSpecificWin(save_pls_data_dir,'B',data0,target_list,window_bins,bi_hand_id_shuffle,time_map,trajectory_map,smooth_sigma,smooth_win);
    end
    
        
    %% 2.1 decode
    smooth_sigma = 1;
    smooth_win = 20;
    for date_i = 1:size(date_list,1)
        date = date_list(date_i,:);
        save_pls_decode_dir = [savedir date '\pls_decode\'];
        if ~exist(save_pls_decode_dir, 'dir')
            mkdir(save_pls_decode_dir);
        end        

        dataset_L = load([savedir date '\pls_data\' 'dataset_window5_L.mat']);
        dataset_L = dataset_L.dataset;
        dataset_R = load([savedir date '\pls_data\' 'dataset_window5_R.mat']);
        dataset_R = dataset_R.dataset;
        dataset_Bi = load([savedir date '\pls_data\' 'dataset_window5_B.mat']);
        dataset_Bi = dataset_Bi.dataset;
        % R_square_list: W-L, A-R, A-Bi; W-R, A-L, A-Bi; W-Bi, A-L, A-R.
        % A: Across; W:Within

        % 2.1.1 decode for single units
        filename = [date '.DatafileID' datafile_id  '.S' session  '.' traject_source];
        save_autoselectedtrial_dir = [save_selectedtrial_dir filename '.Auto\'];    
        load([save_autoselectedtrial_dir date '.DataRasterGo.mat'])
        n_id_thresh = DataRasterGo.FiringRateGo.mean_fr_session>threshFR;
        NeuronNo = DataRasterGo.RasterGo.NeuronNo(n_id_thresh);    
        [R_square_list,MSE_list,Coef_list] = RegressionPred(save_pls_decode_dir,dataset_L,dataset_R,dataset_Bi,NeuronNo,'unit',smooth_sigma,smooth_win);
        [R_square_list_shuffle,MSE_list_shuffle,Coef_list_shuffle] = RegressionPred(save_pls_decode_dir,dataset_L,dataset_R,dataset_Bi,NeuronNo,'unitShuffle',smooth_sigma,smooth_win);

        % 2.1.2 decode for area
        NeuronMask = [];
        for a = 1:size(area,1)
            %NeuronMask = [NeuronMask;(NeuronNo_SDJudge <= 6*sum(area_ch(1:a+1)) & NeuronNo_SDJudge >= 1+6*sum(area_ch(1:a)))];
            NeuronMask = [NeuronMask;(NeuronNo <= 6*sum(area_ch(1:a+1)) & NeuronNo >= 1+6*sum(area_ch(1:a)))];
        end
        NeuronMask = NeuronMask';
        [R_square_list_A,MSE_list_A,Coef_list_A] = RegressionPred(save_pls_decode_dir,dataset_L,dataset_R,dataset_Bi,NeuronMask,'area',smooth_sigma,smooth_win);
        [R_square_list_A_shuffle,MSE_list_A_shuffle,Coef_list_A_shuffle] = RegressionPred(save_pls_decode_dir,dataset_L,dataset_R,dataset_Bi,NeuronMask,'areaShuffle',smooth_sigma,smooth_win);

        % 2.1.3 decode for ensemble neurons
        NeuronMask2 = sum(NeuronMask,2);
        [R_square_list_E,MSE_list_E,Coef_list_E] = RegressionPred(save_pls_decode_dir,dataset_L,dataset_R,dataset_Bi,NeuronMask2,'ensemble',smooth_sigma,smooth_win);
        [R_square_list_E_shuffle,MSE_list_E_shuffle,Coef_list_E_shuffle] = RegressionPred(save_pls_decode_dir,dataset_L,dataset_R,dataset_Bi,NeuronMask2,'ensembleShuffle',smooth_sigma,smooth_win);

        % 2.1.4 decode for SD_Judge_units
        load([savedir '..\Process.Units\' 'Neuron_SDJudge\' date '_Neuron_SDJudge.mat']);
        [R_square_list_SDJudge,MSE_list_SDJudge,Coef_list_SDJudge] = RegressionPred(save_pls_decode_dir,dataset_L,dataset_R,dataset_Bi,NeuronMask_SDJudge,'SDJudge',smooth_sigma,smooth_win);
        [R_square_list_SDJudge_shuffle,MSE_list_SDJudge_shuffle,Coef_list_SDJudge_shuffle] = RegressionPred(save_pls_decode_dir,dataset_L,dataset_R,dataset_Bi,NeuronMask_SDJudge,'SDJudgeShuffle',smooth_sigma,smooth_win);

        % 2.2 classify units with coef_list
        cc_th = 0.1;
        NeuronMask_LOnly = sum(Coef_list(:,1,1:2)>cc_th,3)>0 & sum(Coef_list(:,4,3:4)>cc_th,3)==0;% Left Only
        sum(NeuronMask_LOnly)
    %     [R_square_list_group,MSE_list_group,Coef_list_group] = RegressionPred(save_pls_decode_dir,dataset_L,dataset_R,dataset_Bi,NeuronMask3,'group',smooth_sigma,smooth_win);
    %     [R_square_list_group_shuffle,MSE_list_group_shuffle,Coef_list_group_shuffle] = RegressionPred(save_pls_decode_dir,dataset_L,dataset_R,dataset_Bi,NeuronMask3,'groupShuffle',smooth_sigma,smooth_win);

        NeuronMask_ROnly = sum(Coef_list(:,1,1:2)>cc_th,3)==0 & sum(Coef_list(:,4,3:4)>cc_th,3)>0;% Right Only
        sum(NeuronMask_ROnly)
    %     [R_square_list_groupR,MSE_list_groupR,Coef_list_groupR] = RegressionPred(save_pls_decode_dir,dataset_L,dataset_R,dataset_Bi,NeuronMask4,'groupR',smooth_sigma,smooth_win);
    %     [R_square_list_groupR_shuffle,MSE_list_groupR_shuffle,Coef_list_groupR_shuffle] = RegressionPred(save_pls_decode_dir,dataset_L,dataset_R,dataset_Bi,NeuronMask4,'groupRShuffle',smooth_sigma,smooth_win);

        NeuronMask_BOnly = sum(Coef_list(:,1,1:2)>cc_th,3)>0 & sum(Coef_list(:,4,3:4)>cc_th,3)>0;% Both hand
        sum(NeuronMask_BOnly)
    %     [R_square_list_groupBoth,MSE_list_groupBoth,Coef_list_groupBoth] = RegressionPred(save_pls_decode_dir,dataset_L,dataset_R,dataset_Bi,NeuronMask5,'groupBoth',smooth_sigma,smooth_win);
    %     [R_square_list_groupBoth_shuffle,MSE_list_groupBoth_shuffle,Coef_list_groupBoth_shuffle] = RegressionPred(save_pls_decode_dir,dataset_L,dataset_R,dataset_Bi,NeuronMask5,'groupBothShuffle',smooth_sigma,smooth_win);

        % 2.3 Coef list for SDJudge
        Coef_SDJudge_units = Coef_list(NeuronMask_SDJudge,:,:);
        Coef_SDJudge_units_shuffle = Coef_list_shuffle(NeuronMask_SDJudge,:,:);

        save([save_pls_decode_dir 'SDJudge\'  'NeuronMask_SDJudge.mat'],'Coef_SDJudge_units','Coef_SDJudge_units_shuffle','NeuronNo','NeuronNo_SDJudge','NeuronMask_SDJudge','NeuronMask_LOnly','NeuronMask_ROnly','NeuronMask_BOnly');    
    %     save_pls_decode_dir = [savedir date '\pls_decode\'];
    %     dataset_L = load([savedir date '\pls_data\' 'dataset_window5_L.mat']);
    %     dataset_L = dataset_L.dataset;
    %     dataset_R = load([savedir date '\pls_data\' 'dataset_window5_R.mat']);
    %     dataset_R = dataset_R.dataset;
    %     dataset_Bi = load([savedir date '\pls_data\' 'dataset_window5_B.mat']);
    %     dataset_Bi = dataset_Bi.dataset;
    end
    
    %% 3. plot for 2.1.4

    fold_name = ['groupR\' 'TrainR.HandWithin.Fold3'];
    fold_name = ['area.RM1\' 'TrainR.HandWithin.Fold4'];
    
    for fold_i = 1:5
        fold_name = ['SDJudge\' 'TrainB.HandWithin.Fold' num2str(fold_i)];
        load([save_pls_decode_dir  fold_name '.mat']);
        aa = [test_y,y_predict];
        var_i = 1;    
        var_name = ['posxL';'posyL';'posxR';'posyR'];
        %var_name = ['vxL';'vyL';'vxR';'vyR'];
        trial_plot = unique(trial_mask);
        if ~exist([save_pls_decode_dir fold_name '\'], 'dir')
            mkdir([save_pls_decode_dir fold_name '\']);
        end
        ylim = 0.5;
        for i = 1:size(trial_plot,1)
            plot_id = find(trial_mask == trial_plot(i));
            phase_curr = phase_mask(plot_id);
            fig = figure('Visible', 'off');
            %vx=aa(plot_id,[var_i,var_i+4]);
            posxL=aa(plot_id,[var_i,var_i+8]);
            plot(posxL);
            hold on;
            plot([max(find(phase_curr==2)) max(find(phase_curr==2))],[0.5-ylim 0.5+ylim],'r-');
            hold on;
            plot([max(find(phase_curr==3)) max(find(phase_curr==3))],[0.5-ylim 0.5+ylim],'r-');
            set(gca,'YLim',[0.5-ylim 0.5+ylim]);  
            title(['Trial ' num2str(trial_plot(i)) ';Tar ' num2str(target_mask(plot_id(1))) ';CC: ' num2str(corr(aa(plot_id,var_i),aa(plot_id,var_i+8)))]);
            saveas(fig,[save_pls_decode_dir fold_name '\' 'Trial ' num2str(trial_plot(i)) ';Tar ' num2str(target_mask(plot_id(1))) var_name(var_i,:)  '.tif']);
            clf;

            fig2 = figure('Visible', 'off');
            posyL=aa(plot_id,[var_i+1,var_i+1+8]);
            plot(posyL);
            hold on;
            plot([max(find(phase_curr==2)) max(find(phase_curr==2))],[0.5-ylim 0.5+ylim],'r-');
            hold on;
            plot([max(find(phase_curr==3)) max(find(phase_curr==3))],[0.5-ylim 0.5+ylim],'r-');
            set(gca,'YLim',[0.5-ylim 0.5+ylim]);  
            title(['Trial ' num2str(trial_plot(i)) ';Tar ' num2str(target_mask(plot_id(1))) ';CC: ' num2str(corr(aa(plot_id,var_i+1),aa(plot_id,var_i+1+8)))]);
            saveas(fig2,[save_pls_decode_dir fold_name '\' 'Trial ' num2str(trial_plot(i)) ';Tar ' num2str(target_mask(plot_id(1))) var_name(var_i+1,:)  '.tif']);
            clf;

            fig3 = figure('Visible', 'off');
            posxR=aa(plot_id,[var_i+2,var_i+2+8]);
            plot(posxR);
            hold on;
            plot([max(find(phase_curr==2)) max(find(phase_curr==2))],[0.5-ylim 0.5+ylim],'r-');
            hold on;
            plot([max(find(phase_curr==3)) max(find(phase_curr==3))],[0.5-ylim 0.5+ylim],'r-');
            set(gca,'YLim',[0.5-ylim 0.5+ylim]);  
            title(['Trial ' num2str(trial_plot(i)) ';Tar ' num2str(target_mask(plot_id(1))) ';CC: ' num2str(corr(aa(plot_id,var_i+2),aa(plot_id,var_i+2+8)))]);
            saveas(fig3,[save_pls_decode_dir fold_name '\' 'Trial ' num2str(trial_plot(i)) ';Tar ' num2str(target_mask(plot_id(1))) var_name(var_i+2,:)  '.tif']);
            clf;

            fig4 = figure('Visible', 'off');
            posyR=aa(plot_id,[var_i+3,var_i+3+8]);
            plot(posyR);
            hold on;
            plot([max(find(phase_curr==2)) max(find(phase_curr==2))],[0.5-ylim 0.5+ylim],'r-');
            hold on;
            plot([max(find(phase_curr==3)) max(find(phase_curr==3))],[0.5-ylim 0.5+ylim],'r-');
            set(gca,'YLim',[0.5-ylim 0.5+ylim]);  
            title(['Trial ' num2str(trial_plot(i)) ';Tar ' num2str(target_mask(plot_id(1))) ';CC: ' num2str(corr(aa(plot_id,var_i+3),aa(plot_id,var_i+3+8)))]);
            saveas(fig4,[save_pls_decode_dir fold_name '\' 'Trial ' num2str(trial_plot(i)) ';Tar ' num2str(target_mask(plot_id(1))) var_name(var_i+3,:)  '.tif']);
            clf;

            C = fix(rem(rem(target_mask(plot_id(1)),100),10));
            B = fix((rem(target_mask(plot_id(1)),100) - C)/10);

    %         if B>0 & C==0
    %             traj = [0.26,0.5];
    %             traj_pred = [0.26,0.5];
    %         end
    %         if B==0 & C>0
    %             traj = [0.74,0.5];
    %             traj_pred = [0.74,0.5];
    %         end
    %         
    %             
    %         for pos_i = 1:size(plot_id,1)
    %             pos_curr = traj(end,:) + aa(plot_id(pos_i),var_i:var_i+1);
    %             traj = [traj;pos_curr];
    %         end
    %         
    %         
    %         for pos_i = 1:size(plot_id,1)
    %             pos_curr = traj_pred(end,:) + aa(plot_id(pos_i),var_i+4:var_i+1+4);
    %             traj_pred = [traj_pred;pos_curr];
    %         end

            fig5 = figure('Visible', 'off');
            %plot(traj(:,1),traj(:,2),'bo');
            plot(posxL(:,1),posyL(:,1),'b-'); % true
            hold on;
            plot(posxL(:,2),posyL(:,2),'r-'); % pred
            hold on;
            plot(posxR(:,1),posyR(:,1),'bo'); % true
            hold on;
            %plot(traj_pred(:,1),traj_pred(:,2),'ro');
            plot(posxR(:,2),posyR(:,2),'ro'); % pred
            set(gca,'YLim',[0 1]);  
            set(gca,'XLim',[0 1]);  
            title(['Trial ' num2str(trial_plot(i)) ';Tar ' num2str(target_mask(plot_id(1))) ]);
            saveas(fig5,[save_pls_decode_dir fold_name '\' 'Trial ' num2str(trial_plot(i)) ';Tar ' num2str(target_mask(plot_id(1))) 'Traj'  '.tif']);
            clf;
            save([save_pls_decode_dir fold_name '\' 'Trial ' num2str(trial_plot(i)) ';Tar ' num2str(target_mask(plot_id(1)))  '.mat'],'phase_curr','posxL','posyL','posxR','posyR');
        end
    end
    
    %% 4. single units, within-hand
%     neuron_type = Coef_list(:,[1,4,7],:)>0.1;
%     neuron_type = [a(:,[1,3],[1,2]),a(:,[2,3],[3,4])];
%     
%     % single units: within-hand t test
%     [h_within_L,p_within_L]=ttest(R_square_list(:,1,:),R_square_list_shuffle(:,1,:)); % Left Hand
%     [h_within_R,p_within_R]=ttest(R_square_list(:,4,:),R_square_list_shuffle(:,4,:)); % Right Hand
%     [h_within_B,p_within_B]=ttest(R_square_list(:,7,:),R_square_list_shuffle(:,7,:)); % Bi Hand
%     save([save_pls_decode_dir 'neuron_contribution_p.mat'],'p_within_L','p_within_R','p_within_B','NeuronNo','NeuronMask');
%     % single units: within-hand R-square
%     R_square_LM_WithinL = [reshape(R_square_list(left_Nid_mask,1,:),[sum(left_Nid_mask),4]),reshape(R_square_list_shuffle(left_Nid_mask,1,:),[sum(left_Nid_mask),4])];
%     R_square_LM_WithinR = [reshape(R_square_list(left_Nid_mask,4,:),[sum(left_Nid_mask),4]),reshape(R_square_list_shuffle(left_Nid_mask,4,:),[sum(left_Nid_mask),4])];
%     R_square_LM_WithinBi = [reshape(R_square_list(left_Nid_mask,7,:),[sum(left_Nid_mask),4]),reshape(R_square_list_shuffle(left_Nid_mask,7,:),[sum(left_Nid_mask),4])];
%     R_square_RM_WithinL = [reshape(R_square_list(right_Nid_mask,1,:),[sum(right_Nid_mask),4]),reshape(R_square_list_shuffle(right_Nid_mask,1,:),[sum(right_Nid_mask),4])];
%     R_square_RM_WithinR = [reshape(R_square_list(right_Nid_mask,4,:),[sum(right_Nid_mask),4]),reshape(R_square_list_shuffle(right_Nid_mask,4,:),[sum(right_Nid_mask),4])];
%     R_square_RM_WithinBi = [reshape(R_square_list(right_Nid_mask,7,:),[sum(right_Nid_mask),4]),reshape(R_square_list_shuffle(right_Nid_mask,7,:),[sum(right_Nid_mask),4])];
%     R_square_PMd_WithinL = [reshape(R_square_list(PMd_Nid_mask,1,:),[sum(PMd_Nid_mask),4]),reshape(R_square_list_shuffle(PMd_Nid_mask,1,:),[sum(PMd_Nid_mask),4])];
%     R_square_PMd_WithinR  = [reshape(R_square_list(PMd_Nid_mask,4,:),[sum(PMd_Nid_mask),4]),reshape(R_square_list_shuffle(PMd_Nid_mask,4,:),[sum(PMd_Nid_mask),4])];
%     R_square_PMd_WithinBi = [reshape(R_square_list(PMd_Nid_mask,7,:),[sum(PMd_Nid_mask),4]),reshape(R_square_list_shuffle(PMd_Nid_mask,7,:),[sum(PMd_Nid_mask),4])];
%     save([save_pls_decode_dir 'neuron_contribution_R.mat'],'R_square_LM_WithinL','R_square_LM_WithinR','R_square_LM_WithinBi','R_square_RM_WithinL','R_square_RM_WithinR','R_square_RM_WithinBi','R_square_PMd_WithinL','R_square_PMd_WithinR','R_square_PMd_WithinBi','NeuronNo','NeuronMask');
%     
%     
%     % classification for the neuron type : pie chart
%     % VxL: W-L,W-R,W-B; VyL: W-L,W-R,W-B; VxR: W-L,W-R,W-B; VyR: W-L,W-R,W-B
%     R_square_classify_neurons = reshape(R_square_list(:,[1,4,7],:),[size(R_square_list,1),3*4]) % Neuron * (3,3,3,3)�?4个变量，3种运动状�?
%     R_square_classify_neurons_LM = R_square_classify_neurons(left_Nid_mask,:);
%     R_square_classify_neurons_LM_mask = R_square_classify_neurons_LM>0.1;
%     R_square_classify_neurons_PMd = R_square_classify_neurons(PMd_Nid_mask,:);
%     R_square_classify_neurons_PMd_mask = R_square_classify_neurons_PMd>0.1;
%     R_square_classify_neurons_RM = R_square_classify_neurons(right_Nid_mask,:);
%     R_square_classify_neurons_RM_mask = R_square_classify_neurons_RM>0.1;
%     save([save_pls_decode_dir 'neuron_classify.mat'],'R_square_classify_neurons_LM','R_square_classify_neurons_LM_mask','R_square_classify_neurons_PMd','R_square_classify_neurons_PMd_mask','R_square_classify_neurons_RM','R_square_classify_neurons_RM_mask');

    %% 5. area: within-hand t test

%     [h_within_L_LM,p_within_L_LM]=ttest(R_square_list_A(1,1,:),R_square_list_A_shuffle(1,1,:)); % Left Hand
%     [h_within_R_LM,p_within_R_LM]=ttest(R_square_list_A(1,4,:),R_square_list_A_shuffle(1,4,:)); % Right Hand
%     [h_within_B_LM,p_within_B_LM]=ttest(R_square_list_A(1,7,:),R_square_list_A_shuffle(1,7,:)); % Bi Hand
%     [h_within_L_PMd,p_within_L_PMd]=ttest(R_square_list_A(2,1,:),R_square_list_A_shuffle(2,1,:)); % Left Hand
%     [h_within_R_PMd,p_within_R_PMd]=ttest(R_square_list_A(2,4,:),R_square_list_A_shuffle(2,4,:)); % Right Hand
%     [h_within_B_PMd,p_within_B_PMd]=ttest(R_square_list_A(2,7,:),R_square_list_A_shuffle(2,7,:)); % Bi Hand
%     [h_within_L_RM,p_within_L_RM]=ttest(R_square_list_A(3,1,:),R_square_list_A_shuffle(3,1,:)); % Left Hand
%     [h_within_R_RM,p_within_R_RM]=ttest(R_square_list_A(3,4,:),R_square_list_A_shuffle(3,4,:)); % Right Hand
%     [h_within_B_RM,p_within_B_RM]=ttest(R_square_list_A(3,7,:),R_square_list_A_shuffle(3,7,:)); % Bi Hand

    % save
    
    %% 6. test, plot example raster (after interp)
%     interp_num0 = 20;
%     for n = 1:size(NeuronNo,2)
%         raster_nid = [];
%         for i = 1:size(raster,2)
%             go_id = find(time_map(i,:,2)>0);
%             raster0 = raster(n,i,go_id(1:end-4));
%             raster0 = reshape(raster0,[size(raster0,3),1]);
%             x = (1:size(raster0,1))'; 
%             xi = (1:interp_num0)';
%             data_go = interp1(x,raster0,xi,'linear','extrap')';
%             return_id = find(time_map(i,:,3)>0);
%             raster_r = raster(n,i,return_id(1:end-4));
%             raster_r = reshape(raster_r,[size(raster_r,3),1]);
%             x = (1:size(raster_r,1))'; 
%             xi = (1:interp_num0)';
%             data_return = interp1(x,raster_r,xi,'linear','extrap')';
%             raster_nid = [raster_nid;[data_go,data_return]];                
%         end
%     end
%     
%     mean_fr=[];
%     for  t = 1:15
%         mean_fr = [mean_fr;mean(raster_nid(sum(target_num_curr(1:t))+1:sum(target_num_curr(1:t+1)),:),1)];
%     end
%     fig = figure('Visible', 'off');
%     set(fig, 'Position', [1,1,960,960]);
%     plot(smoothdata(mean_fr',2,'gaussian',15));
%     saveas(fig, fullfile(save_pls_decode_dir, 'mean_fr.jpg') );
%     mean_fr_meanphase = [mean(mean_fr(:,1:20),2),mean(mean_fr(:,21:40),2)];
    
    % save
    
    %% 4. plt pred
%     clf;
%     fig = figure('Visible', 'off');
%     set(fig, 'Position', [1,1,960,960]);
%     plot(parameters(:,3:4));
%     saveas(fig, fullfile(save_pls_data_dir, 'Rsquare.jpg') );


%     for window_bins = 1:15
%         clf;
%         fig = figure('Visible', 'off');
%         set(fig, 'Position', [1,1,960,960]);
%         plot(parameters(parameters(:,1)==window_bins,7:end));
%         saveas(fig, fullfile(save_pls_data_dir, ['window' num2str(window_bins) '.jpg']) );    
%     end
end