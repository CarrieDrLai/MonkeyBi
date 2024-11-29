function [dataset] = RegressionDataSpecificWin(save_dataset,hand,data0,target_list,window_bins,trial_id_shuffle,time_map,trajectory_map,smooth_sigma,smooth_win)
    %% Go data
    neural = zeros([size(data0,1),window_bins,10000]);
    vel = [];
    trajectory = [];
    trial_mask = [];
    target_mask = [];
    phase_mask = [];
    neural_idcount = 0;
    interp_num = 20;
    start_bin =  49;
    del_bins = 5;
    del_lag = 0;

    for t = 1:size(trial_id_shuffle,1)
        trajectory_curr = [trajectory_map(trial_id_shuffle(t),[find(time_map(trial_id_shuffle(t),:,2)>0),find(time_map(trial_id_shuffle(t),:,3)>0)],:)];
        %trajectory_curr = [trajectory_map(trial_id_shuffle(t),[find(time_map(trial_id_shuffle(t),:,2)>0),find(time_map(trial_id_shuffle(t),:,3)>0),find(time_map(trial_id_shuffle(t),:,4)>0),find(time_map(trial_id_shuffle(t),:,5)>0)],:)];
        trajectory_curr = reshape(trajectory_curr,[size(trajectory_curr,2) 4]);
        
        vel_curr = diff(trajectory_curr);
        %vel_curr_cut =
        %[zeros([1,size(vel_curr,2)]);vel_curr(1-del_lag:end-del_bins-del_lag,:)];
        %V0ฒน0
        vel_curr_cut = vel_curr(1-del_lag:end-del_bins-del_lag,:);
        vel_curr_cut = smoothdata(vel_curr_cut,smooth_sigma,'gaussian',smooth_win);
        vel = [vel;vel_curr_cut];
%         vel_curr_cut_interp = [];
%         for vi = 1:4
%             x = (1:size(vel_curr_cut,1))'; 
%             y = vel_curr_cut(:,vi);
%             xi = (1:interp_num)';
%             vel_curr_cut_interp = [vel_curr_cut_interp,interp1(x,y,xi,'linear','extrap')];    
%         end
        
        trajectory_curr_cut = trajectory_curr(2-del_lag:end-del_bins-del_lag,:);
        trajectory_curr_cut = smoothdata(trajectory_curr_cut,smooth_sigma,'gaussian',smooth_win);
        trajectory = [trajectory;trajectory_curr_cut];
        
        if size(max(find(time_map(trial_id_shuffle(t),:,5)>0)),2)>0
            data0_curr = data0(:,trial_id_shuffle(t),start_bin:(max(find(time_map(trial_id_shuffle(t),:,5)>0))-del_bins));
        else
            data0_curr = data0(:,trial_id_shuffle(t),start_bin:(max(find(time_map(trial_id_shuffle(t),:,4)>0))-del_bins));
        end
        data0_curr = reshape(data0_curr,[size(data0_curr,1) size(data0_curr,3)]);
        data0_curr = smoothdata(data0_curr,smooth_sigma,'gaussian',smooth_win);
%         data0_curr_interp = [];
%         for ni = 1:size(data0_curr,1)
%             x = (1:size(data0_curr,2))'; 
%             y = data0_curr(ni,:)';
%             xi = (1:interp_num+4)';
%             data0_curr_interp = [data0_curr_interp;interp1(x,y,xi,'linear','extrap')'];    
%         end
        
        for i = 1: size(vel_curr_cut,1)
            neural_idcount = neural_idcount+1;
            neural(:,:,neural_idcount) = data0_curr(:,i:i+4);
        end
        
        phase_mask0 = zeros([size(vel_curr_cut,1) 1]);
        %max_len = min(size(vel_curr_cut,1),size([find(time_map(trial_id_shuffle(t),:,2)>0),find(time_map(trial_id_shuffle(t),:,3)>0),find(time_map(trial_id_shuffle(t),:,4)>0),find(time_map(trial_id_shuffle(t),:,5)>0)],2));
        max_len = size(vel_curr_cut,1);
        phase_mask0(1:size(find(time_map(trial_id_shuffle(t),:,2)>0),2))=2;
        phase_mask0(size(find(time_map(trial_id_shuffle(t),:,2))>0,2)+1:max_len)=3;
        %phase_mask0(size(find(time_map(trial_id_shuffle(t),:,2)>0),2)+1:size([find(time_map(trial_id_shuffle(t),:,2)>0),find(time_map(trial_id_shuffle(t),:,3)>0)],2))=3;
        %phase_mask0(size([find(time_map(trial_id_shuffle(t),:,2)>0),find(time_map(trial_id_shuffle(t),:,3)>0)],2)+1:max_len)=3;

%         if size(max(find(time_map(trial_id_shuffle(t),:,5)>0)),2)>0
%             phase_mask0(size([find(time_map(trial_id_shuffle(t),:,2)>0),find(time_map(trial_id_shuffle(t),:,3)>0)],2)+1:size([find(time_map(trial_id_shuffle(t),:,2)>0),find(time_map(trial_id_shuffle(t),:,3)>0),find(time_map(trial_id_shuffle(t),:,4)>0)],2))=4;
%             phase_mask0(size([find(time_map(trial_id_shuffle(t),:,2)>0),find(time_map(trial_id_shuffle(t),:,3)>0),find(time_map(trial_id_shuffle(t),:,4)>0)],2)+1:max_len)=5;
%         else
%             phase_mask0(size([find(time_map(trial_id_shuffle(t),:,2)>0),find(time_map(trial_id_shuffle(t),:,3)>0)],2)+1:max_len)=4;
%         end
            
        %[size(phase_mask0,1),size(vel_curr_cut,1)]
        phase_mask = [phase_mask;phase_mask0];     
        trial_mask = [trial_mask;repmat(t,[size(vel_curr_cut,1) 1])];
        target_mask = [target_mask;repmat(target_list(trial_id_shuffle(t)),[size(vel_curr_cut,1) 1])];
    end
    neural = neural(:,:,1:neural_idcount);
    
    %% Return data
%     neural_Return =zeros([size(data0,1),window_bins,10000]);
%     vel_Return = [];
%     trial_mask_Return = [];
%     neural_idcount = 0;
%     for t = 1:size(trial_id_shuffle,1)
%         trajectory_curr = trajectory_map(trial_id_shuffle(t),find(time_map(trial_id_shuffle(t),:,3)>0),:);
%         vel_curr = diff(trajectory_curr);
%         vel_curr = reshape(vel_curr,[size(vel_curr,2) 4]);
%         vel_curr_cut = vel_curr(4:end-5,:);
%         vel_Return = [vel_Return;vel_curr_cut];
%         
%         data0_curr = data0(:,trial_id_shuffle(t),min(find(time_map(trial_id_shuffle(t),:,3)>0)):(max(find(time_map(trial_id_shuffle(t),:,3)>0))-5));
%         data0_curr = reshape(data0_curr,[size(data0_curr,1) size(data0_curr,3)]);
%         for i = 1: size(vel_curr_cut,1)
%             neural_idcount = neural_idcount+1;
%             neural_Return(:,:,neural_idcount) = data0_curr(:,i:i+4);
%         end
%         
%         trial_mask_Return = [trial_mask_Return;repmat(t,[size(vel_curr_cut,1) 1])];
%     end
%     neural_Return = neural_Return(:,:,1:neural_idcount);
    
    %% 5-fold mask
    fold_no = 5;
    trial_no = size(trial_id_shuffle,1);
    fold_trial_id = randperm(trial_no);			
    fold_mask = [];
    %fold_mask_Return = [];
    for t = 1:trial_no
        curr_fold_bins = sum(trial_mask==t);
        fold_mask = [fold_mask;repmat(ceil(find(fold_trial_id==t)/ceil(trial_no/fold_no)),[curr_fold_bins,1])];
%         curr_fold_bins_Return = sum(trial_mask_Return==t);
%         fold_mask_Return = [fold_mask_Return;repmat(ceil(find(fold_trial_id==t)/ceil(trial_no/5)),[curr_fold_bins_Return,1])];
    
    end    
    
    %% save
%     trialmask = struct('trial_mask',trial_mask,'trial_mask_Return',trial_mask_Return);
%     foldmask = struct('fold_mask',fold_mask,'fold_mask_Return',fold_mask_Return);
%     data_Go = struct('neural_Go',neural_Go,'vel_Go',vel_Go);
%     data_Return = struct('neural_Return',neural_Return,'vel_Return',vel_Return);
%     dataset = struct('trialmask',trialmask,'foldmask',foldmask,'data_Go',data_Go,'data_Return',data_Return,'window_bins',window_bins,'trial_id_shuffle',trial_id_shuffle,'hand',hand);
    
    targetmask = struct('target_mask',target_mask);
    trialmask = struct('trial_mask',trial_mask,'phase_mask',phase_mask);
    foldmask = struct('fold_mask',fold_mask);
    data_Go = struct('neural_Go',neural,'vel_Go',vel,'traj_Go',trajectory);
    %data_Return = struct('neural_Return',neural_Return,'vel_Return',vel_Return);
    dataset = struct('targetmask',targetmask,'trialmask',trialmask,'foldmask',foldmask,'data_Go',data_Go,'window_bins',window_bins,'trial_id_shuffle',trial_id_shuffle,'hand',hand);
    save([save_dataset 'dataset_' 'window' num2str(window_bins) '_' hand '.mat'],'dataset');
end