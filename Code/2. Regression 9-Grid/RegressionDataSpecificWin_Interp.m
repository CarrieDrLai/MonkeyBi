function [dataset] = RegressionDataSpecificWin(save_dataset,hand,data0,window_bins,trial_id_shuffle,time_map,trajectory_map)
    %% Go data
    neural_Go = zeros([size(data0,1),window_bins,10000]);
    vel_Go = [];
    trial_mask = [];
    neural_idcount = 0;
    interp_num = 20;
    for t = 1:size(trial_id_shuffle,1)
        trajectory_curr = trajectory_map(trial_id_shuffle(t),find(time_map(trial_id_shuffle(t),:,2)>0),:);
        vel_curr = diff(trajectory_curr);
        vel_curr = reshape(vel_curr,[size(vel_curr,2) 4]);
        vel_curr_cut = vel_curr(3:end-5,:);
        vel_curr_cut_interp = [];
        for vi = 1:4
            x = (1:size(vel_curr_cut,1))'; 
            y = vel_curr_cut(:,vi);
            xi = (1:interp_num)';
            vel_curr_cut_interp = [vel_curr_cut_interp,interp1(x,y,xi,'linear','extrap')];    
        end
        
        vel_Go = [vel_Go;vel_curr_cut_interp];
        
        data0_curr = data0(:,trial_id_shuffle(t),50:(max(find(time_map(trial_id_shuffle(t),:,2)>0))-5));
        data0_curr = reshape(data0_curr,[size(data0_curr,1) size(data0_curr,3)]);
        data0_curr_interp = [];
        for ni = 1:size(data0_curr,1)
            x = (1:size(data0_curr,2))'; 
            y = data0_curr(ni,:)';
            xi = (1:interp_num+4)';
            data0_curr_interp = [data0_curr_interp;interp1(x,y,xi,'linear','extrap')'];    
        end
        
        for i = 1: size(vel_curr_cut_interp,1)
            neural_idcount = neural_idcount+1;
            neural_Go(:,:,neural_idcount) = data0_curr_interp(:,i:i+4);
        end
        
        trial_mask = [trial_mask;repmat(t,[size(vel_curr_cut_interp,1) 1])];
    end
    neural_Go = neural_Go(:,:,1:neural_idcount);
    
    %% Return data
    neural_Return =zeros([size(data0,1),window_bins,10000]);
    vel_Return = [];
    trial_mask_Return = [];
    neural_idcount = 0;
    for t = 1:size(trial_id_shuffle,1)
        trajectory_curr = trajectory_map(trial_id_shuffle(t),find(time_map(trial_id_shuffle(t),:,3)>0),:);
        vel_curr = diff(trajectory_curr);
        vel_curr = reshape(vel_curr,[size(vel_curr,2) 4]);
        vel_curr_cut = vel_curr(4:end-5,:);
        vel_curr_cut_interp = [];
        for vi = 1:4
            x = (1:size(vel_curr_cut,1))'; 
            y = vel_curr_cut(:,vi);
            xi = (1:interp_num)';
            vel_curr_cut_interp = [vel_curr_cut_interp,interp1(x,y,xi,'linear','extrap')];    
        end
        vel_Return = [vel_Return;vel_curr_cut_interp];
        
        data0_curr = data0(:,trial_id_shuffle(t),min(find(time_map(trial_id_shuffle(t),:,3)>0)):(max(find(time_map(trial_id_shuffle(t),:,3)>0))-5));
        data0_curr = reshape(data0_curr,[size(data0_curr,1) size(data0_curr,3)]);
        data0_curr_interp = [];
        for ni = 1:size(data0_curr,1)
            x = (1:size(data0_curr,2))'; 
            y = data0_curr(ni,:)';
            xi = (1:interp_num+4)';
            data0_curr_interp = [data0_curr_interp;interp1(x,y,xi,'linear','extrap')'];    
        end
        for i = 1: size(vel_curr_cut_interp,1)
            neural_idcount = neural_idcount+1;
            neural_Return(:,:,neural_idcount) = data0_curr_interp(:,i:i+4);
        end
        
        trial_mask_Return = [trial_mask_Return;repmat(t,[size(vel_curr_cut_interp,1) 1])];
    end
    neural_Return = neural_Return(:,:,1:neural_idcount);
    
    %% 5-fold mask
    trial_no = size(trial_id_shuffle,1);
    fold_trial_id = randperm(trial_no);			
    fold_mask = [];
    fold_mask_Return = [];
    for t = 1:trial_no
        curr_fold_bins = sum(trial_mask==t);
        fold_mask = [fold_mask;repmat(ceil(find(fold_trial_id==t)/ceil(trial_no/5)),[curr_fold_bins,1])];
        curr_fold_bins_Return = sum(trial_mask_Return==t);
        fold_mask_Return = [fold_mask_Return;repmat(ceil(find(fold_trial_id==t)/ceil(trial_no/5)),[curr_fold_bins_Return,1])];
    end    
    
    %% save
    trialmask = struct('trial_mask',trial_mask,'trial_mask_Return',trial_mask_Return);
    foldmask = struct('fold_mask',fold_mask,'fold_mask_Return',fold_mask_Return);
    data_Go = struct('neural_Go',neural_Go,'vel_Go',vel_Go);
    data_Return = struct('neural_Return',neural_Return,'vel_Return',vel_Return);
    dataset = struct('trialmask',trialmask,'foldmask',foldmask,'data_Go',data_Go,'data_Return',data_Return,'window_bins',window_bins,'trial_id_shuffle',trial_id_shuffle,'hand',hand);
    save([save_dataset 'dataset_' 'window' num2str(window_bins) '_' hand '.Interp.mat'],'dataset');
end