function [R_sqr,MSE,CC] = RegressionFunc(save_pls_decode_dir,dataset,hand,fold_id,current_id,dataset_across1,dataset_across2,savedirkey,var_num,smooth_sigma,smooth_win)
    

    save_dir = [save_pls_decode_dir savedirkey '\'];
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    
    data_curr = dataset.data_Go.neural_Go(current_id,:,:);
    data_curr = reshape(data_curr,[size(data_curr,1)*size(data_curr,2),size(data_curr,3)]);
    data_curr = data_curr';
    if length(savedirkey)>10
        data_curr_shuffle = data_curr(randperm(size(data_curr,1)),:);
        traj_shuffle = dataset.data_Go.traj_Go(randperm(size(data_curr,1)),:);
        vel_shuffle = dataset.data_Go.vel_Go(randperm(size(data_curr,1)),:);
        data_curr = [data_curr_shuffle,traj_shuffle,vel_shuffle];
    else
        data_curr = [data_curr,dataset.data_Go.traj_Go,dataset.data_Go.vel_Go];
    end
    
    mu = mean(data_curr);
    mu(end-7) = 0.26;
    mu(end-6) = 0.5;
    mu(end-5) = 0.74;
    mu(end-4) = 0.5;
    sig = std(data_curr);
    %rr = corrcoef(data_curr);
    %zvalue = zscore(data_curr); 
    zvalue = [];
    for i = 1:size(data_curr,2)
        zvalue = [zvalue,zscore(data_curr(:,i))];
    end
    
    
    fold_mask = (dataset.foldmask.fold_mask == fold_id);
    phase_mask = (dataset.trialmask.phase_mask(fold_mask));
    trial_mask = (dataset.trialmask.trial_mask(fold_mask));
    target_mask =  (dataset.targetmask.target_mask(fold_mask));
    
    test_x = zvalue(fold_mask,1:(end-var_num));
    test_y = data_curr(fold_mask,(end-var_num+1):end);
    %test_y = zvalue(fold_mask,(end-3):end);
    train_x = zvalue(~fold_mask,1:(end-var_num));
    train_y = zvalue(~fold_mask,(end-var_num+1):end); 
    
    % pls 
    [XL, YL, XS, YS, BETA, PCTVAR, MSE, stats] = plsregress(train_x,train_y);
%     xw = train_x\XS; % 自变量提出成分的系数
%     yw = train_y\YS; % 因变量提出成分的系数
%     x_0 = PCTVAR(1,:); % 自变量提出成分的贡献�?
%     y_0 = PCTVAR(2,:); % 因变量提出成分的贡献�?
%     x_1 = sum(x_0); % 累计贡献�?
%     y_1 = sum(y_0); % 累计贡献�?    

    x_num = size(train_x,2);
    y_num = size(train_y,2);
    beta = zeros([x_num+1,y_num]); % 估计值，K�?
    beta(1,:) = mu(x_num+1:end)-mu(1:x_num)./sig(1:x_num)*BETA([2:end],:).*sig(x_num+1:end); % 估计值，b尖，常数�?
    beta([2:x_num+1],:) = (1./sig(1:x_num))'*sig(x_num+1:end).*BETA([2:end],:); % 估计值，k尖，变量系数
    test_len = size(test_x,1);
    b_predict = repmat(beta(1,:),[test_len,1]);
    k_predict = beta(2:end,:);
    y_predict = b_predict +  test_x*k_predict;
    
    y_predict = smoothdata(y_predict,smooth_sigma,'gaussian',smooth_win);
    
    [R_sqr_within,MSE_within,CC_within] = RegressionFuncEvaluation(y_num,test_len,y_predict,test_y);
    
    save([save_dir 'Train' hand '.HandWithin.' 'Fold' num2str(fold_id) '.mat'],'y_predict','test_y','R_sqr_within','MSE_within','CC_within','trial_mask','target_mask','phase_mask','test_x','b_predict','k_predict');
    %plotregression
    % figure;
    %plotregression(y_pre_1,y_test_1)
    %e1 = y_test_1 - y_pre_1;
    %e2 = y_test_2 - y_pre_2; 
    
    %% across-hand 
    [R_sqr_across1,MSE_across1,CC_across1] = RegressionFuncAcrossPred(save_dir,dataset_across1,hand,'1',current_id,fold_id,beta,y_num,var_num,smooth_sigma,smooth_win);
    [R_sqr_across2,MSE_across2,CC_across2] = RegressionFuncAcrossPred(save_dir,dataset_across2,hand,'2',current_id,fold_id,beta,y_num,var_num,smooth_sigma,smooth_win);
    
    R_sqr = [R_sqr_within;R_sqr_across1;R_sqr_across2];
    MSE = [MSE_within;MSE_across1;MSE_across2];
    CC = [CC_within;CC_across1;CC_across2];
    
    %save([save_dir 'Train' hand '.Evaluate.' 'Fold' num2str(fold_id) '.mat'],'R_sqr','MSE','CC');
end