function [R_sqr_across,MSE_across,CC_across] = RegressionFuncAcrossPred(save_dir,dataset_across,hand,across_id,n_id,fold_id,beta,y_num,var_num,smooth_sigma,smooth_win)

    data_curr_across = dataset_across.data_Go.neural_Go(n_id,:,:);
    data_curr_across = reshape(data_curr_across,[size(data_curr_across,1)*size(data_curr_across,2),size(data_curr_across,3)]);
    data_curr_across = data_curr_across';
    data_curr_across = [data_curr_across,dataset_across.data_Go.traj_Go,dataset_across.data_Go.vel_Go];
    zvalue_across = zscore(data_curr_across); 
    
    fold_mask_across = (dataset_across.foldmask.fold_mask == fold_id);
    trial_mask_across = (dataset_across.trialmask.trial_mask(fold_mask_across));
    target_mask_across =  (dataset_across.targetmask.target_mask(fold_mask_across));
    phase_mask_across =  (dataset_across.trialmask.phase_mask(fold_mask_across));
    test_x_across = zvalue_across(fold_mask_across,1:end-var_num);
    test_y_across = data_curr_across(fold_mask_across,end-var_num+1:end);
    
    test_len_across = size(test_x_across,1);
    b_predict_across = repmat(beta(1,:),[test_len_across,1]);
    k_predict_across = beta(2:end,:);
    y_predict_across = b_predict_across +  test_x_across*k_predict_across;
    
    y_predict_across = smoothdata(y_predict_across,smooth_sigma,'gaussian',smooth_win);
    
    [R_sqr_across,MSE_across,CC_across] = RegressionFuncEvaluation(y_num,test_len_across,y_predict_across,test_y_across);
    save([save_dir 'Train' hand '.HandAcross.' across_id 'Fold' num2str(fold_id) '.mat'],'y_predict_across','test_y_across','R_sqr_across','MSE_across','CC_across','trial_mask_across','target_mask_across','phase_mask_across','test_x_across','b_predict_across','k_predict_across');
end
