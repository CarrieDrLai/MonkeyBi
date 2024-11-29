function [R_square_list,MSE_list,Coef_list] = RegressionPredArea(save_pls_decode_dir,dataset_L,dataset_R,dataset_Bi,NeuronMask)

    R_square_list =zeros([size(NeuronMask,1) 9 4]); % Area * (3trainset * 3testset(within&across)) * 4var 
    MSE_list = zeros([size(NeuronMask,1) 9 4]); % Area * (3trainset * 3testset(within&across)) * 4var 
    Coef_list = zeros([size(NeuronMask,1) 9 4]); % Area * (3trainset * 3testset(within&across)) * 4var 
    
    for n_id = 1: size(NeuronMask,1)
        R_sqr_nid = zeros([9,4,5]);
        MSE_nid = zeros([9,4,5]);
        CC_nid = zeros([9,4,5]);
        currid = 
        for fold_id=1:5
            [R_sqr_L,MSE_L,CC_L] = RegressionFunc(save_pls_decode_dir,dataset_L,'L',fold_id,currid,dataset_R,dataset_Bi);
            [R_sqr_R,MSE_R,CC_R] = RegressionFunc(save_pls_decode_dir,dataset_R,'R',fold_id,currid,dataset_L,dataset_Bi);
            [R_sqr_Bi,MSE_Bi,CC_Bi] = RegressionFunc(save_pls_decode_dir,dataset_Bi,'B',fold_id,currid,dataset_L,dataset_R);
            R_sqr_nid(:,:,fold_id) = [R_sqr_L;R_sqr_R;R_sqr_Bi];
            MSE_nid(:,:,fold_id) = [MSE_L;MSE_R;MSE_Bi];
            CC_nid(:,:,fold_id) = [CC_L;CC_R;CC_Bi];
        end
        R_sqr_nid_mean = mean(R_sqr_nid,3);
        MSE_nid_mean = mean(MSE_nid,3);
        CC_nid_mean = mean(CC_nid,3);
        R_square_list(n_id,:,:) = R_sqr_nid_mean;
        MSE_list(n_id,:,:) = MSE_nid_mean;
        Coef_list(n_id,:,:) = CC_nid_mean;           
        save([save_pls_decode_dir 'Nid' num2str(n_id) '.AllEvaluate.mat'],'R_sqr_nid_mean','MSE_nid_mean','CC_nid_mean','R_sqr_nid','MSE_nid','CC_nid');
    end
    save([save_pls_decode_dir 'AllEvaluate.mat'],'R_square_list','MSE_list','Coef_list');
end