function [R_square_list,MSE_list,Coef_list] = RegressionPred(save_pls_decode_dir,dataset_L,dataset_R,dataset_Bi,NeuronMask,keyword,smooth_sigma,smooth_win)

    var_num = size(dataset_L.data_Go.vel_Go,2) + size(dataset_L.data_Go.traj_Go,2);
    R_square_list =zeros([size(NeuronMask,2) 9 var_num]); % Neuron\Area * (3trainset * 3testset(within&across)) * 4var 
    MSE_list = zeros([size(NeuronMask,2) 9 var_num]); % Neuron\Area * (3trainset * 3testset(within&across)) * 4var 
    Coef_list = zeros([size(NeuronMask,2) 9 var_num]); % Neuron\Area * (3trainset * 3testset(within&across)) * 4var
    
    A = ['LM1';'PMd';'RM1';];
    fold_no = max(dataset_L.foldmask.fold_mask);
    for n_id = 1: size(NeuronMask,2)
        R_sqr_nid = zeros([9,var_num,fold_no]);
        MSE_nid = zeros([9,var_num,fold_no]);
        CC_nid = zeros([9,var_num,fold_no]);
        if keyword(1:4) == 'unit'
            current_id = n_id;
            savedirkey = [keyword '.' 'Nid' num2str(n_id)];
        else
            if keyword(1:4) == 'area'
                current_id = find(NeuronMask(:,n_id)==1);
                savedirkey = [ keyword '.' A(n_id,:) ];
                if size(current_id,1) ==0
                    continue
                end
            else
                current_id = find(NeuronMask(:,n_id)==1);
                savedirkey = keyword;
            end
        end
%         if length(keyword)>10
%             savedirkey = [savedirkey '.Shuffle']
%         end
        
        
        for fold_id=1:fold_no
            [R_sqr_L,MSE_L,CC_L] = RegressionFunc(save_pls_decode_dir,dataset_L,'L',fold_id,current_id,dataset_R,dataset_Bi,savedirkey,var_num,smooth_sigma,smooth_win);
            [R_sqr_R,MSE_R,CC_R] = RegressionFunc(save_pls_decode_dir,dataset_R,'R',fold_id,current_id,dataset_L,dataset_Bi,savedirkey,var_num,smooth_sigma,smooth_win);
            [R_sqr_Bi,MSE_Bi,CC_Bi] = RegressionFunc(save_pls_decode_dir,dataset_Bi,'B',fold_id,current_id,dataset_L,dataset_R,savedirkey,var_num,smooth_sigma,smooth_win);
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
        save([save_pls_decode_dir savedirkey '.AllEvaluate.mat'],'R_sqr_nid_mean','MSE_nid_mean','CC_nid_mean','R_sqr_nid','MSE_nid','CC_nid');
    end
    save([save_pls_decode_dir keyword '.AllEvaluate.mat'],'R_square_list','MSE_list','Coef_list');
end