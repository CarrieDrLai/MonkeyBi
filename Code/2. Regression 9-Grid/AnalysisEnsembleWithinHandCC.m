function [left_xy,right_xy,bi_LxyRxy,left_xy_shuffle,right_xy_shuffle,bi_LxyRxy_shuffle,p_list,p_hat_list,h_list]=AnalysisEnsembleWithinHandCC(savedir,date_list)
    
    left_xy = zeros([1,2]);
    right_xy = zeros([1,2]);
    bi_LxyRxy  = zeros([1,4]);
    left_xy_shuffle = zeros([1,2]);
    right_xy_shuffle = zeros([1,2]);
    bi_LxyRxy_shuffle  = zeros([1,4]);

    for d_i = 1:size(date_list,1)
        load([savedir date_list(d_i,:) '\pls_decode' '\SDJudge.AllEvaluate.mat']);
        Coef_list = reshape(Coef_list,[size(Coef_list,2) size(Coef_list,3)]);
        left_xy(d_i,:) = Coef_list(1,1:2);
        right_xy(d_i,:) = Coef_list(4,3:4);
        bi_LxyRxy(d_i,:) = Coef_list(7,1:4);

        load([savedir date_list(d_i,:) '\pls_decode' '\SDJudgeShuffle.AllEvaluate.mat']);
        Coef_list = reshape(Coef_list,[size(Coef_list,2) size(Coef_list,3)]);
        left_xy_shuffle(d_i,:) = Coef_list(1,1:2);
        right_xy_shuffle(d_i,:) = Coef_list(4,3:4);
        bi_LxyRxy_shuffle(d_i,:) = Coef_list(7,1:4);
    end
    
%     [h_Lx,p_Lx] = signrank(left_xy(:,1),left_xy_shuffle(:,1));
%     [h_Ly,p_Ly] = signrank(left_xy(:,2),left_xy_shuffle(:,2));
%     [h_Rx,p_Rx] = signrank(right_xy(:,1),right_xy_shuffle(:,1));
%     [h_Ry,p_Ry] = signrank(right_xy(:,2),right_xy_shuffle(:,2));
%     [h_BLx,p_BLx] = signrank(bi_LxyRxy(:,1),bi_LxyRxy_shuffle(:,1));
%     [h_BLy,p_BLy] = signrank(bi_LxyRxy(:,2),bi_LxyRxy_shuffle(:,2));
%     [h_BRx,p_BRx] = signrank(bi_LxyRxy(:,3),bi_LxyRxy_shuffle(:,3));
%     [h_BRy,p_BRy] = signrank(bi_LxyRxy(:,4),bi_LxyRxy_shuffle(:,4));
    
    test_data = [left_xy,right_xy,bi_LxyRxy];
    test_data_shuffle = [left_xy_shuffle,right_xy_shuffle,bi_LxyRxy_shuffle];
    h_list = zeros([1]);
    p_list = zeros([1]);    
    p_hat_list = zeros([1]);
    for i = 1:8
        %[h_list(i),p_list(i)] = ttest(test_data(:,i),test_data_shuffle(:,i))
        [p_list(i),h_list(i)] = signrank(test_data(:,i),test_data_shuffle(:,i),0.05);
        if p_list(i)<0.001
            p_hat_list(i) = 3;
        else
            if p_list(i)<0.01
                p_hat_list(i) = 2;
            else
                if p_list(i)<0.05
                    p_hat_list(i) = 1;
                end
            end
        end
    end
    
    if ~exist([savedir 'AnalysisDecode\' ], 'dir')
        mkdir([savedir 'AnalysisDecode\']);
    end
    save([savedir 'AnalysisDecode\' 'EnsembleWithinHand_cc.mat'],'date_list','left_xy','right_xy','bi_LxyRxy','left_xy_shuffle','right_xy_shuffle','bi_LxyRxy_shuffle');
    save([savedir 'AnalysisDecode\' 'EnsembleWithinHand_cc_pvalue.mat'],'p_list','p_hat_list','h_list');
end