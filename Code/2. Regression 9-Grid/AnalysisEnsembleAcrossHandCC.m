function [Uni2Bi_Lxy,Uni2Bi_Rxy,Bi2Uni_LxyRxy,Uni2Bi_Lxy_shuffle,Uni2Bi_Rxy_shuffle,Bi2Uni_LxyRxy_shuffle,p_list,p_hat_list,h_list]=AnalysisEnsembleAcrossHandCC(savedir,date_list)
    
    Uni2Bi_Lxy = zeros([1,2]);
    Uni2Bi_Rxy = zeros([1,2]);
    Bi2Uni_LxyRxy  = zeros([1,4]);
    
    Uni2Bi_Lxy_shuffle = zeros([1,2]);
    Uni2Bi_Rxy_shuffle = zeros([1,2]);
    Bi2Uni_LxyRxy_shuffle  = zeros([1,4]);

    for d_i = 1:size(date_list,1)
        load([savedir date_list(d_i,:) '\pls_decode' '\SDJudge.AllEvaluate.mat']);
        Coef_list = reshape(Coef_list,[size(Coef_list,2) size(Coef_list,3)]);
        Uni2Bi_Lxy(d_i,:) = Coef_list(2,1:2);
        Uni2Bi_Rxy(d_i,:) = Coef_list(5,3:4);
        Bi2Uni_LxyRxy(d_i,1:2) = Coef_list(8,1:2);
        Bi2Uni_LxyRxy(d_i,3:4) = Coef_list(9,3:4);

        load([savedir date_list(d_i,:) '\pls_decode' '\SDJudgeShuffle.AllEvaluate.mat']);
        Coef_list = reshape(Coef_list,[size(Coef_list,2) size(Coef_list,3)]);
        Uni2Bi_Lxy_shuffle(d_i,:) = Coef_list(2,1:2);
        Uni2Bi_Rxy_shuffle(d_i,:) = Coef_list(5,3:4);
        Bi2Uni_LxyRxy_shuffle(d_i,1:2) = Coef_list(8,1:2);
        Bi2Uni_LxyRxy_shuffle(d_i,3:4) = Coef_list(9,3:4);
    end
    
    [h_Lx,p_Lx] = ttest(Uni2Bi_Lxy(:,1),Uni2Bi_Lxy_shuffle(:,1));
    [h_Ly,p_Ly] = ttest(Uni2Bi_Lxy(:,2),Uni2Bi_Lxy_shuffle(:,2));
    [h_Rx,p_Rx] = ttest(Uni2Bi_Rxy(:,1),Uni2Bi_Rxy_shuffle(:,1));
    [h_Ry,p_Ry] = ttest(Uni2Bi_Rxy(:,2),Uni2Bi_Rxy_shuffle(:,2));
    [h_BLx,p_BLx] = ttest(Bi2Uni_LxyRxy(:,1),Bi2Uni_LxyRxy_shuffle(:,1));
    [h_BLy,p_BLy] = ttest(Bi2Uni_LxyRxy(:,2),Bi2Uni_LxyRxy_shuffle(:,2));
    [h_BRx,p_BRx] = ttest(Bi2Uni_LxyRxy(:,3),Bi2Uni_LxyRxy_shuffle(:,3));
    [h_BRy,p_BRy] = ttest(Bi2Uni_LxyRxy(:,4),Bi2Uni_LxyRxy_shuffle(:,4));
    
    test_data = [Uni2Bi_Lxy,Uni2Bi_Rxy,Bi2Uni_LxyRxy];
    test_data_shuffle = [Uni2Bi_Lxy_shuffle,Uni2Bi_Rxy_shuffle,Bi2Uni_LxyRxy_shuffle];
    h_list = zeros([1]);
    p_list = zeros([1]);    
    p_hat_list = zeros([1]);
    for i = 1:8
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
    save([savedir 'AnalysisDecode\' 'EnsembleAcrosssHand_cc.mat'],'date_list','Uni2Bi_Lxy','Uni2Bi_Rxy','Bi2Uni_LxyRxy','Uni2Bi_Lxy_shuffle','Uni2Bi_Rxy_shuffle','Bi2Uni_LxyRxy_shuffle');
    save([savedir 'AnalysisDecode\' 'EnsembleAcrosssHand_cc_pvalue.mat'],'p_list','p_hat_list','h_list');
end