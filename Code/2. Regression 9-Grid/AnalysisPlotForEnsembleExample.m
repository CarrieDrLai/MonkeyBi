function [temp_within,temp_across,cc_temp_within,cc_temp_across] = AnalysisEnsembleWithinHandCC(savedir)

    %% plot ensemble decode data for example

    savedir =  [data_dir_root 'Decode\'];

    %W-BR  A-R
    load([savedir '2024' '-' '03' '-' '25' '\pls_decode\SDJudge\' 'TrainB.HandWithin.Fold5.mat'])
    BRx=[test_y(:,3),y_predict(:,3)];
    BRy=[test_y(:,4),y_predict(:,4)];
    load([savedir '2024' '-' '03' '-' '25' '\pls_decode\SDJudge\' 'TrainB.HandAcross.2Fold5.mat'])
    Rx=[test_y_across(:,3),y_predict_across(:,3)];
    Ry=[test_y_across(:,4),y_predict_across(:,4)];

    temp_within = BRy(1:300,:);
    cc_temp_within = corr(temp_within(:,1),temp_within(:,2))
    temp_across = Ry(300:599,:);
    cc_temp_across = corr(temp_across(:,1),temp_across(:,2))
end