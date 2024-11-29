function [Within,Across_Uni_2_Bi,Across_Bi_2_Uni,aa,p_hat_list]=AnalysisUnitsCC(savedir,date_list,area,area_ch)
    
    for a = 1:size(area,1)
        Within = zeros([1,8]);          % WLx-WBLx; WLy-WBLy: WRx-WBRx; WRy-WBRy:
        Across_Uni_2_Bi = zeros([1,8]); % WLx-ABLx; WLy-ABLy: WRx-ABRx; WRy-ABRy:
        Across_Bi_2_Uni = zeros([1,8]); % WBLx-ALx; WBLy-ALy: WBRx-ARx; WBRy-ARy:     
        N_num = 0;
        date_mask = [];
        p_within_list = zeros([1,4]);
        h_within_list = zeros([1,4]);
        p_Across_Uni_list = zeros([1,4]);
        h_Across_Uni_list = zeros([1,4]);
        p_Across_Bi_list = zeros([1,4]);
        h_Across_Bi_list = zeros([1,4]);
        
        for d_i = 1:size(date_list,1)
            load([savedir date_list(d_i,:) '\pls_decode' '\SDJudge\NeuronMask_SDJudge.mat']);
            mask = find(NeuronNo_SDJudge <= 6*sum(area_ch(1:a+1)) & NeuronNo_SDJudge >= 1+6*sum(area_ch(1:a))) ;
            N_num = N_num+size(mask,2);
            
            WLx_WBLx = reshape(Coef_SDJudge_units(mask,[1,7],1),[size(mask,2) 2]);
            WRx_WBRx = reshape(Coef_SDJudge_units(mask,[4,7],3),[size(mask,2) 2]);
            WLy_WBLy = reshape(Coef_SDJudge_units(mask,[1,7],2),[size(mask,2) 2]);
            WRy_WBRy = reshape(Coef_SDJudge_units(mask,[4,7],4),[size(mask,2) 2]);
            
            Within(N_num-size(mask,2)+1:N_num,1:2) = WLx_WBLx;
            Within(N_num-size(mask,2)+1:N_num,5:6) = WRx_WBRx;
            Within(N_num-size(mask,2)+1:N_num,3:4) = WLy_WBLy;
            Within(N_num-size(mask,2)+1:N_num,7:8) = WRy_WBRy;
           
            WLx_ABLx = reshape(Coef_SDJudge_units(mask,[1,3],1),[size(mask,2) 2]);
            WRx_ABRx = reshape(Coef_SDJudge_units(mask,[4,6],3),[size(mask,2) 2]);
            WLy_ABLy = reshape(Coef_SDJudge_units(mask,[1,3],2),[size(mask,2) 2]);
            WRy_ABRy = reshape(Coef_SDJudge_units(mask,[4,6],4),[size(mask,2) 2]);
            
            Across_Uni_2_Bi(N_num-size(mask,2)+1:N_num,1:2) = WLx_ABLx;
            Across_Uni_2_Bi(N_num-size(mask,2)+1:N_num,5:6) = WRx_ABRx;
            Across_Uni_2_Bi(N_num-size(mask,2)+1:N_num,3:4) = WLy_ABLy;
            Across_Uni_2_Bi(N_num-size(mask,2)+1:N_num,7:8) = WRy_ABRy;     
            
            WBLx_ALx = reshape(Coef_SDJudge_units(mask,[7,8],1),[size(mask,2) 2]);
            WBRx_ARx = reshape(Coef_SDJudge_units(mask,[7,9],3),[size(mask,2) 2]);
            WBLy_ALy = reshape(Coef_SDJudge_units(mask,[7,8],2),[size(mask,2) 2]);
            WBRy_ARy = reshape(Coef_SDJudge_units(mask,[7,9],4),[size(mask,2) 2]);
            
            Across_Bi_2_Uni(N_num-size(mask,2)+1:N_num,1:2) = WBLx_ALx;
            Across_Bi_2_Uni(N_num-size(mask,2)+1:N_num,5:6) = WBRx_ARx;
            Across_Bi_2_Uni(N_num-size(mask,2)+1:N_num,3:4) = WBLy_ALy;
            Across_Bi_2_Uni(N_num-size(mask,2)+1:N_num,7:8) = WBRy_ARy;      
            date_mask = [date_mask;repmat(d_i,[size(mask,2) 1])];
        end
        if ~exist([savedir 'AnalysisDecode\' ], 'dir')
            mkdir([savedir 'AnalysisDecode\']);
        end
        
        cumul_Within = [];
        cumul_U2B = [];
        cumul_B2U = [];
        for i = -0.8:0.1:0.7
            cumul_Within = [cumul_Within;sum((Within<=i+0.1))];
            cumul_U2B = [cumul_U2B;sum((Across_Uni_2_Bi<=i+0.1))];
            cumul_B2U = [cumul_B2U;sum((Across_Bi_2_Uni<=i+0.1))];
        end
        cumul_Within = cumul_Within./size(Within,1);
        cumul_U2B = cumul_U2B./size(Across_Uni_2_Bi,1);
        cumul_B2U = cumul_B2U./size(Across_Bi_2_Uni,1);
        cumul_xaxis = -0.8:0.1:0.7;
        save([savedir 'AnalysisDecode\' 'Units_cc_' area(a,:) '.mat'],'date_list','Within','Across_Uni_2_Bi','Across_Bi_2_Uni','cumul_xaxis','cumul_Within','cumul_U2B','cumul_B2U','date_mask');
        
        p_hat_list = zeros([4,3]);
        for i = 1:4
            [p_within_list(i),h_within_list(i)] = signrank(Within(:,2*i-1),Within(:,2*i),0.05);
            [p_Across_Uni_list(i),h_Across_Uni_list(i)] = signrank(Across_Uni_2_Bi(:,2*i-1),Across_Uni_2_Bi(:,2*i),0.05);
            [p_Across_Bi_list(i),h_Across_Bi_list(i)] = signrank(Across_Bi_2_Uni(:,2*i-1),Across_Bi_2_Uni(:,2*i),0.05);
            
            aa=[p_within_list(i),p_Across_Uni_list(i),p_Across_Bi_list(i)];
            for aaa = 1:3
                if aa(aaa)<0.001
                    p_hat_list(i,aaa) = 3;
                else
                    if aa(aaa)<0.01
                        p_hat_list(i,aaa) = 2;
                    else
                        if aa(aaa)<0.05
                            p_hat_list(i,aaa) = 1;
                        end
                    end
                end
            end
        end
        save([savedir 'AnalysisDecode\' 'Units_cc_pvalue_' area(a,:) '.mat'],'p_within_list','p_Across_Uni_list','p_Across_Bi_list','h_within_list','h_Across_Uni_list','h_Across_Bi_list','p_hat_list');  
    end
end