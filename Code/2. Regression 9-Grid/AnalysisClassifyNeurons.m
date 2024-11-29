function [type_list]=AnalysisClassifyNeurons(savedir,area,cc_thresh);
    for a = 1:size(area,1)
        load([savedir 'AnalysisDecode\' 'Units_cc_' area(a,:)  '.mat']);
        type_list = zeros([size(date_mask,1),8]); % (Uni-Bi, Bi-Uni); *4
        for var_i = 1:4
            neuron_uni_within_only = Across_Uni_2_Bi(:,2*var_i-1)>cc_thresh;
            neuron_bi_across_only = Across_Uni_2_Bi(:,2*var_i)>cc_thresh;
            type_list(:,2*var_i-1) = neuron_bi_across_only * 2 + neuron_uni_within_only;

            neuron_bi_within_only = Across_Bi_2_Uni(:,2*var_i-1)>cc_thresh;
            neuron_uni_across_only = Across_Bi_2_Uni(:,2*var_i)>cc_thresh;
            type_list(:,2*var_i) = 4 + neuron_uni_across_only * 2 + neuron_bi_within_only;
        end
        
        type_neuron_num = zeros([4,8]);
        for var_i = 1:4
            temp = [];
            for i = 1:4
                temp=[temp;sum(type_list(:,var_i*2-1)==i-1);];
            end
            type_neuron_num(:,var_i*2-1) = temp;
        end
        for var_i = 1:4
            temp = [];
            for i = 1:4
                temp=[temp;sum(type_list(:,var_i*2)==4+i-1);];
            end
            type_neuron_num(:,var_i*2) = temp;
        end        
        save([savedir 'AnalysisDecode\' 'classfy_neuron_type_' area(a,:) '_thresh_' num2str(cc_thresh) '.mat' ],'type_list','cc_thresh','type_neuron_num');
    end
end