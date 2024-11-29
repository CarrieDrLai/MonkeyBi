function [readNeurons,nsp_ch] = InitializationReadNeuron(NEV,notmember_list)

    electrode = NEV.Data.Spikes.Electrode;
    neuron_id_mask =  (NEV.Data.Spikes.Unit<7) & (NEV.Data.Spikes.Unit>0);
    neuron_id0 =  NEV.Data.Spikes.Unit(neuron_id_mask);
    electrode0 = NEV.Data.Spikes.Electrode(neuron_id_mask);
    neuron_id = (electrode0 - 1) * 6 + uint16(neuron_id0);
    readNeurons = unique(neuron_id);
%     temp = readNeurons;
%     readNeurons=[];
%     for i = 1:size(temp,2)
%         if mod(temp(i),6)>0
%             readNeurons=[readNeurons,temp(i)];
%         end
%     end

    if size(NEV.ElectrodesInfo,2)<256
        nsp_ch = 128;
        readNeurons = readNeurons(find(readNeurons<=128*6));
    else
        nsp_ch = 256;
        readNeurons = readNeurons(find(readNeurons<=256*6));
    end
% 
%     del_nid_mask = mod(readNeurons,6);
%     del_nid = find(del_nid_mask>1);
%     del_nid = del_nid((readNeurons(del_nid) - readNeurons(del_nid-1))~=1);
%     n_id = 1:1:size(readNeurons,2);
%     n_id = n_id(~ismember(n_id,del_nid));
%     readNeurons = readNeurons(n_id);
%     readElectrodes = unique(ceil(double(readNeurons)/6));

    readNeurons = readNeurons(~ismember(readNeurons,notmember_list));

end