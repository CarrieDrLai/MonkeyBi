function [trial_list,trial_no,phase_no,spikes_50Hz,start_point,start_point_30k] = InitializationcompileNEVNS6file(NEV,readNeurons,nsp_ch,start_prepare_id,f,f0,f_joystick)
    %%%%%% Basic Info
    %%%%%% Analog 

    %%%%%% Neural 
    timeStamps = NEV.Data.Spikes.TimeStamp;
    numOfTimestamps = timeStamps(end);
    %f = NEV.MetaTags.SampleRes;
    electrode = NEV.Data.Spikes.Electrode;
    neuron_id = (electrode - 1) * 6 + uint16(NEV.Data.Spikes.Unit);
%     readNeurons = unique(neuron_id);
%     temp = readNeurons;
%     readNeurons=[];
%     for i = 1:size(temp,2)
%         if mod(temp(i),6)>0
%             readNeurons=[readNeurons,temp(i)];
%         end
%     end
% 
%     if size(NEV.ElectrodesInfo,2)<256
%         nsp_ch = 128;
%         readNeurons = readNeurons(find(readNeurons<=128*6));
%     else
%         nsp_ch = 256;
%         readNeurons = readNeurons(find(readNeurons<=256*6));
%     end
% 
%     del_nid_mask = mod(readNeurons,6);
%     del_nid = find(del_nid_mask>1);
%     del_nid = del_nid((readNeurons(del_nid) - readNeurons(del_nid-1))~=1);
%     n_id = 1:1:size(readNeurons,2);
%     n_id = n_id(~ismember(n_id,del_nid));
%     readNeurons = readNeurons(n_id);
%     readElectrodes = unique(ceil(double(readNeurons)/6));
%     %notmember_list = [187,259,313,379,385,433,445,446,547,571,583,637,643,751,757,758,763,811,812,813,817,937,985,1027,1028,1029,1045,1051,1052,1053,1069,1070];
%     notmember_list=[];
%     %notmember_list =[561,573,1047,1071];
%     %notmember_list =[765,639]; % 20221107
%     %notmember_list =[735]; % 20230313
%     %floor(1+double(readNeurons)/6)
%     readNeurons = readNeurons(~ismember(readNeurons,notmember_list));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%   1 Com Marker   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    joystick_timestamp = NEV.Data.SerialDigitalIO.TimeStamp; %30k
    joystick_timestamp = joystick_timestamp./(30000/f);      %downsample
    joystick_second = NEV.Data.SerialDigitalIO.TimeStampSec;
    joystick_comwrite = NEV.Data.SerialDigitalIO.UnparsedData';

    % 99: Go out & Hold Judge; 100: Reward 102:Return center & Hold Judge 104: Reward 105:Intertrial 
    phase1 = find(joystick_comwrite==99); 
    phase1 = phase1(start_prepare_id:end);
    phase_end = find(joystick_comwrite==max(joystick_comwrite));% intertrial
    phase_end = phase_end(phase_end>phase1(1));

    if phase1(start_prepare_id)<4
        start_point = joystick_timestamp(phase1(1)-2);     %f
        start_com_id = phase1(1)-2;
    else
        start_point = joystick_timestamp(phase1(2)-2);  
        start_com_id = phase1(2)-2;
    end

    start_trial_no = joystick_comwrite(min(find(joystick_comwrite==10))-1)-48;
    
    trial=[];
    trial_com_id=[];
    for i = start_trial_no:max(size(phase1,2),size(phase_end,2))
        A = double(floor(double(i)/100));
        B = double(floor((double(i) - 100*A)/10));
        C = double(i - 100 * A - 10 * B);
        if A+B==0 % C = 1:9
            temp = find(joystick_comwrite==C+48);
            if size(temp,2)>0
                trial_com_id = [trial_com_id,double(min(temp))];
            end
        else %A+B>0
            if A==0 
                temp = find(joystick_comwrite==B+48);
                if size(temp,2)>0
                    if temp(end)+1>size(joystick_comwrite,2)
                        temp = temp(1:end-1);
                    end
                    if size(temp,2)>0 
                        temp_id = find(joystick_comwrite(temp+1)==C+48);
                        trial_com_id = [trial_com_id,double(temp(min(temp_id)))];
                    end
                end
            else
                temp = find(joystick_comwrite==A+48);
                if size(temp,2)>0
                    if temp(end)+2>size(joystick_comwrite,2)
                        temp = temp(1:end-1);
                    end
                    if size(temp,2)>0
                        temp_id = (joystick_comwrite(temp+1)==B+48);
                        temp_id2 = (joystick_comwrite(temp+2)==C+48);
                        temp2 = temp(temp_id & temp_id2);
                        trial_com_id = [trial_com_id,double(min(temp2))];
                    end
                end
            end
        end
        trial = [trial,double(100 * A + 10 * B + C)];
    end
       
    trial_id = trial_com_id(trial_com_id>=start_com_id);
    %end_point = joystick_timestamp(temp2(end));
    trial_no_timestamp = joystick_timestamp(trial_id);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%   1 Com Marker   %%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%   2 Trial & Phase  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    trial_no_50Hz = (trial_no_timestamp-start_point+1)/(f/f0)+1;
    phase_time_50Hz = (joystick_timestamp-start_point+1)/(f/f0)+1;
    %trial_list =[];
    % k=1;
    %for i = 1:(size(trial_no_50Hz,2))
    %    while rem(k,10)~= mod(trial(i),10)
    %        k=k+1;
    %    end
    %    trial_list = [trial_list;k];
    %end
    %trial_list = 1:1:max(trial);
    % 99: Go out & Hold Judge; 100: Reward 102:Return center & Hold Judge 104: Reward 105:Intertrial 
    phase_reward_center = find(joystick_comwrite==100);% reward
    phase_reward_center = phase_reward_center(phase_reward_center>phase1(1) & phase_reward_center< phase_end(end));
    phase_go = find(joystick_comwrite==102);% phase go =102 
    phase_go = phase_go(phase_go>phase1(1) & phase_go< phase_end(end));
    phase_reward = find(joystick_comwrite==104);% reward
    phase_reward = phase_reward(phase_reward>phase1(1) & phase_reward< phase_end(end));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%²¹³äÈ±Ê§trial
    %trial_id_lose = find(diff(trial)>1 & diff(trial)<10)+1;
    %trial_id_lose = find(diff(trial_com_id)>40);
    %if size(trial_id_lose,2)>0
    %    for i = 1:size(trial_id_lose,2)
    %        lose_next_pre = phase1(find(trial_com_id(trial_id_lose(i))<find(joystick_comwrite==99) & find(joystick_comwrite==99) < trial_com_id(trial_id_lose(i)+1)));
    %        lose_next_pre = lose_next_pre(end);
    %        lose_last_inter = phase_end(find(trial_com_id(trial_id_lose(i))<find(joystick_comwrite==104) & find(joystick_comwrite==104) < trial_com_id(trial_id_lose(i)+1)));
    %        lose_last_inter = lose_last_inter(1);
    %        if size(lose_next_pre,1) == 1
    %            trial_no_50Hz = [trial_no_50Hz,phase_time_50Hz(lose_next_pre)];
    %        else
    %            if size(lose_last_inter,1) == 1
    %                trial_no_50Hz = [trial_no_50Hz,phase_time_50Hz(lose_last_inter)];
    %            end
    %        end
    %    end
    %end
    trial_no_50Hz = sort(trial_no_50Hz);
    trial_list = 1:1:size(trial_no_50Hz,2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%²¹³äÈ±Ê§trial

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%trialno mask & phaseno mask
    % 99: Go out & Hold Judge; 100: Reward 102:Return center & Hold Judge 104: Reward 105:Intertrial 
    phase_no = zeros(1,trial_no_50Hz(end));
    %phase_no(1,phase_time_50Hz(phase1))=1;%del %%prepare_center 
    %phase_no(1,phase_time_50Hz(phase_pre_joystick))=2;%del  %%prepare_joystick
    phase_no(1,phase_time_50Hz(phase1))=3;%99  %%hold_center
    phase_no(1,phase_time_50Hz(phase_reward_center))=4;%  %%reward_center = 100
    %phase_no(1,phase_time_50Hz(phase_pre))=5;           
    phase_no(1,phase_time_50Hz(phase_go))=6;   %%go = 102 
    phase_no(1,phase_time_50Hz(phase_reward))=7;%104
    phase_no(1,phase_time_50Hz(phase_end))=8;%105
    phase_id = find(phase_no>0);
    for i = 1:(size(phase_id,2)-1)
        phase_no(1,phase_id(i):(phase_id(i+1)-1))=phase_no(phase_id(i));
    end
    phase_no = phase_no(1,1:end-1);

    trial_no = zeros(1,phase_id(end)-1);
    for i = 1:(size(trial_no_50Hz,2)-1)
        trial_no(1,trial_no_50Hz(i):(trial_no_50Hz(i+1)-1))=trial_list(i);
        %timestamp_pre_joystick = joystick_comwrite(phase_pre_joystick(i));
    end
    trial_no(1,trial_no_50Hz(size(trial_no_50Hz,2)):end) = trial_list(end);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%   2 Trial & Phase  %%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%  3  Single Unit Spike  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numOfNeurons = size(readNeurons,2);
    spikes_50Hz = [];
    start_point_30k = start_point/(f/30000);
    for n_id = 1: numOfNeurons
        spike0 = zeros(numOfTimestamps, 1);
        spike0(timeStamps(neuron_id == readNeurons(n_id)),1) = 1;
        spike_nid_50Hz=[];
        for t = 1:(size(trial_no,2))
            spike_nid_50Hz = [spike_nid_50Hz,sum(spike0((start_point_30k+(t-1)*(30000/f0)+1):(start_point_30k+t*30000/f0),:))];
        end
        %spike_nid_50Hz = [spike_nid_50Hz,sum(spike0((start_point_30k+(size(trial_no,2)-1)*(30000/f0)+1):(start_point_30k+size(trial_no,2)*(30000/f0)),:))];
        spikes_50Hz = [spikes_50Hz;spike_nid_50Hz];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%  3  Single Unit Spike  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%  4  Single Unit Waveform  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %waveforms_50Hz = zeros(numOfNeurons,1,size(spike_nid_50Hz,2));
    %waveforms_50Hz = zeros(numOfNeurons,46,size(spike_nid_50Hz,2));
    %start_point_30k = start_point/(f/30000);
%     for t = 1:(size(trial_no,2)-1)
%         id1 = start_point_30k+(t-1)*(30000/f0)+1;
%         id2 = (start_point_30k+t*30000/f0);
%         time_id = find(timeStamps>id1 & timeStamps<id2);
%         nid_curr_id = neuron_id(time_id);
%         nid_curr = unique(nid_curr_id);
%         if size(nid_curr,2)>0
%             for i = 1:size(nid_curr,2)
%                 n = find(readNeurons == nid_curr(i));
%                 if size(n,2)>0
%                     waveforms0 = mean(NEV.Data.Spikes.Waveform(:,time_id(nid_curr_id==nid_curr(i))),2);
%                     waveforms_50Hz(n,:,t) = waveforms0;
%                 end
%             end
%         end
%     end
    % temp = waveforms_50Hz(13,:,:);
    % plot(reshape(temp,[46,size(temp,3)]));
    clear NEV;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%  4  Single Unit Waveform  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Combine NEV data
end