function [trial_list,trial_no,phase_no,spikes_50Hz,readNeurons,nsp_ch,start_point] = compileNEVNS6file(NEV,start_prepare_id,f0,f_joystick)
    %%%%%% Basic Info
    %%%%%% Analog 


    %%%%%% Neural 
    timeStamps = NEV.Data.Spikes.TimeStamp;
    numOfTimestamps = timeStamps(end);
    f = NEV.MetaTags.SampleRes;
    electrode = NEV.Data.Spikes.Electrode;
    neuron_id = (electrode - 1) * 6 + uint16(NEV.Data.Spikes.Unit);
    readNeurons = unique(neuron_id);
    temp = readNeurons;
    readNeurons=[];
    for i = 1:size(temp,2)
        if mod(temp(i),6)>0
            readNeurons=[readNeurons,temp(i)];
        end
    end

    if size(NEV.ElectrodesInfo,2)<256
        nsp_ch = 128;
        readNeurons = readNeurons(find(readNeurons<=128*6));
    else
        nsp_ch = 256;
        readNeurons = readNeurons(find(readNeurons<=256*6));
    end

    del_nid_mask = mod(readNeurons,6);
    del_nid = find(del_nid_mask>1);
    del_nid = del_nid((readNeurons(del_nid) - readNeurons(del_nid-1))~=1);
    n_id = 1:1:size(readNeurons,2);
    n_id = n_id(~ismember(n_id,del_nid));
    readNeurons = readNeurons(n_id);
    readElectrodes = unique(ceil(double(readNeurons)/6));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%   1 Com Marker   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    joystick_timestamp = NEV.Data.SerialDigitalIO.TimeStamp;
    joystick_second = NEV.Data.SerialDigitalIO.TimeStampSec;
    joystick_comwrite = NEV.Data.SerialDigitalIO.UnparsedData';

    % 99: Go out & Hold Judge; 100: Reward 102:Return center & Hold Judge 104: Reward 105:Intertrial 
    phase1 = find(joystick_comwrite==99); 
    phase1 = phase1(start_prepare_id:end-1);
    phase_end = find(joystick_comwrite==max(joystick_comwrite));% intertrial
    phase_end = phase_end(phase_end>phase1(1));

    if phase1(start_prepare_id)>2
        start_point = joystick_timestamp(phase1(1)-2);
        start_com_id = phase1(1)-2;
    else
        start_point = joystick_timestamp(phase1(2)-2);
        start_com_id = phase1(2)-2;
    end

    no = find(joystick_comwrite==10);
    trial=[];
    trial_com_id=[];
    for i = 1:size(no,2)
        id_curr = 0;
        if joystick_comwrite(no(i)-1)<99
            if i == 1
                length = no(i)-1;
            else
                length = (no(i)-1) - no(i-1);
            end
            if length>0
                for j = 1:length
                    id_curr = id_curr + (joystick_comwrite(no(i)-length-1+j) - 48) * 10^(length-j);
                end
                trial = [trial,uint32(id_curr)];
                %trial_no_timestamp = [trial_no_timestamp,joystick_timestamp(no(i)-length)];
                trial_com_id = [trial_com_id,no(i)-length];
            end
        end
    end % 注意 1108session1,不停发送trialno & phase % 1109之后，每个phases发一次trial
    if (unique(uint32(trial)))<size(trial,2)
        [C,ia,ic] = unique(trial);
    end
    trial_id = trial_com_id(trial_com_id>=start_com_id);
    start_trial_no = trial(1);
    end_point = joystick_timestamp(trial_id(end));
    trial_no_timestamp = joystick_timestamp(trial_id);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%   1 Com Marker   %%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%   2 Trial & Phase  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    trial_no_50Hz = (trial_no_timestamp-start_point+1)/(f/f0)+1;
    phase_time_50Hz = (joystick_timestamp-start_point+1)/(f/f0)+1;
    %trial_list =[];
    k=1;
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%补充缺失trial
    %trial_id_lose = find(diff(trial)>1 & diff(trial)<10)+1;
    trial_id_lose = find(diff(trial_com_id)>16);
    if size(trial_id_lose,2)>0
        for i = 1:size(trial_id_lose,2)
            lose_next_pre = phase1(find(trial_com_id(trial_id_lose(i))<find(joystick_comwrite==99) & find(joystick_comwrite==99) < trial_com_id(trial_id_lose(i)+1)));
            lose_next_pre = lose_next_pre(end);
            lose_last_inter = phase_end(find(trial_com_id(trial_id_lose(i))<find(joystick_comwrite==104) & find(joystick_comwrite==104) < trial_com_id(trial_id_lose(i)+1)));
            lose_last_inter = lose_last_inter(1);
            if size(lose_next_pre,1) == 1
                trial_no_50Hz = [trial_no_50Hz,phase_time_50Hz(lose_next_pre)];
            else
                if size(lose_last_inter,1) == 1
                    trial_no_50Hz = [trial_no_50Hz,phase_time_50Hz(lose_last_inter)];
                end
            end
        end
    end
    trial_no_50Hz = sort(trial_no_50Hz);
    trial_list = 1:1:size(trial_no_50Hz,2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%补充缺失trial

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

    trial_no = zeros(1,phase_time_50Hz(end));
    for i = 1:(size(trial_no_50Hz,2)-1)
        trial_no(1,trial_no_50Hz(i):(trial_no_50Hz(i+1)-1))=trial_list(i);
        %timestamp_pre_joystick = joystick_comwrite(phase_pre_joystick(i));
    end
    trial_no(1,trial_no_50Hz(size(trial_no_50Hz,2)):end) = trial_list(end);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%   2 Trial & Phase  %%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%  3  Single Unit   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numOfNeurons = size(readNeurons,2);
    spikes_50Hz = [];
    for n_id = 1: numOfNeurons
        spike0 = zeros(numOfTimestamps, 1);
        spike0(timeStamps(neuron_id == readNeurons(n_id)),1) = 1;
        spike_nid_50Hz=[];
        for t = 1:(size(trial_no,2)-1)
            spike_nid_50Hz = [spike_nid_50Hz,sum(spike0((start_point+(t-1)*(f/f0)+1):(start_point+t*f/f0),:))];
        end
        spike_nid_50Hz = [spike_nid_50Hz,sum(spike0((start_point+(size(trial_no,2)-1)*(f/f0)+1):end,:))];
        spikes_50Hz = [spikes_50Hz;spike_nid_50Hz];
    end
    clear NEV;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%  3  Single Unit   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Combine NEV data
end