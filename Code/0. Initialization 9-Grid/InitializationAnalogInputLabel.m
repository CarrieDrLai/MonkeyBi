function [sensor,trajectory,stable,trajectory2] = AnalogInputLabel(raw_data)

    for input_id = size(raw_data.ElectrodesInfo,2)-15 :size(raw_data.ElectrodesInfo,2)
        
        if raw_data.ElectrodesInfo(input_id).Label(1:6) == 'sensor' 
            if raw_data.ElectrodesInfo(input_id).Label(1:7) == 'sensoro'
                continue
            else
                if raw_data.ElectrodesInfo(input_id).Label(1:7) == 'sensorS'
                    continue
                else
                    sensor = raw_data.Data(input_id,:);
                    [input_id]
                    break
                end
            end
        end
        %{
        if NS6.ElectrodesInfo(input_id).Label(1:2)=='x '
            x = NS6.Data(input_id,:);
        end
        if NS6.ElectrodesInfo(input_id).Label(1:2)=='y '
            y = NS6.Data(input_id,:);
        end
        if NS6.ElectrodesInfo(input_id).Label(1:12)=='sensorStable'
            stable = NS6.Data(input_id,:);
        end
        if NS6.ElectrodesInfo(input_id).Label(1:2)=='x2'
            x2 = NS6.Data(input_id,:);
        end
        if NS6.ElectrodesInfo(input_id).Label(1:2)=='y2'
            y2 = NS6.Data(input_id,:);
        end
        %}
        
    end
    x = raw_data.Data(input_id+1,:);
    y = raw_data.Data(input_id+2,:);
    stable = raw_data.Data(input_id+3,:);
    x2 = raw_data.Data(input_id+4,:);
    y2 = raw_data.Data(input_id+5,:);
    
    sensor = sensor./32768;
    trajectory = [x;y]./32768;
    stable = stable./32768;
    trajectory2 = [x2;y2]./32768;
end