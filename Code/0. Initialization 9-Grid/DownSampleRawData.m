function [trajectory_downsample] = DownSampleRawData(trajectory_30k,trial_no,start_point,f,f0,f_joystick)

    traject_start_point = uint32(f_joystick*(double(start_point)/double(f)));
    temp = trajectory_30k(:,traject_start_point:traject_start_point+uint32(double(size(trial_no,2)-1)/f0*f_joystick-1));
    temp = reshape(temp,[2,uint32(f_joystick)/uint32(f0),size(temp,2)/(uint32(f_joystick)/uint32(f0))]);
    trajectory_downsample = mean(temp,2);
    trajectory_downsample = reshape(trajectory_downsample,[2,size(trial_no,2)-1]);
end