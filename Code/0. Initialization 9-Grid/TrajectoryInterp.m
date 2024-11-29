function new_array = TrajectoryInterp(target_len,x)

        x2=zeros([1,size(x,2)*2-1]);
        x2(2:2:end) = 0.5*(x(2:end)+x(1:end-1));
        x2(1:2:end) = x;
        x3=zeros([1,size(x2,2)*2-1]);
        x3(2:2:end) = 0.5*(x2(2:end)+x2(1:end-1));
        x3(1:2:end) = x2;
        x4=zeros([1,size(x3,2)*2-1]);
        x4(2:2:end) = 0.5*(x3(2:end)+x3(1:end-1));
        x4(1:2:end) = x3;
        x5=zeros([1,size(x4,2)*2-1]);
        x5(2:2:end) = 0.5*(x4(2:end)+x4(1:end-1));
        x5(1:2:end) = x4;
        x6=zeros([1,size(x5,2)*2-1]);
        x6(2:2:end) = 0.5*(x5(2:end)+x5(1:end-1));
        x6(1:2:end) = x5;        
        
        x2=zeros([1,size(x6,2)*2-1]);
        x2(2:2:end) = 0.5*(x6(2:end)+x6(1:end-1));
        x2(1:2:end) = x6;
        x3=zeros([1,size(x2,2)*2-1]);
        x3(2:2:end) = 0.5*(x2(2:end)+x2(1:end-1));
        x3(1:2:end) = x2;
        x4=zeros([1,size(x3,2)*2-1]);
        x4(2:2:end) = 0.5*(x3(2:end)+x3(1:end-1));
        x4(1:2:end) = x3;
        x5=zeros([1,size(x4,2)*2-1]);
        x5(2:2:end) = 0.5*(x4(2:end)+x4(1:end-1));
        x5(1:2:end) = x4;
        x6=zeros([1,size(x5,2)*2-1]);
        x6(2:2:end) = 0.5*(x5(2:end)+x5(1:end-1));
        x6(1:2:end) = x5;                
        
        if target_len >size(x6,2)
            x7=zeros([1,size(x6,2)*2-1]);
            x7(2:2:end) = 0.5*(x6(2:end)+x6(1:end-1));
            x7(1:2:end) = x6;      
            x8=zeros([1,size(x7,2)*2-1]);
            x8(2:2:end) = 0.5*(x7(2:end)+x7(1:end-1));
            x8(1:2:end) = x7;     
            x9=zeros([1,size(x8,2)*2-1]);
            x9(2:2:end) = 0.5*(x8(2:end)+x8(1:end-1));
            x9(1:2:end) = x8;
            new_array = x9(1:floor(size(x9,2)/target_len):end);
            new_array = new_array(1:target_len);
        else
            new_array = x6(1:floor(size(x6,2)/target_len):end);
            new_array = new_array(1:target_len);
        end
end