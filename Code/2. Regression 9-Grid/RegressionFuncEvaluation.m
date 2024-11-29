function [R_sqr,MSE,CC] = RegressionFuncEvaluation(y_num,test_len,y_predict,test_y)
    R_sqr = [];
    MSE = [];
    CC = [];
    for var = 1: y_num
        R_sqr = [R_sqr,(test_len*sum(y_predict(:,var).*test_y(:,var))-sum(test_y(:,var))*sum(y_predict(:,var)))^2/((test_len*sum((test_y(:,var)).^2)-sum(test_y(:,var))^2)*(test_len*sum((y_predict(:,var)).^2)-sum(y_predict(:,var))^2))];
        MSE = [MSE,1/test_len * sum((test_y(:,var)-y_predict(:,var)).^2)];
        CC = [CC,corr(y_predict(:,var),(test_y(:,var)))];
    end
end