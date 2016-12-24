function [ac, cp, c_ac] = control_aero(N,M,c_r,c,sweep,y)
% %     %TEST
%     N = 3;
%     M = N-1;
%     c_r = 10;
%     c = [5 10 5];
%     sweep = 0;
%     y = [-10 0 10];
    
    %array de c i de y passat a punts mitjos
    c1 = c(1:M);
    c2 = c(2:N);
        
    y1 = y(1:M); %N punts M panells
    y2 = y(2:N); 
    
    y_ac = y1 + abs(y2 - y1)/2;
    y_cp = y_ac;
    
    c_ac = c1 + (c2 - c1)/2; 
        
    x_ac = c_r/4 + tan(sweep)*y_ac;
    x_cp = x_ac + c_ac/2; 
    
    ac = [x_ac; y_ac; zeros(1,M)];
    cp = [x_cp; y_cp; zeros(1,M)];

end