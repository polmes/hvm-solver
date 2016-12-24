function [BV, TV] = bv_tv(ac, alpha, beta, delta_y, y, M)
    %Vinf = [u_inf v_inf w_inf]; %potser l'usuari hauria de donar alpha
    %BV i TV son matrius tipus (3,2,M)
    
    BV = zeros(3,2,M); 
    
    for i = 1:1:M
        
       BV1x = ac(1,i); 
       BV1y = ac(2,i) - abs(delta_y(i))/2;
       BV1z = 0;
       BV1 = [BV1x; BV1y; BV1z];
       
       BV2x = BV1x;
       BV2y = BV1y + abs(delta_y(i));
       BV2z = 0;
       BV2 = [BV2x; BV2y; BV2z];
       
       
       BV(:,:,i) = [BV1 BV2]; 
       
    end
    
    
    
    TV = zeros(3,2,M); 
    
    for i = 1:M
        
       TV1x = 20;
       TV1y = y(i)+tan(beta)*20;
       TV1z = tan(alpha)*20;
       TV1 = [TV1x; TV1y; TV1z];
       
       TV2x = 20;
       TV2y = y(i+1)+tan(beta)*20;
       TV2z = tan(alpha)*20;
       TV2 = [TV2x; TV2y; TV2z];
       
       %es podria complicar si tenim en compte que els TV de les puntes han
       %d'estar més lluny? 
       
       TV(:,:,i) = [TV1 TV2]; 
       
    end

end
