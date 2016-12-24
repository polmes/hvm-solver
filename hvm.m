% function [CL,Cl,CM,CDi] = hvm()
% user parameters

    % b = 1;
    clear variables;
    lambda = 0.5;
    AR = 10;
    S = 1/AR;
    sweep = 0; % deg2rad(10);
    M = 100;
    alpha_cr = deg2rad(-1.8);
    alpha_ct = deg2rad(-1.8);
    epsilon_cr = 0; % deg2rad(0);
    epsilon_ct = 0; % deg2rad(5);
    % flap_inf = 0.4;
    % flap_sup = 0.7;
    % x_h = 0.8;
    % eta = deg2rad(15);

    alpha = deg2rad(4);
    beta = deg2rad(0);
    Vinf = [cos(alpha)*cos(beta) cos(alpha)*sin(beta) sin(alpha)]; % unit vector
    % Globals
    N = M+1;

    % wing parameters

    cr = 2/(AR*(lambda+1));
    ct = lambda*cr;

    % control_aero vindrà aquí quan tot sigui una f

    % discretization

    y = linspace(-0.5,0.5,N);
    % full cos too (ap. 5)

    c = zeros(1,N);
    c(y>=0) = cr+(ct-cr)*2.*y(y>=0);
    c(y<0) = cr-(ct-cr)*2.*y(y<0);

    [ac,cp,c_ac] = control_aero(N,M,cr,c,sweep,y);
    cp_y = cp(2,:);

    % aerodynamic twist

    alpha_l0(cp_y>=0) = alpha_cr+(alpha_ct-alpha_cr)*2.*cp_y(cp_y>=0);
    alpha_l0(cp_y<0) = alpha_cr-(alpha_ct-alpha_cr)*2.*cp_y(cp_y<0);

    % geometric twist

    epsilon(cp_y>=0) = epsilon_cr+(epsilon_ct-epsilon_cr)*2.*cp_y(cp_y>=0);
    epsilon(cp_y<0) = epsilon_cr-(epsilon_ct-epsilon_cr)*2.*cp_y(cp_y<0);

    % flaps

    % theta_h = acos(1-2*x_h*c_ac); % c_ac la te feta el freeman
    % delta_alpha_l0 = zeros(1,2); % el 2 ha de ser M
    % reg_flap = (y > flap_inf & y < flap_sup | y < -flap_inf & y > -flap_sup);
    % delta_alpha_l0(reg_flap) = -(1-(theta_h(reg_flap)/pi)+(sin(theta_h(reg_flap))/pi))*(eta);
    % alpha_l0 = alpha_l0+delta_alpha_l0;


    % normal vectors

    nx = sin(epsilon - alpha_l0);
    ny = zeros(1,M);
    nz = cos(epsilon - alpha_l0);

    n = [nx;ny;nz];

    % delta_y array

    delta_y = y(2:N) - y(1:M);

    [BV,TV] = bv_tv(ac, alpha, beta, delta_y, y, M);
    
    BV1 = [ac(1,:); ac(2,:) - abs(delta_y)/2; ac(3,:)];
    BV2 = [ac(1,:); ac(2,:) + abs(delta_y)/2; ac(3,:)];
    
    y1 = y(1:M);
    y2 = y(2:N);
    tmp = 20 * ones(1,M);
    TV1 = [tmp; y1 + tmp  * tan(beta); tmp  * tan(alpha)];
    TV2 = [tmp; y2 + tmp  * tan(beta); tmp  * tan(alpha)];
    
    in = combvec(1:M,1:M);
    at = TV1(:,in(1,:));
    bt = BV1(:,in(1,:));
    ct = BV2(:,in(1,:));
    dt = TV2(:,in(1,:));
    
    vt = uvwt(cp(:,in(2,:)),at,bt) + uvwt(cp(:,in(2,:)),bt,ct) + uvwt(cp(:,in(2,:)),ct,dt);
    At = reshape(dot(vt,n(:,in(2,:)),1),[M M])';
    
    A = zeros(M);
    v_array = zeros(3,M,M);
    w_wake = zeros(M);
    for i = 1:M % control points
        for j = 1:M % wakes
            % A to B, B to C, C to D
            aa = TV(:,1,j);
            bb = BV(:,1,j);
            cc = BV(:,2,j);
            dd = TV(:,2,j);

            % Should be able to simplify by using array of v's and a,b,c,d
            % calculated from in(:,:) matrix. Keep uvw() function
            v = uvw(cp(:,i),aa,bb) + uvw(cp(:,i),bb,cc) + uvw(cp(:,i),cc,dd);
            v_array(:,i,j) = v;
            v_wake = uvw(ac(:,i),aa,bb) + uvw(ac(:,i),cc,dd);
            w_wake(i,j) = v_wake(3);
%             v_array(:,i) = v_array(:,i) + v;
            
            A(i,j) = v' * n(:,i);
        end
    end
    RHS = -(Vinf * n)';

    gamma = linsolve(A,RHS)';
    gammat = linsolve(At,RHS)';
    
    Cl = 2 ./ (delta_y .* c_ac) .* (gamma .* delta_y);
    CL = (2/S) * (delta_y * gamma'); % 1/AR = S
    
    CLt = (2/S) * (delta_y * gammat');
    
    % induced drag
%     alpha_i = -v_array(3,:) / Vinf(1);
%     CDi = -(2/S) * sum(gamma .* delta_y * sum(gamma .* v_wake(3,:)));
    
    CDi = 0;
    for i = 1:M
        CDj = 0;
        for j = 1:M
            CDj = CDj + gamma(j) * w_wake(i,j);
        end
        CDi = CDi + CDj * gamma(i) * delta_y(i);
    end
    CDi = -(2/S) * CDi;

    % coefficient moment

    c_aero = (2*cr/3) * ((1+lambda+lambda^2)/(1+lambda));
    CM = -2/(S * c_aero) * sum(gamma .* ac(1,:) .* delta_y) * cos(alpha); % LE
% end
