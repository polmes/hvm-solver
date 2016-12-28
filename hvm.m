function [CL,Cl,CM,CD] = hvm(AR,lambda,sweep,alpha_cr,alpha_ct,epsilon_cr,epsilon_ct,alpha,beta,M,dist,flap,x_h,y_inf,y_sup,eta)
    if (nargin == 0)
        clc;
        
        disp('WING PARAMETERS');
        AR = input('Aspect Ratio: ');
        lambda = input('Taper Ratio: ');
        sweep = input('Quarter-chord Sweep Angle (º): ');
        
        disp('AERODYNAMIC TWIST');
        alpha_cr = input('Root airfoil alpha_l0 (º): ');
        alpha_ct = input('Tip airfoil alpha_l0 (º): ');
        
        disp('GEOMETRIC TWIST');
        epsilon_cr = input('Root airfoil epsilon (º): ');
        epsilon_ct = input('Tip airfoil epsilon (º): ');
        
        disp('AERODYNAMIC CONDITIONS');
        alpha = input('Angle of attack (º): ');
        beta = input('Cross-wind angle (º): ');
        
        disp('SIMULATION PARAMETERS');
        M = input('Number of panels: ');
        dist = input('(a) Uniform\n(b) Full cosine\nDistribution for geometry discretization (choose from the above a or b): ','s');
        
        disp('FLAP ANALYSIS');
        flap = input('Include a flap? (y/n) ','s');
    elseif (nargin < 12)
        error('Not enough input arguments');
    elseif (nargin > 16)
        error('Too many input arguments');
    end
    
    % Aspect Ratio
    if (AR < 0)
        error('Invalid Aspect Ratio');
    elseif (AR < 5)
        warning('Horshoe Vortex Method may not provide accurate results for such a low aspect ratio');
    end
    
    % Taper Ratio
    if (lambda < 0 || lambda > 1)
        error('Invalid Taper Ratio');
    end
    
    % Root and tip chords
    c_r = 2 / (AR * (lambda + 1));
    c_t = lambda * c_r;
    
    % Angles
    sweep = deg2rad(sweep);
    alpha_cr = deg2rad(alpha_cr);
    alpha_ct = deg2rad(alpha_ct);
    epsilon_cr = deg2rad(epsilon_cr);
    epsilon_ct = deg2rad(epsilon_ct);
    alpha = deg2rad(alpha);
    beta = deg2rad(beta);
    
    % Global variables
    N = M + 1;
    S = 1/AR; % area, assuming b = 1
    Vinf = [cos(alpha)*cos(beta) cos(alpha)*sin(beta) sin(alpha)]; % unit vector

    % Geometry discretization distribution
    if (dist == 'a')
        y = linspace(-0.5,0.5,N); % uniform
    elseif (dist == 'b')
        k = 0:M;
        y = 1/2 * (1 - cos(k/M * pi)) - 0.5; % full cosine
    else
        error('Invalid type of distribution');
    end
    
    % Chords
    c = zeros(1,N);
    right = (y >= 0);
    c(right) = c_r + (c_t - c_r) * 2 .* y(right);
    left = (y < 0);
    c(left) = c_r - (c_t - c_r) * 2 .* y(left);
    
    % Panel lengths
    y1 = y(1:M);
    y2 = y(2:N);
    delta_y = y2 - y1;
    
    % Mid Chords
    c1 = c(1:M);
    c2 = c(2:N);
    c_AC = c1 + (c2 - c1)/2;
    
    % Aerodynamic Centers
    y_AC = y1 + delta_y/2;
    x_AC = c_r/4 + tan(sweep) * abs(y_AC);
    AC = [x_AC; y_AC; zeros(1,M)];
    
    % Control Points
    y_CP = y_AC;
    x_CP = x_AC + c_AC/2; 
    CP = [x_CP; y_CP; zeros(1,M)];
    CP_right = (y_CP >= 0);
    CP_left = (y_CP < 0);
    
    % Aerodynamic twist
    alpha_l0 = zeros(1,M);
    alpha_l0(CP_right) = alpha_cr + (alpha_ct - alpha_cr) * 2 .* y_CP(CP_right);
    alpha_l0(CP_left) = alpha_cr - (alpha_ct - alpha_cr) * 2 .* y_CP(CP_left);

    % Geometric twist
    epsilon = zeros(1,M);
    epsilon(CP_right) = epsilon_cr + (epsilon_ct - epsilon_cr) * 2 .* y_CP(CP_right);
    epsilon(CP_left) = epsilon_cr - (epsilon_ct - epsilon_cr) * 2. * y_CP(CP_left);

    % Flap parameters
    if (flap == 'y' || flap == 'Y')
        flap = true;
        
        if (nargin == 0)
            x_h = input('Flap hinge x position (in tenths of chord): ');
            y_inf = input('Flap hinge y start position (in tenths of semispan): ');
            y_sup = input('Flap hinge y end position (in tenths of semispan): ');
            eta = input('Flap deflection angle (º): ');
        elseif (nargin < 16)
            error('Not enough input arguments for flap analysis');
        end
    else
        flap = false;
    end
    
    % Flap analysis
    if (flap)
        % Fix input parameters
        x_h = x_h / 10;
        y_inf = y_inf / 20;
        y_sup = y_sup / 20;
        eta = deg2rad(eta);
        
        % Thin Airfoil Theory
        theta_h = acos(1 - 2 * x_h .* c_AC);
        delta_alpha_l0 = zeros(1,M);
        reg_flap = (y_CP >= y_inf & y_CP <= y_sup | y_CP <= -y_inf & y_CP >= -y_sup);
        delta_alpha_l0(reg_flap) = -(1 - theta_h(reg_flap)/pi + sin(theta_h(reg_flap))/pi) * eta;
        
        % Correction factor
        x_factor = deg2rad([10,20,30,40,50,60,70]);
        y_factor = [0.8,0.7,0.53,0.45,0.4,0.36,0.34];
        coef_factor = polyfit(x_factor,y_factor,6);
        factor = @(x) coef_factor(1)*x.^6 + coef_factor(2)*x.^5 + coef_factor(3)*x.^4+ coef_factor(4)*x.^3+ coef_factor(5)*x.^2 + coef_factor(6)*x + + coef_factor(7);
        delta_alpha_l0 = delta_alpha_l0 * factor(eta);
        alpha_l0 = alpha_l0 + delta_alpha_l0;
    end

    % Normal vectors
    n = [sin(epsilon - alpha_l0); zeros(1,M); cos(epsilon - alpha_l0)];
    
    % Bounded Vortices
    BV1 = [AC(1,:); AC(2,:) - delta_y/2; AC(3,:)];
    BV2 = [AC(1,:); AC(2,:) + delta_y/2; AC(3,:)];
    
    % Trailing Vortices
    tmp = 20 * ones(1,M);
    TV1 = [tmp; y1 + tmp  * tan(beta); tmp  * tan(alpha)];
    TV2 = [tmp; y2 + tmp  * tan(beta); tmp  * tan(alpha)]; 

    % i's and j's
    ij = ndgrid(1:M,1:M); % combvec equivalent
    in = zeros(2,M*M); % indices
    in(1,:) = reshape(ij,[1 M*M]);
    in(2,:) = reshape(ij',[1 M*M]);
    
    % A to B, B to C, C to D
    aa = TV1(:,in(1,:));
    bb = BV1(:,in(1,:));
    cc = BV2(:,in(1,:));
    dd = TV2(:,in(1,:));
    
    % Calculate induced velocity at P by vortex line from 1 to 2
    function v = uvw(pp,p1,p2)
        r0 = p2 - p1;
        r1 = pp - p1;
        r2 = pp - p2;
        r1xr2 = cross(r1,r2);
        v = 1/(4*pi) * bsxfun(@times, dot(r0, bsxfun(@rdivide, r1, sqrt(sum(r1.^2))) - bsxfun(@rdivide, r2, sqrt(sum(r2.^2))), 1), bsxfun(@rdivide, r1xr2, dot(r1xr2, r1xr2, 1)));
    end
    
    % Induced velocity at CP (i) by horshoe vortex (j)
    v = uvw(CP(:,in(2,:)),aa,bb) + uvw(CP(:,in(2,:)),bb,cc) + uvw(CP(:,in(2,:)),cc,dd);
    wake = uvw(AC(:,in(2,:)),aa,bb) + uvw(AC(:,in(2,:)),cc,dd);
    
    % System of equations
    A = reshape(dot(v,n(:,in(2,:)),1),[M M])'; % influence coefficients
    RHS = -(Vinf * n)'; % normal velocity
    gamma = linsolve(A,RHS)'; % circulation
    
    % Lift coefficient
    CL = (2/S) * (delta_y * gamma'); % whole wing value
    Cl = 2 ./ (delta_y .* c_AC) .* (gamma .* delta_y); % distribution
    
    % Pitching moment coefficient about the leading edge
    c_aero = (2/3 * c_r) * (1 + lambda + lambda^2) / (1 + lambda);
    CM = -2/(S * c_aero) * sum(gamma .* AC(1,:) .* delta_y) * cos(alpha);
    
    % Induced drag coefficient
    w_wake = reshape(wake(3,:),[M M])';
    CD = -(2/S) * sum(gamma .* delta_y .* (gamma * w_wake));
    
    % Plots
    figure;
    plot(y_CP,Cl);
    title('Lift distribution along the span');
    xlabel('$$\frac{y}{b}$$','Interpreter','latex');
    ylabel('$$C_l$$','Interpreter','latex');
end
