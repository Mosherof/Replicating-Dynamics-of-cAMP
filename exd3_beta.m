%=========================================================
% FULL SYSTEM ODE FUNCTION
%=========================================================
function dydt = fullSystem(t, y, p)
    dydt = zeros(15, 1); 
    L = p.L;

    % 1. dCAMPi
    dydt(1) = -y(1)*p.kCAMP_ + y(2)*p.kCAMP;

    % 2. alpha GTP
    dydt(2) = -y(2)*p.kHyd + y(4)*p.kHyd_ + y(6)*p.thetaA*p.v_*p.kGTP + y(8)*p.kGTP;

    % 3. beta gamma
    dydt(3) = -y(3)*y(4)*p.kGRA + y(5)*p.kGRA_ + y(6)*p.thetaA*p.v_*p.kGTP + y(8)*p.kGTP; 

    % 4. alpha GDP
    dydt(4) = y(2)*p.kHyd + y(4)*(-p.kHyd_ - y(3)*p.kGRA) + y(5)*p.kGRA_; 

    % 5. G
    dydt(5) = y(3)*y(4)*p.kGRA + y(5)*(-p.kG*y(15) - p.v*p.kG*y(12) - p.mu*p.kG*y(13) ...
        - p.mu*p.v*p.kG*y(10) - p.kGRA_) + y(6)*p.mu_*p.v_*p.kG_ + y(7)*p.v_*p.kG_ ...
        + y(8)*p.mu_*p.kG_ + y(9)*p.kG_; 

    % 6. LR*G
    dydt(6) = y(5)*y(10)*p.mu*p.v*p.kG + y(6)*(-p.mu_*p.v_*p.kG_ - p.ski_*p.v_*p.kL_ ...
        - p.mu_*p.ski_*p.kAct_ - p.thetaA*p.v_*p.kGTP) + y(7)*p.mu*p.ski*p.kAct ...
        + y(8)*p.ski*p.v*p.kL*L;

    % 7. LRG 
    dydt(7) = y(5)*y(12)*p.v*p.kG + y(6)*p.mu_*p.ski_*p.kAct_ + y(7)*(-p.v_*p.kG_ - ...
        p.v_*p.kL_ - p.mu*p.ski*p.kAct)+ y(9)*p.v*p.kL*L; 
    
    % 8. R*G
    dydt(8) = y(5)*y(13)*p.mu*p.kG + y(6)*p.ski_*p.v_*p.kL_ + y(8)*(-p.mu_*p.kG_ - ...
        p.ski*p.v*p.kL*L - p.mu_*p.kAct_ - p.kGTP) + y(9)*p.mu*p.kAct; 

    % 9. RG
    dydt(9) = y(5)*y(15)*p.kG + y(7)*p.v_*p.kL_ + y(8)*p.mu_*p.kAct_ + y(9)*(-p.kG_ - ...
        p.v*p.kL*L - p.mu*p.kAct); 

    % 10. LR*
    dydt(10) =  y(6)*(p.mu_*p.v_*p.kG_ + p.thetaA*p.v_*p.kGTP)+ ... 
        y(10)*(-p.ski_*p.kAct_ - p.ski_*p.kL_- p.mu*p.v*p.kG*y(5)...
        - p.kDS) + y(12)*p.ski*p.kAct + y(13)*(p.ski*p.kL*L);

    % 11. LR DS
    dydt(11) = y(10)*p.kDS - y(11)*p.delta_*p.kL_ + y(14)*p.delta*p.kL*L; 

    % 12. LR
    dydt(12) =  y(7)*p.v_*p.kG_ + y(10)*p.ski_*p.kAct_ + y(12)*(-p.kL_ - ...
        p.ski*p.kAct - p.v*p.kG*y(5)) + p.kL*L*y(15); 

    % 13. R*
    dydt(13) = y(8)*(p.mu_*p.kG_ + p.kGTP) + y(10)*p.ski_*p.kL_ + ...
        y(13)*(-p.kAct_ - p.ski*p.kL*L - p.mu*p.kG*y(5)) + y(15)*p.kAct; 

    % 14. R DS
    dydt(14) = p.delta_*p.kL_*y(11) - p.delta*p.kL*L*y(14); 

    % 15. R
    dydt(15) =  y(9)*p.kG_ + y(12)*p.kL_ + y(13)*p.kAct_ + y(15)*(-p.kL*L - p.kAct - p.kG*y(5)) ;
end

%=========================================================
% INITIAL CONDITIONS FUNCTION
%=========================================================
function y0 = initialConditions(p)
    y0 = zeros(15,1);
    y0(15) = p.RTot;  
    y0(5)  = p.GTot;  
end

%=========================================================
% MAIN SCRIPT
%=========================================================
clear; clc;

% PARAMETERS
p.kL_ = 0.0134* 60;
p.kL = 1.922 * 10^7* 60;
p.kAct_ = 1.790 *10^4* 60;
p.kAct = 0.2438* 60;
p.kG_ = 0.0289* 60;
p.kG = 7.780 * 10 ^ 4* 60; 
p.delta_ = 1.486 *10^(-16);
p.delta = 2.552 *10^(-9);
p.ski_ = 0.9976;
p.ski = 2.2611 *10 ^4; 
p.mu_ = 0.0384;
p.mu = 1.9944;
p.kGTP = 0.0993* 60;
p.kDS = 6.341e-04* 60;
p.v_ = 0.1683;
p.v = 0.7571;
p.thetaA = 4.2125;
p.kHyd_ =2.492*10^(-15)* 60;
p.kHyd = 0.0634* 60;
p.kCAMP_ = 0.1268* 60;
p.kCAMP = 4.8893* 60;
p.RTot = 6.648 * (10^-9);
p.GTot = 2.231*10^(-10);
p.kRA_ = 6.516*10^(-19)* 60;
p.kRA = 2.029*10^(6)* 60;
p.kGRA_ = p.kRA_;
p.kGRA = p.kRA;

a = 5.518 * 10 ^ 11;
tspan = [0 60]; % 60 minutes
L_2_values = [1e-11, 1e-10, 1e-9, 3e-9, 1e-8, 1e-7];

figure; hold on;

% FIRST SOLVER: L = 0
p.L = 0;
y0 = initialConditions(p);
options = odeset('MaxStep',1);
[t, y] = ode15s(@(t,y) fullSystem(t,y,p), tspan, y0, options);

% SECOND SOLVER: varying L using steady state from first run
for L = L_2_values
    p.L = L;
    y1 = y(end, :)';  % use steady state from first solver
    [t2, y2] = ode15s(@(t,y2) fullSystem(t,y2,p), tspan, y1, options);
    plot(t2, a .*y2(:,1), 'LineWidth', 3);  % plot cAMP vs time
end

xlabel('t (mins)');
xticks([10 20 30 40 50 60])
xticklabels({'10', '20', '30', '40', '50', '60'})
ylabel('S(t) = a[cAMP]');
ylim([0 40]);
legend('L=1e-11','L=1e-10','L=1e-9','L=3e-9','L=1e-8','L=1e-7');
title('B-cell - ExD3');

f = gcf;
exportgraphics(f, 'exd3 beta.png', 'Resolution', 300);