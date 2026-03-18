% STOCHASTIC EXTENSION: Hybrid ODE/Gillespie SSA
% Implementing Hybrid-Stochastic Approach to Bridge et al. model.
%
% Hybrid Method:
%     1. Run the full ODE to get the aGTP(t) timecourse
%        (upstream signalling is deterministic / high-copy)
%     2. Use aGTP(t) as a time-varying production rate to
%        drive a 2-reaction Gillespie simulation of cAMP:
%          R1: aGTP -> aGTP + cAMP   rate = kCAMP  * aGTP(t)
%          R2: cAMP -> 0             rate = kCAMP_ * cAMP
%     3. Repeat 500 times per ligand, sample endpoint cAMP

%   Volume was chosen so cAMP endpoint ~1000 molecules.
%   S(t) at endpoint ~ 6.5, a = 5.518e11
%   [cAMP]_end = 6.5/5.518e11 = 1.18e-11 M
%   V = 1000 / (1.18e-11 * 6.022e23) = 1.41e-10 L
%=========================================================

clear; clc; rng(42); % set seed

% Constants
NA     = 6.022e23;
a_sens = 5.518e11; % biosensor scaling factor
V_cell = 1.41e-10; % L  (gives ~1000 cAMP molecules)
scale  = NA * V_cell; % Molarity to molecule count

% Simulation settings
N_sim  = 500; % sims per ligand
t_end  = 60; % minutes as per paper
L_conc = 100e-9; % 100 nM (matches model init conds)

fprintf('=== Hybrid ODE/Gillespie Extension ===\n');
fprintf('Volume = %.2e L\n', V_cell);
fprintf('Expected cAMP molecules at endpoint: ~%.0f\n\n', ...
    (6.5/a_sens) * scale); % confirmed around 1000

% PARAMETER SETS  (Table 4, neuron/Cell A column)
% ExF1
pF.kL_    = 0.2485   *60;  pF.kL     = 2.233e7  *60;
pF.kAct_  = 347.67   *60;  pF.kAct   = 0.0015   *60;
pF.kG_    = 5.103e-6 *60;  pF.kG     = 4.726e6  *60;
pF.delta_ = 1.885e-4;      pF.delta  = 0.0390;
pF.ski_   = 1.8551;        pF.ski    = 1.2346e4;
pF.mu_    = 0.0345;        pF.mu     = 0.5206;
pF.kGTP   = 0.0144   *60;  pF.kDS    = 1.299e-5 *60;
pF.v_     = 0.9644;        pF.v      = 0.9390;
pF.thetaA = 0.8668;
pF.kHyd_  = 9.126e-4 *60;  pF.kHyd   = 0.1393   *60;
pF.kCAMP_ = 0.0964   *60;  pF.kCAMP  = 0.1236   *60;
pF.kGRA_  = 3.144e-19*60;  pF.kGRA   = 8.213e3  *60;
pF.RTot   = 6.591e-10;     pF.GTot   = 1.405e-9;
pF.L      = L_conc;

% ExD3
pD.kL_    = 0.0134   *60;  pD.kL     = 1.922e7  *60;
pD.kAct_  = 1.790e4  *60;  pD.kAct   = 0.2438   *60;
pD.kG_    = 0.0289   *60;  pD.kG     = 7.780e4  *60;
pD.delta_ = 1.944e-4;      pD.delta  = 0.0311;
pD.ski_   = 0.9776;        pD.ski    = 2.2611e4;
pD.mu_    = 0.0384;        pD.mu     = 1.9944;
pD.kGTP   = 0.0993   *60;  pD.kDS    = 0.0018   *60;
pD.v_     = 0.1683;        pD.v      = 0.7571;
pD.thetaA = 4.2125;
pD.kHyd_  = 2.492e-15*60;  pD.kHyd   = 0.0634   *60;
pD.kCAMP_ = 0.1268   *60;  pD.kCAMP  = 4.8893   *60;
pD.kGRA_  = 6.516e-19*60;  pD.kGRA   = 2.029e6  *60;
pD.RTot   = 6.591e-10;     pD.GTot   = 1.405e-9;
pD.L      = L_conc;

% ODE pre-equilibration (L=0, 600 min) similar paper
fprintf('ODE pre-equilibration\n');
ode_opts = odeset('MaxStep',1,'RelTol',1e-8,'AbsTol',1e-20);

pF0=pF; pF0.L=0;
y0=zeros(15,1); y0(15)=pF.RTot; y0(5)=pF.GTot;
[~,ytmp] = ode15s(@(t,y)odefull(t,y,pF0),[0 600],y0,ode_opts);
ic_F = ytmp(end,:)';

pD0=pD; pD0.L=0;
y0=zeros(15,1); y0(15)=pD.RTot; y0(5)=pD.GTot;
[~,ytmp] = ode15s(@(t,y)odefull(t,y,pD0),[0 600],y0,ode_opts);
ic_D = ytmp(end,:)';

% Full ODE at 100 nM to get aGTP(t) timecourse
fprintf('Full ODE at 100 nM to get aGTP(t)\n');
t_query = linspace(0, t_end, 3600);   % get value every second

[~, yF] = ode15s(@(t,y)odefull(t,y,pF), t_query, ic_F, ode_opts);
[~, yD] = ode15s(@(t,y)odefull(t,y,pD), t_query, ic_D, ode_opts);

aGTP_F = yF(:,2);   % aGTP timecourse, ExF1 (Molar)
aGTP_D = yD(:,2);   % aGTP timecourse, ExD3 (Molar)

% Convert to molecule counts ONCE beforehand
aGTP_F_mols = max(aGTP_F * scale, 0);
aGTP_D_mols = max(aGTP_D * scale, 0);

% ODE based cAMP signal
S_ode_F = a_sens * yF(:,1);
S_ode_D = a_sens * yD(:,1);

fprintf('  ODE ExF1 endpoint S(60) = %.4f\n', S_ode_F(end));
fprintf('  ODE ExD3 endpoint S(60) = %.4f\n\n', S_ode_D(end));

% Hybrid Gillespie on cAMP only
fprintf('Hybrid Gillespie (%d trajectories per ligand)\n', N_sim);

cAMP_F = zeros(N_sim,1);
cAMP_D = zeros(N_sim,1);

% Initial cAMP molecule count from ODE Init Cond
cAMP0_F = max(round(ic_F(1) * scale), 0);
cAMP0_D = max(round(ic_D(1) * scale), 0);
fprintf('  Initial cAMP molecules: ExF1=%d, ExD3=%d\n\n', cAMP0_F, cAMP0_D);

fprintf('ExF1:\n');
for i = 1:N_sim
    if mod(i,100)==0, fprintf('  %d/%d\n',i,N_sim); end
    cAMP_F(i) = hybrid_gillespie(t_query, aGTP_F_mols, cAMP0_F, ...
                                  pF.kCAMP, pF.kCAMP_, t_end);
end

fprintf('ExD3:\n');
for i = 1:N_sim
    if mod(i,100)==0, fprintf('  %d/%d\n',i,N_sim); end
    cAMP_D(i) = hybrid_gillespie(t_query, aGTP_D_mols, cAMP0_D, ...
                                  pD.kCAMP, pD.kCAMP_, t_end);
end

% Convert molecule # back to signal units
S_F = a_sens * (cAMP_F / scale);
S_D = a_sens * (cAMP_D / scale);

% Stats
fprintf('\n========== RESULTS ==========\n');
fprintf('ExF1:  mean=%.4f  std=%.4f  CV=%.1f%%\n', ...
    mean(S_F), std(S_F), 100*std(S_F)/mean(S_F));
fprintf('ExD3:  mean=%.4f  std=%.4f  CV=%.1f%%\n', ...
    mean(S_D), std(S_D), 100*std(S_D)/mean(S_D));
fprintf('ODE ExF1 endpoint: %.4f\n', S_ode_F(end));
fprintf('ODE ExD3 endpoint: %.4f\n', S_ode_D(end));
fprintf('Mean difference (ExF1-ExD3): %.4f\n', mean(S_F)-mean(S_D));

[h_ks, p_ks, ks_stat] = kstest2(S_F, S_D);
fprintf('\nTwo-sample KS test:\n');
fprintf('  KS statistic = %.4f,  p = %.4f\n', ks_stat, p_ks);
if h_ks == 0
    fprintf('  Distributions NOT significantly different (p>0.05)\n');
    fprintf('  --> Convergence holds at single-cell level.\n');
else
    fprintf('  Distributions ARE significantly different (p<=0.05)\n');
    fprintf('  --> ODE convergence is a mean-field artifact.\n');
end

% Endpoint distribution figure
figure('Position',[100 100 800 520]); hold on;
h1 = histogram(S_F, 40, 'FaceColor',[0.18 0.55 0.90], ...
    'FaceAlpha',0.65,'EdgeColor','none');
h2 = histogram(S_D, 40, 'FaceColor',[0.95 0.35 0.25], ...
    'FaceAlpha',0.65,'EdgeColor','none');
xline(mean(S_F),'Color',[0.05 0.25 0.75],'LineWidth',2.5);
xline(mean(S_D),'Color',[0.75 0.10 0.05],'LineWidth',2.5);
yl = ylim;
text(mean(S_F), yl(2)*0.90, sprintf('ExF1\n\\mu=%.2f', mean(S_F)), ...
    'Color',[0.05 0.25 0.75],'FontSize',11,'FontWeight','bold', ...
    'HorizontalAlignment','center');
text(mean(S_D), yl(2)*0.78, sprintf('ExD3\n\\mu=%.2f', mean(S_D)), ...
    'Color',[0.75 0.10 0.05],'FontSize',11,'FontWeight','bold', ...
    'HorizontalAlignment','center');
xlabel('S(t = 60 min) = a \cdot [cAMP_i]','FontSize',13);
ylabel('Number of simulated cells','FontSize',13);
title({sprintf('Single-Cell cAMP Distribution at t = 60 min'), ...
    sprintf('Neuron, [L]=100 nM, N=%d per ligand | KS p=%.3f', N_sim, p_ks)}, ...
    'FontSize',12);
legend([h1 h2],{'ExF1 (G protein-biased)','ExD3 (balanced)'}, ...
    'FontSize',11,'Location','best');
box on; set(gca,'FontSize',11);
exportgraphics(gcf,'fig1_endpoint_distributions.png','Resolution',300);

% FIGURE 2: ODE mean-field trajectories + stochastic endpoint
figure('Position',[100 100 800 480]); hold on;
plot(t_query, S_ode_F, '-', 'Color',[0.18 0.55 0.90], ...
    'LineWidth',3,'DisplayName','ExF1 (ODE mean-field)');
plot(t_query, S_ode_D, '-', 'Color',[0.95 0.35 0.25], ...
    'LineWidth',3,'DisplayName','ExD3 (ODE mean-field)');
errorbar(62, mean(S_F), std(S_F), 'o', ...
    'Color',[0.05 0.25 0.75],'MarkerFaceColor',[0.05 0.25 0.75], ...
    'MarkerSize',9,'LineWidth',2.5, ...
    'DisplayName',sprintf('ExF1 stoch \\mu=%.2f \\sigma=%.2f', mean(S_F),std(S_F)));
errorbar(62, mean(S_D), std(S_D), 's', ...
    'Color',[0.75 0.10 0.05],'MarkerFaceColor',[0.75 0.10 0.05], ...
    'MarkerSize',9,'LineWidth',2.5, ...
    'DisplayName',sprintf('ExD3 stoch \\mu=%.2f \\sigma=%.2f', mean(S_D),std(S_D)));
xline(60,'--k','LineWidth',1,'HandleVisibility','off');
text(61,1,'stoch.','FontSize',9,'Color',[0.4 0.4 0.4]);
xlabel('t (mins)','FontSize',13);
ylabel('S(t) = a \cdot [cAMP_i]','FontSize',13);
title({'Deterministic ODE vs Stochastic Ensemble','Neuron, [L]=100 nM'},'FontSize',13);
legend('FontSize',10,'Location','northeast');
xlim([0 67]); ylim([0 40]);
xticks([0 10 20 30 40 50 60]);
box on; set(gca,'FontSize',11);
exportgraphics(gcf,'fig2_ode_vs_stochastic.png','Resolution',300);

fprintf('\nDone. Figures saved.\n');

%=========================================================
% Hybrid Gillespie for cAMP only
% R1: production  likelihood = kCAMP  * aGTP_mols(t)
% R2: degradation likelihood = kCAMP_ * cAMP_mols
%=========================================================
function cAMP_end = hybrid_gillespie(t_vec, aGTP_mols_vec, cAMP0, ...
                                      kCAMP, kCAMP_, t_end)

    cAMP  = cAMP0;
    t     = 0;
    dt    = t_vec(2) - t_vec(1); % same step size
    N_t   = length(t_vec);

    while t < t_end
        idx = min(floor(t / dt) + 1, N_t);
        aGTP_mols = aGTP_mols_vec(idx);

        a1 = kCAMP  * aGTP_mols;
        a2 = kCAMP_ * cAMP;
        a0 = a1 + a2;

        if a0 == 0, break; end

        % Time to next event
        tau = -log(rand) / a0;
        t   = t + tau;
        if t > t_end, break; end

        % Randomness to simulate a random chance of next reaction occurring
        % based one # of cAMP and aGTP
        if rand * a0 < a1
            cAMP = cAMP + 1;
        else
            cAMP = max(cAMP - 1, 0);
        end
    end
    cAMP_end = cAMP;
end

% Full ODE based on paper
function dydt = odefull(t, y, p)
    dydt = zeros(15,1);
    L = p.L;
    dydt(1)  = -y(1)*p.kCAMP_ + y(2)*p.kCAMP;
    dydt(2)  = -y(2)*p.kHyd + y(4)*p.kHyd_ + y(6)*p.thetaA*p.v_*p.kGTP + y(8)*p.kGTP;
    dydt(3)  = -y(3)*y(4)*p.kGRA + y(5)*p.kGRA_ + y(6)*p.thetaA*p.v_*p.kGTP + y(8)*p.kGTP;
    dydt(4)  =  y(2)*p.kHyd + y(4)*(-p.kHyd_ - y(3)*p.kGRA) + y(5)*p.kGRA_;
    dydt(5)  =  y(3)*y(4)*p.kGRA + y(5)*(-p.kG*y(15) - p.v*p.kG*y(12) ...
                - p.mu*p.kG*y(13) - p.mu*p.v*p.kG*y(10) - p.kGRA_) ...
                + y(6)*p.mu_*p.v_*p.kG_ + y(7)*p.v_*p.kG_ ...
                + y(8)*p.mu_*p.kG_ + y(9)*p.kG_;
    dydt(6)  =  y(5)*y(10)*p.mu*p.v*p.kG + y(6)*(-p.mu_*p.v_*p.kG_ ...
                - p.ski_*p.v_*p.kL_ - p.mu_*p.ski_*p.kAct_ - p.thetaA*p.v_*p.kGTP) ...
                + y(7)*p.mu*p.ski*p.kAct + y(8)*p.ski*p.v*p.kL*L;
    dydt(7)  =  y(5)*y(12)*p.v*p.kG + y(6)*p.mu_*p.ski_*p.kAct_ ...
                + y(7)*(-p.v_*p.kG_ - p.v_*p.kL_ - p.mu*p.ski*p.kAct) + y(9)*p.v*p.kL*L;
    dydt(8)  =  y(5)*y(13)*p.mu*p.kG + y(6)*p.ski_*p.v_*p.kL_ ...
                + y(8)*(-p.mu_*p.kG_ - p.ski*p.v*p.kL*L - p.mu_*p.kAct_ - p.kGTP) ...
                + y(9)*p.mu*p.kAct;
    dydt(9)  =  y(5)*y(15)*p.kG + y(7)*p.v_*p.kL_ + y(8)*p.mu_*p.kAct_ ...
                + y(9)*(-p.kG_ - p.v*p.kL*L - p.mu*p.kAct);
    dydt(10) =  y(6)*(p.mu_*p.v_*p.kG_ + p.thetaA*p.v_*p.kGTP) ...
                + y(10)*(-p.ski_*p.kAct_ - p.ski_*p.kL_ - p.mu*p.v*p.kG*y(5) - p.kDS) ...
                + y(12)*p.ski*p.kAct + y(13)*p.ski*p.kL*L;
    dydt(11) =  y(10)*p.kDS - y(11)*p.delta_*p.kL_ + y(14)*p.delta*p.kL*L;
    dydt(12) =  y(7)*p.v_*p.kG_ + y(10)*p.ski_*p.kAct_ ...
                + y(12)*(-p.kL_ - p.ski*p.kAct - p.v*p.kG*y(5)) + p.kL*L*y(15);
    dydt(13) =  y(8)*(p.mu_*p.kG_ + p.kGTP) + y(10)*p.ski_*p.kL_ ...
                + y(13)*(-p.kAct_ - p.ski*p.kL*L - p.mu*p.kG*y(5)) + y(15)*p.kAct;
    dydt(14) =  p.delta_*p.kL_*y(11) - p.delta*p.kL*L*y(14);
    dydt(15) =  y(9)*p.kG_ + y(12)*p.kL_ + y(13)*p.kAct_ ...
                + y(15)*(-p.kL*L - p.kAct - p.kG*y(5));
end