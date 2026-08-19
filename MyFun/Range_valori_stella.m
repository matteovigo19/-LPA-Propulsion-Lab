clc;
clear;
close all;

%% 1. INPUT DI PRE-DESIGN E TERMODINAMICA CEA
thrust = 55000;         % Spinta desiderata [N]
eps = 200;              % Rapporto d'espansione ugello (Ae/At)
pch_bar = 20;           % Pressione di camera [bar]
pch = pch_bar * 1e5;    % Pressione di camera [Pa]
O_F0_design = 2.5;      % O/F nominale
GOX0_design = 500;      % Flusso massico ossidante [kg/(m^2 s)]
t_burn = 300;           % Tempo di combustione [s]
g0 = 9.80665;           % Accelerazione di gravità [m/s^2]

% --- Parametri Geometria Stella a 6 Punte ---
N_pts = 6;              % Numero di punte della stella
K_p   = 1.8;            % Fattore di forma perimetro (P_stella = K_p * P_cerchio)
LD_min = 1.5;           % Limite inferiore L/D_ext
LD_max = 4.0;           % Limite superiore L/D_ext

% --- Proprietà Propellente (HTPB / LOX) ---
rho_f = 920;            % Densità combustibile [kg/m^3]
a_rf = 0.027 * 1e-3;    % Coeff. regressione [m/s / (kg/m^2 s)^n]
n_rf = 0.75;            % Esponente regressione [-]

% --- Estrazione Dati Termodinamici da CEA ---
load("CEA_functions.mat");
Tch = T_fun_of_p(O_F0_design, pch_bar);
k   = k_fun_of_p(O_F0_design, pch_bar);
R   = R_fun_of_p(O_F0_design, pch_bar);

K2 = k * ((2 / (k + 1))^((k + 1) / (k - 1)));
cstar = sqrt(R * Tch / K2);

eps_fun = @(Me, k_val) (1 / Me) * sqrt(((1 + 0.5 * (k_val - 1) * (Me^2)) / (0.5 * (k_val + 1)))^((k_val + 1) / (k_val - 1)));
Me = fzero(@(Me) eps - eps_fun(Me, k), 2);
pe = pch / ((1 + 0.5 * (k - 1) * (Me^2))^(k / (k - 1)));
Te = Tch / (1 + 0.5 * (k - 1) * (Me^2));
ue = Me * sqrt(k * R * Te);

mdot_tot = thrust / (ue + pe * eps * (cstar / pch));
At = cstar * mdot_tot / pch;

fox = mdot_tot / (O_F0_design + 1);
mox = O_F0_design * fox;


%% 2. GEOMETRIA INIZIALE STELLA A 6 PUNTE (t = 0 s)
Ap0 = mox / GOX0_design;               % Area porto iniziale [m^2]
r_eq0 = sqrt(Ap0 / pi);                % Raggio equivalente circolare [m]
Pb0 = (2 * pi * r_eq0) * K_p;          % Perimetro maggiorato dalla stella a 6 punte [m]

rf0 = a_rf * GOX0_design^n_rf;         % Rateo regressione iniziale [m/s]
L_grain = fox / (rho_f * Pb0 * rf0);   % Lunghezza grano [m]

fprintf('=== GEOMETRIA STELLA A 6 PUNTE (t = 0 s) ===\n');
fprintf('Numero punte        : %d\n', N_pts);
fprintf('Area porto (Ap0)    : %.5f m^2\n', Ap0);
fprintf('Perimetro iniziale  : %.4f m (Fattore K_p = %.1f)\n', Pb0, K_p);
fprintf('Lunghezza grano (L) : %.4f m\n\n', L_grain);


%% 3. MAPPA L/D, SHIFT O/F E OTTIMIZZAZIONE
O_F_vec = linspace(1.5, 6, 60);
GOX_vec = linspace(200, 800, 60);
t_vec_opt = linspace(0, t_burn, 100);
dt_step = t_burn / (length(t_vec_opt) - 1);

LD_est = nan(length(O_F_vec), length(GOX_vec));
delta_OF_mat = zeros(length(O_F_vec), length(GOX_vec));

for i = 1:length(O_F_vec)
    for j = 1:length(GOX_vec)
        O_F0 = O_F_vec(i);
        GOX0 = GOX_vec(j);
        
        % TermodinamicaCEA locale
        Tch0 = T_fun_of_p(O_F0, pch_bar);
        k0   = k_fun_of_p(O_F0, pch_bar);
        R0   = R_fun_of_p(O_F0, pch_bar);
        K2_0 = k0 * ((2 / (k0 + 1))^((k0 + 1) / (k0 - 1)));
        cstar0 = sqrt(R0 * Tch0 / K2_0);
        
        Me0 = fzero(@(Me) eps - eps_fun(Me, k0), 2);
        pe0 = pch / ((1 + 0.5 * (k0 - 1) * Me0^2)^(k0 / (k0 - 1)));
        Te0 = Tch0 / (1 + 0.5 * (k0 - 1) * Me0^2);
        ue0 = Me0 * sqrt(k0 * R0 * Te0);
        
        mdot_0 = thrust / (ue0 + pe0 * eps * (cstar0 / pch));
        fox0   = mdot_0 / (O_F0 + 1);
        mox_k  = O_F0 * fox0;
        
        % Inizializzazione Geometria Stella
        Ap0_k  = mox_k / GOX0;
        r_eq_k = sqrt(Ap0_k / pi);
        Pb0_k  = (2 * pi * r_eq_k) * K_p;
        rf0_k  = a_rf * GOX0^n_rf;
        L_k    = fox0 / (rho_f * Pb0_k * rf0_k);
        
        % Integrazione Temporale Numerica per Stella a 6 Punte
        y_web = 0;
        OF_t_k = zeros(size(t_vec_opt));
        
        for nt = 1:length(t_vec_opt)
            r_eff = r_eq_k + y_web;
            Ap_t  = Ap0_k + Pb0_k * y_web + pi * (y_web^2);
            % Transizione geometrica: il perimetro recede fino al valore circolare equivalente
            Pb_t  = max(2 * pi * r_eff, Pb0_k - 2 * N_pts * y_web);
            
            Gox_t = mox_k / Ap_t;
            rf_t  = a_rf * (Gox_t)^n_rf;
            fox_t = rho_f * Pb_t * L_k * rf_t;
            
            OF_t_k(nt) = mox_k / fox_t;
            y_web = y_web + rf_t * dt_step; % Avanzamento web
        end
        
        % Diametro esterno finale e L/D_ext
        D_ext_k = 2 * (r_eq_k + y_web);
        LD_est(i,j) = L_k / D_ext_k;
        delta_OF_mat(i,j) = max(OF_t_k) - min(OF_t_k);
    end
end

% --- Ottimo Assoluto ---
[min_val_global, min_idx_global] = min(delta_OF_mat(:));
[best_i_g, best_j_g] = ind2sub(size(delta_OF_mat), min_idx_global);
best_OF0_global  = O_F_vec(best_i_g);
best_GOX0_global = GOX_vec(best_j_g);

% --- Ottimo Vincolato nell'area L/D in [3, 5] ---
mask_valid = (LD_est >= LD_min) & (LD_est <= LD_max);
delta_OF_constrained = delta_OF_mat;
delta_OF_constrained(~mask_valid) = Inf;

[min_val_constr, min_idx_constr] = min(delta_OF_constrained(:));
[best_i_c, best_j_c] = ind2sub(size(delta_OF_constrained), min_idx_constr);
best_OF0_constr  = O_F_vec(best_i_c);
best_GOX0_constr = GOX_vec(best_j_c);

fprintf('=== RISULTATI OTTIMIZZAZIONE (GRANO A STELLA 6 PUNTE) ===\n');
fprintf('1) OTTiMO ASSOLUTO:\n   O/F_0 = %.2f, G_ox,0 = %.0f kg/(m^2 s) -> Shift = %.3f, L/D_ext = %.2f\n', ...
    best_OF0_global, best_GOX0_global, min_val_global, LD_est(best_i_g, best_j_g));
fprintf('2) OTTiMO VINCOLATO (L/D_ext in [%.1f, %.1f]):\n   O/F_0 = %.2f, G_ox,0 = %.0f kg/(m^2 s) -> Shift = %.3f, L/D_ext = %.2f\n\n', ...
    LD_min, LD_max, best_OF0_constr, best_GOX0_constr, min_val_constr, LD_est(best_i_c, best_j_c));


%% 4. GRAFICO MAPPA VALORI ACCETTABILI ED EVIDENZIAZIONE REGIONE
figure('Name', 'Mappa Shift O/F - Stella a 6 Punte', 'Position', [150 150 950 650]);
[X, Y] = meshgrid(GOX_vec, O_F_vec);

% Superficie del Delta O/F
contourf(X, Y, delta_OF_mat, 30, 'LineColor', 'none');
colormap parula; c = colorbar; c.Label.String = '\Delta O/F (Shift Massimo)';
hold on;

% Oscuramento trasparente delle regioni NON ammissibili (L/D fuori range)
invalid_mask = double(~mask_valid);
contourf(X, Y, invalid_mask, [0.5 0.5], 'FaceColor', [0.2 0.2 0.2], 'FaceAlpha', 0.45, 'LineColor', 'none');

% Isocline L/D generali
[C_ld, h_ld] = contour(X, Y, LD_est, [2 2.5 3 3.5 4 4.5 5 6 8], 'k:', 'LineWidth', 1.0);
clabel(C_ld, h_ld, 'Color', 'k', 'FontSize', 8);

% Bordi Rossi Spessi per L/D_ext ammissibile [3, 5]
[C_b, h_b] = contour(X, Y, LD_est, [LD_min LD_max], 'r', 'LineWidth', 2.5);
clabel(C_b, h_b, 'Color', 'r', 'FontSize', 11, 'FontWeight', 'bold');

% Marcatori dei Punti di Ottimo
p_global = plot(best_GOX0_global, best_OF0_global, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
p_constr = plot(best_GOX0_constr, best_OF0_constr, 'gd', 'MarkerSize', 13, 'MarkerFaceColor', 'g', 'LineWidth', 1.5);

title('Grano a Stella (6 Punte): Mappa \Delta O/F con Regione Ammissibile L/D_{ext} \in [1.5, 4]');
xlabel('G_{ox,0} [kg/(m^2 s)]'); ylabel('O/F_0');
legend([p_global, p_constr], ...
    {sprintf('Ottimo Assoluto (O/F=%.2f, G=%.0f)', best_OF0_global, best_GOX0_global), ...
     sprintf('Ottimo Vincolato L/D (O/F=%.2f, G=%.0f)', best_OF0_constr, best_GOX0_constr)}, ...
    'Location', 'northeast');
grid on;


%% 5. SIMULAZIONE TEMPORALE E PLOT PARAMETRI (2x3)
design_points = [ ...
    O_F0_design, GOX0_design;          % Nominal iniziale
    best_OF0_global, best_GOX0_global;   % Ottimo Assoluto
    best_OF0_constr, best_GOX0_constr];  % Ottimo Vincolato

n_dp = size(design_points, 1);
t_sim = linspace(0, t_burn, 300);
dt_sim = t_burn / (length(t_sim) - 1);
results = struct();

for idp = 1:n_dp
    O_F0_sim = design_points(idp, 1);
    GOX0_sim = design_points(idp, 2);
    
    Tch0 = T_fun_of_p(O_F0_sim, pch_bar);
    k0   = k_fun_of_p(O_F0_sim, pch_bar);
    R0   = R_fun_of_p(O_F0_sim, pch_bar);
    K2_0 = k0 * ((2 / (k0 + 1))^((k0 + 1) / (k0 - 1)));
    cstar0 = sqrt(R0 * Tch0 / K2_0);
    
    Me0 = fzero(@(Me) eps - eps_fun(Me, k0), 2);
    pe0 = pch / ((1 + 0.5 * (k0 - 1) * Me0^2)^(k0 / (k0 - 1)));
    Te0 = Tch0 / (1 + 0.5 * (k0 - 1) * Me0^2);
    ue0 = Me0 * sqrt(k0 * R0 * Te0);
    
    mdot_0 = thrust / (ue0 + pe0 * eps * (cstar0 / pch));
    At_sim = cstar0 * mdot_0 / pch;
    
    fox0_sim = mdot_0 / (O_F0_sim + 1);
    mox_sim  = O_F0_sim * fox0_sim;
    
    Ap0_sim  = mox_sim / GOX0_sim;
    r_eq_sim = sqrt(Ap0_sim / pi);
    Pb0_sim  = (2 * pi * r_eq_sim) * K_p;
    rf0_sim  = a_rf * GOX0_sim^n_rf;
    L_sim    = fox0_sim / (rho_f * Pb0_sim * rf0_sim);
    
    % Evoluzione temporale della stella
    y_web_sim = 0;
    OF_t_sim  = zeros(size(t_sim));
    Gox_t_sim = zeros(size(t_sim));
    mdot_t_sim = zeros(size(t_sim));
    
    for idx = 1:length(t_sim)
        r_eff = r_eq_sim + y_web_sim;
        Ap_t  = Ap0_sim + Pb0_sim * y_web_sim + pi * (y_web_sim^2);
        Pb_t  = max(2 * pi * r_eff, Pb0_sim - 2 * N_pts * y_web_sim);
        
        Gox_t_sim(idx) = mox_sim / Ap_t;
        rf_t  = a_rf * (Gox_t_sim(idx))^n_rf;
        fox_t = rho_f * Pb_t * L_sim * rf_t;
        
        mdot_t_sim(idx) = mox_sim + fox_t;
        OF_t_sim(idx)   = mox_sim / fox_t;
        
        y_web_sim = y_web_sim + rf_t * dt_sim;
    end
    
    % Calcolo punto fisso dinamica di camera
    pch_t = zeros(size(t_sim)); Tch_t = pch_t; thrust_t = pch_t; Isp_t = pch_t;
    pch_guess = pch_bar;
    
    for idx = 1:length(t_sim)
        OF_i = OF_t_sim(idx);
        for iter = 1:5
            Tch_i = T_fun_of_p(OF_i, pch_guess);
            k_i   = k_fun_of_p(OF_i, pch_guess);
            R_i   = R_fun_of_p(OF_i, pch_guess);
            K2_i  = k_i * ((2 / (k_i + 1))^((k_i + 1) / (k_i - 1)));
            cstar_i = sqrt(R_i * Tch_i / K2_i);
            pch_new_bar = mdot_t_sim(idx) * cstar_i / At_sim / 1e5;
            if abs(pch_new_bar - pch_guess) < 1e-4
                pch_guess = pch_new_bar;
                break;
            end
            pch_guess = pch_new_bar;
        end
        
        pch_t(idx) = pch_guess * 1e5;
        Tch_t(idx) = Tch_i;
        
        Me_i = fzero(@(Me) eps - eps_fun(Me, k_i), 2);
        pe_i = pch_t(idx) / ((1 + 0.5 * (k_i - 1) * Me_i^2)^(k_i / (k_i - 1)));
        Te_i = Tch_i / (1 + 0.5 * (k_i - 1) * Me_i^2);
        ue_i = Me_i * sqrt(k_i * R_i * Te_i);
        
        thrust_t(idx) = mdot_t_sim(idx) * ue_i + pe_i * eps * At_sim;
        Isp_t(idx)    = thrust_t(idx) / (mdot_t_sim(idx) * g0);
    end
    
    results(idp).O_F0 = O_F0_sim; results(idp).GOX0 = GOX0_sim; results(idp).t = t_sim;
    results(idp).pch_bar = pch_t / 1e5; results(idp).Tch = Tch_t; results(idp).OF = OF_t_sim;
    results(idp).Gox = Gox_t_sim; results(idp).thrust = thrust_t; results(idp).Isp = Isp_t;
end

% Subplot 2x3
figure('Name', 'Analisi Temporale Stella 6 Punte', 'Position', [100 100 1200 700]);
labels = {'Nominale (O/F=2.5, G=500)', ...
          sprintf('Ottimo Assoluto (O/F=%.2f, G=%.0f)', best_OF0_global, best_GOX0_global), ...
          sprintf('Ottimo Vincolato L/D (O/F=%.2f, G=%.0f)', best_OF0_constr, best_GOX0_constr)};

subplot(2,3,1); hold on; for i=1:n_dp, plot(results(i).t, results(i).pch_bar, 'LineWidth', 1.5); end
title('Pressione di Camera'); xlabel('t [s]'); ylabel('p_c [bar]'); grid on; legend(labels, 'Location', 'best');

subplot(2,3,2); hold on; for i=1:n_dp, plot(results(i).t, results(i).Tch, 'LineWidth', 1.5); end
title('Temperatura di Camera'); xlabel('t [s]'); ylabel('T_c [K]'); grid on;

subplot(2,3,3); hold on; for i=1:n_dp, plot(results(i).t, results(i).OF, 'LineWidth', 1.5); end
title('Rapporto Miscela (O/F)'); xlabel('t [s]'); ylabel('O/F [-]'); grid on;

subplot(2,3,4); hold on; for i=1:n_dp, plot(results(i).t, results(i).Gox, 'LineWidth', 1.5); end
title('Flusso Massico Ossidante'); xlabel('t [s]'); ylabel('G_{ox} [kg/(m^2 s)]'); grid on;

subplot(2,3,5); hold on; for i=1:n_dp, plot(results(i).t, results(i).thrust, 'LineWidth', 1.5); end
title('Spinta'); xlabel('t [s]'); ylabel('T [N]'); grid on;

subplot(2,3,6); hold on; for i=1:n_dp, plot(results(i).t, results(i).Isp, 'LineWidth', 1.5); end
title('Impulso Specifico'); xlabel('t [s]'); ylabel('I_{sp} [s]'); grid on;