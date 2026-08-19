clc;
clear;
close all;

%% 1. INSERIMENTO VALORI PRE-DESIGN E CALCOLI DA CEA
thrust = 55000;         % Spinta desiderata [N]
eps = 200;              % Rapporto d'espansione ugello (Ae/At)
pch_bar = 20;           % Pressione di camera [bar]
pch = pch_bar * 1e5;    % Pressione di camera [Pa]
O_F0_design = 2.5;      % Rapporto O/F nominale di progetto
GOX0_design = 500;      % Flusso massico di ossidante [kg/(m^2 s)]
t_burn = 300;           % Tempo di combustione [s]
g0 = 9.80665;           % Accelerazione di gravità [m/s^2]

% --- Limiti Accettabili per L/D ---
LD_min = 3.0;
LD_max = 5.0;

% --- Proprietà Propellente (HTPB / LOX) ---
rho_f = 920;            % Densità combustibile HTPB [kg/m^3]
a_rf = 0.027 * 1e-3;    % Coefficiente legge di regressione [m/s / (kg/m^2 s)^n]
n_rf = 0.75;            % Esponente legge di regressione [-]

% --- Estrazione Dati Termodinamici da CEA ---
load("CEA_functions.mat");
Tch = T_fun_of_p(O_F0_design, pch_bar);
k = k_fun_of_p(O_F0_design, pch_bar);
R = R_fun_of_p(O_F0_design, pch_bar);

K2 = k * ((2 / (k + 1))^((k + 1) / (k - 1)));
cstar = sqrt(R * Tch / K2);

eps_fun = @(Me, k_val) (1 / Me) * sqrt(((1 + 0.5 * (k_val - 1) * (Me^2)) / (0.5 * (k_val + 1)))^((k_val + 1) / (k_val - 1)));
Me = fzero(@(Me) eps - eps_fun(Me, k), 2);
pe = pch / ((1 + 0.5 * (k - 1) * (Me^2))^(k / (k - 1)));
Te = Tch / (1 + 0.5 * (k - 1) * (Me^2));
ue = Me * sqrt(k * R * Te);

mdot_tot = thrust / (ue + pe * eps * (cstar / pch));
At = cstar * mdot_tot / pch;
dt = 2 * sqrt(At / pi);

fox = mdot_tot / (O_F0_design + 1);
mox = O_F0_design * fox;


%% 2. CALCOLO DELLA LUNGHEZZA DEL GRANO (L)
Ap0 = mox / GOX0_design;
r0 = sqrt(Ap0 / pi);
D0 = 2 * r0;
Pb0 = pi * D0;
rf0 = a_rf * GOX0_design^n_rf;

L_grain = fox / (rho_f * Pb0 * rf0);

fprintf('=== GEOMETRIA INIZIALE DEL GRANO (t = 0 s) ===\n');
fprintf('Diametro porto (D0) : %.4f m\n', D0);
fprintf('Lunghezza grano (L) : %.4f m\n', L_grain);
fprintf('L/D0 iniziale       : %.2f\n\n', L_grain / D0);


%% 3. CALCOLO DIAMETRO ESTERNO (ANALITICO + MAKE MESH)
r_end = (r0^(2*n_rf + 1) + (2*n_rf + 1) * a_rf * (mox / pi)^n_rf * t_burn)^(1 / (2*n_rf + 1));
D_ext = 2 * r_end;
LD_ext_nominal = L_grain / D_ext;

fprintf('=== REGRESSIONE FINALE (t = %d s) ===\n', t_burn);
fprintf('Raggio finale (r_end) : %.4f m\n', r_end);
fprintf('Diametro esterno D_ext: %.4f m\n', D_ext);
fprintf('L/D_ext finale        : %.2f\n\n', LD_ext_nominal);

% --- Discretizzazione Mesh ---
type = "cylinder";
n_mesh = 500;
vars.geometry.diametro = D0;

save('prevars.mat', 'L_grain', 'D0', 'mox', 'eps', 'Ap0', 'Pb0', 'pch_bar', 'type', '-mat');
[coord, coordmesh, coord_v, idx_v, idx_v_interni, idx_v_esterni] = make_mesh0("cylinder", n_mesh, vars, type);

figure('Name', 'Mesh Geometria Iniziale');
plot(coordmesh(:,1), coordmesh(:,2), 'b.');
axis equal; grid on;
title('Discretizzazione Mesh del Porto Circolare (t = 0 s)');
xlabel('x [m]'); ylabel('y [m]');


%% 4. ED 5. MAPPA L/D, SHIFT DI O/F E DOPPIA OTTIMIZZAZIONE
O_F_vec = linspace(1.5, 4.5, 60);
GOX_vec = linspace(500, 800, 60);
t_vec_opt = linspace(0, t_burn, 100);

LD_est = nan(length(O_F_vec), length(GOX_vec));
delta_OF_mat = zeros(length(O_F_vec), length(GOX_vec));

for i = 1:length(O_F_vec)
    for j = 1:length(GOX_vec)
        O_F0 = O_F_vec(i);
        GOX0 = GOX_vec(j);
        
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
        
        Ap0_k = mox_k / GOX0;
        r0_k  = sqrt(Ap0_k / pi);
        Pb0_k = 2 * pi * r0_k;
        rf0_k = a_rf * GOX0^n_rf;
        L_k   = fox0 / (rho_f * Pb0_k * rf0_k);
        
        % L/D_ext finale
        r_end_k = (r0_k^(2*n_rf + 1) + (2*n_rf + 1) * a_rf * (mox_k / pi)^n_rf * t_burn)^(1 / (2*n_rf + 1));
        LD_est(i,j) = L_k / (2 * r_end_k);  %% POSSIBILE SISTEMARE DIMENSIONE
        
        % Shift di O/F
        r_t_k   = (r0_k^(2*n_rf + 1) + (2*n_rf + 1) * a_rf * (mox_k / pi)^n_rf .* t_vec_opt).^(1 / (2*n_rf + 1));
        Pb_t_k  = 2 * pi * r_t_k;
        Ap_t_k  = pi * r_t_k.^2;
        Gox_t_k = mox_k ./ Ap_t_k;
        fox_t_k = rho_f .* Pb_t_k .* L_k .* (a_rf .* Gox_t_k.^n_rf);
        OF_t_k  = mox_k ./ fox_t_k;
        
        delta_OF_mat(i,j) = max(OF_t_k) - min(OF_t_k);
    end
end

% --- 1. OTTIMO ASSOLUTO (Senza Vincoli) ---
[min_val_global, min_idx_global] = min(delta_OF_mat(:));
[best_i_g, best_j_g] = ind2sub(size(delta_OF_mat), min_idx_global);
best_OF0_global  = O_F_vec(best_i_g);
best_GOX0_global = GOX_vec(best_j_g);

% --- 2. OTTIMO VINCOLATO (Dentro L/D in [3, 5]) ---
mask_valid = (LD_est >= LD_min) & (LD_est <= LD_max);
delta_OF_constrained = delta_OF_mat;
delta_OF_constrained(~mask_valid) = Inf; % Esclude punti fuori range

[min_val_constr, min_idx_constr] = min(delta_OF_constrained(:));
[best_i_c, best_j_c] = ind2sub(size(delta_OF_constrained), min_idx_constr);
best_OF0_constr  = O_F_vec(best_i_c);
best_GOX0_constr = GOX_vec(best_j_c);

fprintf('=== CONFRONTO OTTIMIZZAZIONE ===\n');
fprintf('1) OTTiMO ASSOLUTO (Senza vincolo L/D):\n');
fprintf('   O/F_0 = %.2f, G_ox,0 = %.0f kg/(m^2 s) -> Shift = %.3f, L/D_ext = %.2f\n', ...
    best_OF0_global, best_GOX0_global, min_val_global, LD_est(best_i_g, best_j_g));

fprintf('2) OTTiMO VINCOLATO (L/D tra %.1f e %.1f):\n', LD_min, LD_max);
fprintf('   O/F_0 = %.2f, G_ox,0 = %.0f kg/(m^2 s) -> Shift = %.3f, L/D_ext = %.2f\n\n', ...
    best_OF0_constr, best_GOX0_constr, min_val_constr, LD_est(best_i_c, best_j_c));


% --- PLOT MAPPA CON AREA AMMISSIBILE ED EVIDENZIAZIONE DEI PUNTI ---
figure('Name', 'Mappa Shift O/F e Regione L/D Ammissibile', 'Position', [150 150 950 650]);
[X, Y] = meshgrid(GOX_vec, O_F_vec);

% Fondo con la scala di colori del Delta O/F
contourf(X, Y, delta_OF_mat, 30, 'LineColor', 'none');
colormap parula; c = colorbar; c.Label.String = '\Delta O/F (Shift Massimo)';
hold on;

% Oscuramento trasparente dell'area NON ammissibile
invalid_mask = double(~mask_valid);
contourf(X, Y, invalid_mask, [0.5 0.5], 'FaceColor', [0.2 0.2 0.2], 'FaceAlpha', 0.45, 'LineColor', 'none');

% Isocline standard di L/D
[C_ld, h_ld] = contour(X, Y, LD_est, [2 2.5 3 4 5 6 8 10], 'k:', 'LineWidth', 1.0);
clabel(C_ld, h_ld, 'Color', 'k', 'FontSize', 8);

% Bordi rossi spessi dell'AREA CONSENTITA (L/D in [3, 5])
[C_b, h_b] = contour(X, Y, LD_est, [LD_min LD_max], 'r', 'LineWidth', 2.5);
clabel(C_b, h_b, 'Color', 'r', 'FontSize', 11, 'FontWeight', 'bold');

% Plot Punti Ottimi
p_global = plot(best_GOX0_global, best_OF0_global, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
p_constr = plot(best_GOX0_constr, best_OF0_constr, 'gd', 'MarkerSize', 13, 'MarkerFaceColor', 'g', 'LineWidth', 1.5);

title('Mappa \Delta O/F con Regione L/D_{ext} \in [3, 5] ed Evidenziazione Ottimi');
xlabel('G_{ox,0} [kg/(m^2 s)]'); ylabel('O/F_0');
legend([p_global, p_constr], ...
    {sprintf('Ottimo Assoluto (O/F=%.2f, G=%.0f)', best_OF0_global, best_GOX0_global), ...
     sprintf('Ottimo Vincolato L/D (O/F=%.2f, G=%.0f)', best_OF0_constr, best_GOX0_constr)}, ...
    'Location', 'northeast');
grid on;


%% 6. SIMULAZIONE TEMPORALE E PLOT DEI PARAMETRI (2x3)

% COMMENTARE PER GIOCARE CON I VALORI


% design rispetto a parametri calcolati ottimali

design_points = [ ...
    O_F0_design, GOX0_design;        % Nominal di partenza
    best_OF0_global, best_GOX0_global; % Ottimo Assoluto
    best_OF0_constr, best_GOX0_constr];% Ottimo Vincolato L/D

%design con valori inseriti
% 
% OF_vect = [2,2,2]';
% GOX0_vect = [500,400,300]';
% design_points = [OF_vect,GOX0_vect ];


n_dp = size(design_points, 1);
t_sim = linspace(0, t_burn, 300);
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
    
    Ap0_sim = mox_sim / GOX0_sim;
    r0_sim  = sqrt(Ap0_sim / pi);
    Pb0_sim = 2 * pi * r0_sim;
    rf0_sim = a_rf * GOX0_sim^n_rf;
    L_sim   = fox0_sim / (rho_f * Pb0_sim * rf0_sim);
    
    r_t_sim   = (r0_sim^(2*n_rf + 1) + (2*n_rf + 1) * a_rf * (mox_sim / pi)^n_rf .* t_sim).^(1 / (2*n_rf + 1));
    Ap_t_sim  = pi * r_t_sim.^2;
    Pb_t_sim  = 2 * pi * r_t_sim;
    Gox_t_sim = mox_sim ./ Ap_t_sim;
    fox_t_sim = rho_f * Pb_t_sim * L_sim .* (a_rf * Gox_t_sim.^n_rf);
    mdot_t_sim = mox_sim + fox_t_sim;
    OF_t_sim  = mox_sim ./ fox_t_sim;
    
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
    results(idp).Itot = trapz(t_sim, thrust_t);
end

figure('Name', 'Analisi Temporale Parametri', 'Position', [100 100 1200 700]);
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