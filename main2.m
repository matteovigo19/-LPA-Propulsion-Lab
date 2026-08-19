clear;
close all;
clc;

%% STEP 0: data

tic;

% ================= GEOMETRY DATA =================

load("prevars.mat");

type = "star";

multiplyer = 80;

vars.geometry.ntips = n_tips; % da predesign

% Mesh in SI: metri
vars.geometry.diametro_est = diam_e;    % [m]   da predesign
vars.geometry.diametro_int = diam_i;     % [m]   da predesign

vars.geometry.cyl_tol_rel = 0.03;
vars.geometry.cyl_tol_abs = 0.15e-3;   % 0.15 mm


n = multiplyer * vars.geometry.ntips;

% ================= FUEL DATA =================
% FUEL: HTPB
% regression rate: rf = a*GOx^n

rho_f = 920;                  % [kg/m^3]
a_rf  = 0.027;                % [(mm/s)/((kg/m^2 s)^n)]
n_rf  = 0.75;                 % [-]

% Conversione in SI
a_rf = a_rf * 1e-3;           % [(m/s)/((kg/m^2 s)^n)]

% ================= COMBUSTION DATA =================
% OXIDIZER: O2

%mox = 0.015;                  % [kg/s]  da predesign

load("CEA_functions.mat");

OFmin = 1.1;                  % [-]
OFmax = 19.5;                 % [-]

pmin = 101325;                % [Pa]
pmax = 99e5;                  % [Pa]

tmax = 300;                   % [s]
pamb = 0;                     % [Pa]

% ================= ENGINE DATA =================

ext_diameter = 280.1*2e-3;          % [m]
chamber_length = L; %da predesign       % [m]

throat_diameter = 0.15;       % [m]
%eps = 2; da predesign

%At = 0.25*pi*(throat_diameter^2);   % [m^2]

% Camera coerente in SI
D_camera = ext_diameter;      % [m]

% =========== end of data =================


%% ============================================================
%  STEP 1: CREATE INITIAL MESH
% ============================================================

[coordstar_cart, coord_mesh, coord_v, idx_v] = ...
    make_mesh0(type, n, vars, "cartesiane");

r_v = vecnorm(coord_mesh(idx_v,:), 2, 2);

r_thr = 0.5 * (min(r_v) + max(r_v));

idx_v_interni = idx_v(r_v < r_thr);
idx_v_esterni = idx_v(r_v >= r_thr);


%% ============================================================
%  CHECK INDICI
% ============================================================

n_mesh = size(coord_mesh,1);

fprintf("\nDEBUG MESH INIZIALE:\n");
fprintf("n_mesh = %d\n", n_mesh);
fprintf("max(idx_v) = %d\n", max(idx_v));
fprintf("max(idx_v_interni) = %d\n", max(idx_v_interni));
fprintf("max(idx_v_esterni) = %d\n", max(idx_v_esterni));

if max(idx_v) > n_mesh || max(idx_v_interni) > n_mesh || max(idx_v_esterni) > n_mesh
    error("Indici punte non coerenti con coord_mesh.");
end


%% ============================================================
%  OUTPUT MESH INIZIALE
% ============================================================

figure;
plot(coord_mesh(:,1)*1e3, coord_mesh(:,2)*1e3, 'b.-', ...
    'DisplayName', 'mesh iniziale');
hold on;

plot(coord_mesh(idx_v,1)*1e3, coord_mesh(idx_v,2)*1e3, 'ko', ...
     'MarkerSize', 8, ...
     'LineWidth', 1.5, ...
     'DisplayName', 'tutti i vertici');

plot(coord_mesh(idx_v_interni,1)*1e3, coord_mesh(idx_v_interni,2)*1e3, 'ro', ...
     'MarkerSize', 10, ...
     'LineWidth', 2, ...
     'DisplayName', 'punte interne');

plot(coord_mesh(idx_v_esterni,1)*1e3, coord_mesh(idx_v_esterni,2)*1e3, 'gs', ...
     'MarkerSize', 8, ...
     'LineWidth', 1.5, ...
     'DisplayName', 'punte esterne');

axis equal;
grid on;
legend('Location','best');
xlabel('x [mm]');
ylabel('y [mm]');
title('Controllo indici vertici');


%% ============================================================
%  STEP 2: INITIALIZE CHAMBER, STEADY STATE
% ============================================================

[perim0, Ap0] = eval_mesh(coord_mesh, "cartesian");

tol=1e-10;
if Ap-Ap0>tol
    fprintf("incongreunza area predesign e evalmesh")
    return
end

tol=1e-10;
if Pb-perim0>tol
    fprintf("incongreunza perimetro predesign e evalmesh")
    return
end

%CHECK SU PERIM0 E AP0 controlla coerenza con predesign

fprintf("\nGeometria iniziale:\n");
fprintf("perimeter = %.6e m\n", perim0);
fprintf("area      = %.6e m^2\n", Ap0);

Gox_test = mox / Ap0;
rf_test = a_rf * Gox_test^n_rf;

fprintf("Gox iniziale stimato = %.6f kg/(m^2 s)\n", Gox_test);
fprintf("rf iniziale stimata  = %.6e m/s = %.6f mm/s\n", ...
    rf_test, rf_test*1e3);
fprintf("regressione stimata in tmax = %.6e m = %.6f mm\n", ...
    rf_test*tmax, rf_test*tmax*1e3);


%% === DIMENSIONAMENTO AREA DI GOLA IN BASE ALLA PRESSIONE ===
pc_target_bar = pch_bar;               % [bar] INSERISCI QUI LA TUA PRESSIONE INIZIALE VOLUTA
pc_target = pc_target_bar * 1e5;  % [Pa]

% 1. Calcolo del flusso iniziale (indipendente dalla pressione)
Gox0 = mox / Ap0;
rf0 = a_rf * Gox0^n_rf;
Ab0 = perim0 * chamber_length;
mdot_f0 = rho_f * Ab0 * rf0;

% 2. Portata totale e rapporto di miscela
mdot_tot = mox + mdot_f0;
OF0 = mox / mdot_f0;

% 3. Calcolo del c* (usando le funzioni CEA interpolate alla pressione target)
Tc0 = T_fun_of_p(OF0, pc_target_bar);
k0 = k_fun_of_p(OF0, pc_target_bar);
R0 = R_fun_of_p(OF0, pc_target_bar);

% Formula del coefficiente di Vandenkerckhove (al quadrato)
K2 = k0 * ((2/(k0+1))^((k0+1)/(k0-1))); 
cstar0 = sqrt(R0 * Tc0 / K2);

% 4. Area e diametro di gola
At = (mdot_tot * cstar0) / pc_target;
throat_diameter = 2 * sqrt(At / pi);


%% Compact data in vars

% geometry
vars.geometry.port_area = Ap0;
vars.geometry.burning_perimeter = perim0;
vars.geometry.grain_length = chamber_length;
vars.geometry.throat_area = At;

vars.geometry.idx_v = idx_v;
vars.geometry.idx_v_interni = idx_v_interni;
vars.geometry.idx_v_esterni = idx_v_esterni;
vars.geometry.D_camera = D_camera;
vars.geometry.type = type;
vars.geometry.wall_tol = 1e-3 * D_camera;

% fuel
vars.fuel.a_rf = a_rf;
vars.fuel.n_rf = n_rf;
vars.fuel.rho_f = rho_f;

% combustion
vars.combustion.mdot_ox = mox;
vars.combustion.Tc_fun = T_fun_of_p;
vars.combustion.R_fun = R_fun_of_p;
vars.combustion.k_fun = k_fun_of_p;

vars.combustion.OFmin = OFmin;
vars.combustion.OFmax = OFmax;

vars.combustion.pcmin_bar = pmin * 1e-5;
vars.combustion.pcmax_bar = pmax * 1e-5;


%% ============================================================
%  STEADY STATE INITIAL PRESSURE
% ============================================================

Ap     = vars.geometry.port_area;
perimB = vars.geometry.burning_perimeter;
L      = vars.geometry.grain_length;

% Gox_test = mox / Ap;
% rf_test = a_rf * Gox_test^n_rf;
% Ab_test = perimB * L;
% mdot_f_test = rho_f * Ab_test * rf_test;
% OF_test = mox / mdot_f_test;
% 
% fprintf("\nDEBUG PRE-FZERO:\n");
% fprintf("Ap      = %.6e m^2\n", Ap);
% fprintf("perimB  = %.6e m\n", perimB);
% fprintf("Gox     = %.6f kg/m2s\n", Gox_test);
% fprintf("rf      = %.6e m/s = %.6f mm/s\n", rf_test, rf_test*1e3);
% fprintf("mdot_f  = %.6e kg/s\n", mdot_f_test);
% fprintf("O/F     = %.6f\n", OF_test);

p0 = fzero(@(pc) Z_chamber_stst(pc, vars), [pmin, pmax]);

[~, properties0] = Z_chamber_stst(p0, vars);

OF0 = properties0.O_F;
gox0 = properties0.Gox;

fprintf("\nInitial values:\n");
fprintf("\tpressure = %.1f bar\n", p0*1e-5);
fprintf("\tGOx      = %.1f kg/m2s\n", gox0);
fprintf("\tO/F      = %.2f\n", OF0);


%% ============================================================
%  STEP 3: INTEGRATE COUPLED CHAMBER + MESH SYSTEM
% ============================================================

fine_ode_boolean = false;

plot_case = "refined";
% plot_case = "ode45";

Y0 = [coord_mesh(:,1); coord_mesh(:,2)];

y0 = [p0; Y0(:)];

options = odeset( ...
    'RelTol', 1e-8, ...
    'AbsTol', 1e-10, ...
    'Events', @(t,y) chamber_full_event(t, y, vars, plot_case));

t0 = 0;
tf = tmax;

% Forzo output temporali regolari, così i plot mostrano davvero più istanti
n_output = 10;
tspan = linspace(t0, tf, n_output);

odefun = @(t,y) ode_coupled( ...
    t, ...
    y, ...
    vars, ...
    fine_ode_boolean, ...
    plot_case);

[t_ode, y_in_time, t_event, y_event, i_event] = ode15s(odefun, tspan, y0, options);

if ~isempty(t_event)

    fprintf("\nSimulazione arrestata: camera raggiunta.\n");
    fprintf("Tempo evento = %.6f s\n", t_event(end));

end

p_time = y_in_time(:,1);      % [Pa]
Y_ode  = y_in_time(:,2:end);  % mesh grezza nel tempo


%% ============================================================
%  RAFFINAMENTO A POSTERIORI DELLE MESH SALVATE
% ============================================================
% ode45 salva la mesh grezza.
%
% Qui ricostruisco la mesh raffinata a posteriori, usando la stessa logica
% geometrica della refine. Questa è la mesh da plottare se vuoi vedere
% l'evoluzione corretta rispetto alla geometria raffinata.
%
% In questa versione il raffinamento camera è lasciato commentato:
% se vuoi includere anche il vincolo della camera nel post-processing,
% decommenta le righe indicate dentro il ciclo.

n_steps = length(t_ode);

Y_ode_refined = zeros(size(Y_ode));

for k = 1:n_steps

    Yk_raw = Y_ode(k,:)';

    Pk_raw = stato_to_punti_locale(Yk_raw);

    % Raffinamento geometrico della stella

    Pk_ref = refine_mesh_v3(Pk_raw, idx_v_interni);

    % ------------------------------------------------------------
    % RAFFINAMENTO CAMERA OPZIONALE
    % ------------------------------------------------------------
    % Se vuoi visualizzare anche il vincolo della camera nel post-processing,
    % decommenta la riga seguente.
    %
     [Pk_ref, ~] = refine_mesh_camera(Pk_ref, D_camera);
    %
    % ------------------------------------------------------------

    Y_ode_refined(k,:) = punti_to_stato_locale(Pk_ref)';

end


%% ============================================================
%  BLOCCO OPZIONALE COMPLETO: POST-PROCESSING CAMERA
% ============================================================
% Questo blocco è completamente commentato.
%
% Serve se vuoi costruire una seconda versione della mesh post-processata,
% nella quale:
%
%   1. viene applicata refine_mesh_v3;
%   2. ogni punto che raggiunge/supera la camera viene proiettato
%      esattamente sulla camera.
%
% Per usarlo:
%
%   1. decommenta tutto il blocco;
%   2. nella selezione dati per il plot imposta:
%
%          Y_plot = Y_ode_refined_camera;
%
%      invece di:
%
%          Y_plot = Y_ode_refined;
%
% ------------------------------------------------------------
%
% Y_ode_refined_camera = zeros(size(Y_ode));
%
% for k = 1:n_steps
%
%     Yk_raw = Y_ode(k,:)';
%
%     Pk_raw = stato_to_punti_locale(Yk_raw);
%
%     % Raffinamento geometrico della stella
%     Pk_ref = refine_mesh_v3(Pk_raw, idx_v_interni);
%
%     % Proiezione sulla camera
%     [Pk_ref, idx_wall] = refine_mesh_camera(Pk_ref, D_camera);
%
%     Y_ode_refined_camera(k,:) = punti_to_stato_locale(Pk_ref)';
%
%     if any(idx_wall)
%         fprintf("Camera raggiunta al passo k = %d, t = %.6f s, punti corretti = %d\n", ...
%             k, t_ode(k), nnz(idx_wall));
%     end
%
% end
%
%% ============================================================


%% ============================================================
%  DEBUG RAFFINAMENTO
% ============================================================

k_debug = round(length(t_ode)/2);

Yk_raw = Y_ode(k_debug,:)';
Yk_ref = Y_ode_refined(k_debug,:)';

fprintf("Norma correzione refine al passo k = %d: %.6e\n", ...
    k_debug, norm(Yk_ref - Yk_raw));


%% ============================================================
%  SELEZIONE DATI PER PLOT
% ============================================================

switch plot_case

    case "refined"

        Y_plot = Y_ode_refined;
        t_plot = t_ode;
        plot_label = "mesh raffinata a posteriori";

        % ------------------------------------------------------------
        % Se hai decommentato il blocco camera sopra, puoi usare:
        %
        % Y_plot = Y_ode_refined_camera;
        % plot_label = "mesh raffinata a posteriori con vincolo camera";
        %
        % ------------------------------------------------------------

    case "ode45"

        Y_plot = Y_ode;
        t_plot = t_ode;
        plot_label = "mesh grezza ode45";

    otherwise

        error('plot_case non valido. Usa "refined" oppure "ode45".');

end

n_common = min(size(Y_plot,1), length(t_plot));

Y_plot = Y_plot(1:n_common,:);
t_plot = t_plot(1:n_common);
p_plot = p_time(1:n_common);

fprintf("\nNumero punti usati per il plot = %d\n", n_common);


%% ============================================================
%  CALCOLO AREA, PERIMETRO, GOX, RF, MDOT_F E O/F NEL TEMPO
% ============================================================

n_steps_geom = size(Y_plot, 1);

area_t_m2 = NaN(n_steps_geom, 1);
area_t_mm2 = NaN(n_steps_geom, 1);

perimetro_t_m = NaN(n_steps_geom, 1);
perimetro_t_mm = NaN(n_steps_geom, 1);

Gox_t = NaN(n_steps_geom, 1);
rf_t = NaN(n_steps_geom, 1);

Ab_t = NaN(n_steps_geom, 1);
mdot_f_t = NaN(n_steps_geom, 1);
mdot_in_t = NaN(n_steps_geom, 1);

OF_t = NaN(n_steps_geom, 1);

L      = vars.geometry.grain_length;
mox    = vars.combustion.mdot_ox;
a_rf   = vars.fuel.a_rf;
n_rf   = vars.fuel.n_rf;
rho_f  = vars.fuel.rho_f;

for k = 1:n_steps_geom

    Yk = Y_plot(k,:)';
    Pk = stato_to_punti_locale(Yk);

    [perimetro_t_m(k), area_t_m2(k)] = eval_mesh(Pk, "cartesian");

    perimetro_t_mm(k) = perimetro_t_m(k) * 1e3;
    area_t_mm2(k) = area_t_m2(k) * 1e6;

    if area_t_m2(k) > 0 && perimetro_t_m(k) > 0

        % Flusso ossidante
        Gox_t(k) = mox / area_t_m2(k);

        % Velocità di regressione
        rf_t(k) = a_rf * Gox_t(k)^n_rf;

        % Area bruciante
        Ab_t(k) = perimetro_t_m(k) * L;

        % Portata di fuel
        mdot_f_t(k) = rho_f * Ab_t(k) * rf_t(k);

        % Portata totale entrante
        mdot_in_t(k) = mdot_f_t(k) + mox;

        % Rapporto ossidante/fuel
        if mdot_f_t(k) > 0
            OF_t(k) = mox / mdot_f_t(k);
        else
            OF_t(k) = NaN;
        end

    end

end


%% ============================================================
%  PLOT PRESSIONE CAMERA
% ============================================================

figure;
hold on;
grid on;

plot(t_plot, p_plot*1e-5, ...
    'LineWidth', 1.8, ...
    'DisplayName', 'p_c');

xlabel('Tempo [s]');
ylabel('Pressione camera [bar]');
title("Evoluzione pressione camera - caso " + plot_label);
legend('Location', 'best');


%% ============================================================
%  PLOT EVOLUZIONE AREA E PERIMETRO
% ============================================================

figure;
hold on;
grid on;

yyaxis left
plot(t_plot, area_t_mm2, ...
    'LineWidth', 1.8, ...
    'DisplayName', 'Area');
ylabel('Area [mm^2]');

yyaxis right
plot(t_plot, perimetro_t_mm, ...
    'LineWidth', 1.8, ...
    'DisplayName', 'Perimetro');
ylabel('Perimetro [mm]');

xlabel('Tempo [s]');
title("Evoluzione area e perimetro - caso " + plot_label);
legend('Location', 'best');


%% ============================================================
%  PLOT EVOLUZIONE GOX
% ============================================================

figure;
hold on;
grid on;

plot(t_plot, Gox_t, ...
    'LineWidth', 1.8, ...
    'DisplayName', 'G_{ox}');

xlabel('Tempo [s]');
ylabel('G_{ox} [kg/(m^2 s)]');
title("Evoluzione G_{ox} - caso " + plot_label);
legend('Location', 'best');


%% ============================================================
%  PLOT EVOLUZIONE O/F
% ============================================================

figure;
hold on;
grid on;

plot(t_plot, OF_t, ...
    'LineWidth', 1.8, ...
    'DisplayName', 'O/F');

if isfield(vars.combustion, "OFmin")
    yline(vars.combustion.OFmin, 'k--', ...
        'LineWidth', 1.2, ...
        'DisplayName', 'O/F min CEA');
end

if isfield(vars.combustion, "OFmax")
    yline(vars.combustion.OFmax, 'k--', ...
        'LineWidth', 1.2, ...
        'DisplayName', 'O/F max CEA');
end

xlabel('Tempo [s]');
ylabel('O/F [-]');
title("Evoluzione O/F - caso " + plot_label);
legend('Location', 'best');


%% ============================================================
%  PLOT CONFRONTO GOX E O/F
% ============================================================

figure;
hold on;
grid on;

yyaxis left
plot(t_plot, Gox_t, ...
    'LineWidth', 1.8, ...
    'DisplayName', 'G_{ox}');
ylabel('G_{ox} [kg/(m^2 s)]');

yyaxis right
plot(t_plot, OF_t, ...
    'LineWidth', 1.8, ...
    'DisplayName', 'O/F');
ylabel('O/F [-]');

xlabel('Tempo [s]');
title("Confronto G_{ox} e O/F - caso " + plot_label);
legend('Location', 'best');


%% ============================================================
%  PLOT EVOLUZIONE PORTATE
% ============================================================

figure;
hold on;
grid on;

plot(t_plot, mox*ones(size(t_plot)), ...
    'LineWidth', 1.8, ...
    'DisplayName', '\dot{m}_{ox}');

plot(t_plot, mdot_f_t, ...
    'LineWidth', 1.8, ...
    'DisplayName', '\dot{m}_{f}');

plot(t_plot, mdot_in_t, ...
    'LineWidth', 1.8, ...
    'DisplayName', '\dot{m}_{in}');

xlabel('Tempo [s]');
ylabel('Portata [kg/s]');
title("Evoluzione portate - caso " + plot_label);
legend('Location', 'best');

%% ============================================================
%  PLOT VARIAZIONE GLOBALE MESH
% ============================================================

delta_mesh = zeros(n_steps_geom,1);

Y_initial = Y_plot(1,:)';

for k = 1:n_steps_geom

    Yk = Y_plot(k,:)';
    delta_mesh(k) = norm(Yk - Y_initial);

end

figure;
hold on;
grid on;

plot(t_plot, delta_mesh*1e3, ...
    'LineWidth', 1.8, ...
    'DisplayName', '||Y(t)-Y(0)||');

xlabel('Tempo [s]');
ylabel('Norma variazione mesh [mm]');
title('Variazione globale della mesh nel tempo');
legend('Location', 'best');


%% ============================================================
%  PLOT CONFRONTO - MESH NEL TEMPO
% ============================================================

figure;
hold on;
grid on;
axis equal;

n_steps_plot = size(Y_plot, 1);

n_plot = 10;
n_plot = min(n_plot, n_steps_plot);

idx_plot = round(linspace(1, n_steps_plot, n_plot));
idx_plot = unique(idx_plot, 'stable');

fprintf("\nIndici plottati:\n");
disp(idx_plot);

fprintf("Tempi plottati:\n");
disp(t_plot(idx_plot).');

colors = lines(length(idx_plot));

all_x = [];
all_y = [];

for jj = 1:length(idx_plot)

    k = idx_plot(jj);

    Yk = Y_plot(k,:)';
    Pk = stato_to_punti_locale(Yk) * 1e3;

    all_x = [all_x; Pk(:,1)];
    all_y = [all_y; Pk(:,2)];

end

for jj = 1:length(idx_plot)

    k = idx_plot(jj);

    Yk = Y_plot(k,:)';
    Pk = stato_to_punti_locale(Yk) * 1e3;

    plot(Pk(:,1), Pk(:,2), ...
        'LineWidth', 1.5, ...
        'Color', colors(jj,:), ...
        'DisplayName', sprintf('t = %.3f s', t_plot(k)));

end

P_iniziale = stato_to_punti_locale(Y_plot(1,:)') * 1e3;
P_finale = stato_to_punti_locale(Y_plot(end,:)') * 1e3;

plot(P_iniziale(:,1), P_iniziale(:,2), ...
    'k--', ...
    'LineWidth', 1.4, ...
    'DisplayName', 'Mesh iniziale');

plot(P_finale(:,1), P_finale(:,2), ...
    'k-', ...
    'LineWidth', 2.0, ...
    'DisplayName', 'Mesh finale');

plot(P_finale(idx_v_interni,1), ...
     P_finale(idx_v_interni,2), ...
     'ro', ...
     'MarkerSize', 8, ...
     'LineWidth', 1.5, ...
     'DisplayName','Punte interne finali');

plot(P_finale(idx_v_esterni,1), ...
     P_finale(idx_v_esterni,2), ...
     'gs', ...
     'MarkerSize', 8, ...
     'LineWidth', 1.5, ...
     'DisplayName','Punte esterne finali');

% Camera in mm
R_camera_mm = D_camera * 1e3 / 2;

theta_cam = linspace(0, pi, 600);

x_camera = R_camera_mm*cos(theta_cam);
y_camera = R_camera_mm*sin(theta_cam);

plot(x_camera, y_camera, 'k--', ...
    'LineWidth', 1.2, ...
    'DisplayName', 'Camera di combustione');

% Assi basati sulla mesh, non sulla camera
xmin = min(all_x);
xmax = max(all_x);
ymin = min(all_y);
ymax = max(all_y);

dx = xmax - xmin;
dy = ymax - ymin;

pad = 0.20 * max(dx, dy);

if pad == 0
    pad = 1;
end

xlim([xmin - pad, xmax + pad]);
ylim([ymin - pad, ymax + pad]);

xlabel('x [mm]');
ylabel('y [mm]');
title("Evoluzione della mesh - caso " + plot_label);
legend('Location','bestoutside');


%% ============================================================
%  OUTPUT FINALE
% ============================================================

fprintf("\nValori finali:\n");
fprintf("Pressione finale = %.6f bar\n", p_plot(end)*1e-5);
fprintf("Area finale = %.6f mm^2\n", area_t_mm2(end));
fprintf("Area finale = %.6e m^2\n", area_t_m2(end));
fprintf("Perimetro finale = %.6f mm\n", perimetro_t_mm(end));
fprintf("Perimetro finale = %.6e m\n", perimetro_t_m(end));
fprintf("Gox finale = %.6f kg/(m^2 s)\n", Gox_t(end));
fprintf("Variazione mesh finale = %.6f mm\n", delta_mesh(end)*1e3);


%% ============================================================
%  ANIMAZIONE EVOLUZIONE TEMPORALE DELLA MESH
% ============================================================

figure;
hold on;
grid on;
axis equal;

xlabel('x [mm]');
ylabel('y [mm]');
title("Animazione regressione mesh - caso " + plot_label);

n_steps_anim = size(Y_plot, 1);

n_frame = 200;

idx_anim = round(linspace(1, n_steps_anim, min(n_frame, n_steps_anim)));
idx_anim = unique(idx_anim, 'stable');

all_x_anim = [];
all_y_anim = [];

for jj = 1:length(idx_anim)

    k = idx_anim(jj);

    Yk = Y_plot(k,:)';
    Pk = stato_to_punti_locale(Yk) * 1e3;

    all_x_anim = [all_x_anim; Pk(:,1)];
    all_y_anim = [all_y_anim; Pk(:,2)];

end

xmin = min(all_x_anim);
xmax = max(all_x_anim);
ymin = min(all_y_anim);
ymax = max(all_y_anim);

dx = xmax - xmin;
dy = ymax - ymin;

pad = 0.20 * max(dx, dy);

if pad == 0
    pad = 1;
end

xlim([xmin - pad, xmax + pad]);
ylim([ymin - pad, ymax + pad]);

plot(x_camera, y_camera, 'k--', ...
    'LineWidth', 1.2, ...
    'DisplayName', 'Camera di combustione');

Y0_plot = Y_plot(1,:)';
P0_plot = stato_to_punti_locale(Y0_plot) * 1e3;

plot(P0_plot(:,1), P0_plot(:,2), 'k:', ...
    'LineWidth', 1.2, ...
    'DisplayName', 'Mesh iniziale');

Y_first = Y_plot(idx_anim(1),:)';
P_first = stato_to_punti_locale(Y_first) * 1e3;

h_mesh = plot(P_first(:,1), P_first(:,2), 'b-', ...
    'LineWidth', 1.8, ...
    'DisplayName', 'Mesh corrente');

h_tips_int = plot(P_first(idx_v_interni,1), ...
                  P_first(idx_v_interni,2), ...
                  'ro', ...
                  'MarkerSize', 7, ...
                  'LineWidth', 1.5, ...
                  'DisplayName', 'Punte interne');

h_tips_ext = plot(P_first(idx_v_esterni,1), ...
                  P_first(idx_v_esterni,2), ...
                  'gs', ...
                  'MarkerSize', 7, ...
                  'LineWidth', 1.5, ...
                  'DisplayName', 'Punte esterne');

h_time = text(0.02, 0.95, '', ...
    'Units', 'normalized', ...
    'FontSize', 12, ...
    'FontWeight', 'bold');

legend('Location', 'bestoutside');

for jj = 1:length(idx_anim)

    k = idx_anim(jj);

    Yk = Y_plot(k,:)';
    Pk = stato_to_punti_locale(Yk) * 1e3;

    set(h_mesh, ...
        'XData', Pk(:,1), ...
        'YData', Pk(:,2));

    set(h_tips_int, ...
        'XData', Pk(idx_v_interni,1), ...
        'YData', Pk(idx_v_interni,2));

    set(h_tips_ext, ...
        'XData', Pk(idx_v_esterni,1), ...
        'YData', Pk(idx_v_esterni,2));

    time_string = sprintf('t = %.3f s | p_c = %.2f bar', ...
        t_plot(k), p_plot(k)*1e-5);

    set(h_time, 'String', time_string);

    drawnow;
    pause(0.02);

end

toc
%% ============================================================
%  FUNZIONE LOCALE: STATO -> PUNTI
% ============================================================

function P = stato_to_punti_locale(Y)

    Y = Y(:);

    n = length(Y)/2;

    x = Y(1:n);
    y = Y(n+1:end);

    P = [x y];

end


%% ============================================================
%  FUNZIONE LOCALE: PUNTI -> STATO
% ============================================================

function Y = punti_to_stato_locale(P)

    x = P(:,1);
    y = P(:,2);

    Y = [x; y];

end


%% ============================================================
%  FUNZIONE LOCALE: PROIEZIONE SULLA CAMERA
% ============================================================
% Questa funzione non viene chiamata finché il blocco camera resta commentato.
% Serve solo se decommenti il post-processing opzionale sulla camera.

function [P_ref, idx_wall] = refine_mesh_camera(P, D_camera)

    P_ref = P;

    R_camera = D_camera / 2;

    r = hypot(P(:,1), P(:,2));

    tol_wall = 1e-12 * R_camera;

    idx_wall = r >= R_camera - tol_wall;

    for i = 1:size(P,1)

        if idx_wall(i) && r(i) > eps

            P_ref(i,:) = R_camera * P(i,:) / r(i);

        end

    end

end


%% ============================================================
%  FUNZIONE LOCALE: EVENTO CAMERA RAGGIUNTA
% ============================================================
% Ferma la simulazione quando tutti i punti della mesh fisica/evaluated
% hanno raggiunto il raggio della camera.
%
% Nota:
%   uso la stessa logica del plot_case:
%   - se plot_case = "refined", prima ricostruisco la mesh raffinata;
%   - se plot_case = "ode45", uso la mesh grezza.
%
% L'evento non controlla il diametro, ma il raggio:
%
%   r_i >= R_camera = D_camera/2

function [value, isterminal, direction] = chamber_full_event(~, y, vars, plot_case)

    Y = y(2:end);
    Y = Y(:);

    P_raw = stato_to_punti_locale(Y);

    D_camera = vars.geometry.D_camera;
    R_camera = D_camera / 2;

    idx_v_interni = vars.geometry.idx_v_interni;

    %% ============================================================
    %  1. RICOSTRUZIONE DELLA MESH DA CONTROLLARE
    % ============================================================

    switch plot_case

        case "refined"

    % Prima correzione sulla mesh grezza
    [P_raw, ~] = refine_mesh_camera(P_raw, D_camera);

    % Raffinamento
    P_check = refine_mesh_v3(P_raw, idx_v_interni);

    % Seconda correzione dopo il raffinamento
    [P_check, ~] = refine_mesh_camera(P_check, D_camera);

        case "ode45"

    [P_check, ~] = refine_mesh_camera(P_raw, D_camera);

        otherwise

            error("plot_case non valido in chamber_full_event.");

    end

    %% ============================================================
    %  2. CONTROLLO RAGGI
    % ============================================================

    r = hypot(P_check(:,1), P_check(:,2));

    % Tolleranza di arresto.
    % Non usarla troppo piccola: se è troppo piccola, ode45 potrebbe
    % non intercettare mai esattamente il bordo.
    if isfield(vars.geometry, "wall_event_tol")
        wall_event_tol = vars.geometry.wall_event_tol;
    else
        wall_event_tol = 1e-3;   % [m] = 0.001 mm
    end

    % Voglio fermarmi quando anche il punto più interno è arrivato
    % abbastanza vicino alla camera.
    %
    % value < 0 : non tutti i punti sono arrivati
    % value = 0 : tutti i punti sono arrivati entro tolleranza
    % value > 0 : tutti i punti hanno raggiunto/superato la camera

    value = min(r) - (R_camera - wall_event_tol);

    % Ferma integrazione
    isterminal = 1;

    % Più robusto di +1, perché con correzioni geometriche/refine
    % l'attraversamento può non essere monotono perfetto.
    direction = 0;

end