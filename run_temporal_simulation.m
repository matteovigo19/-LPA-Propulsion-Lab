function sim = run_temporal_simulation(pre, settings)
% RUN_TEMPORAL_SIMULATION
% Esegue la simulazione temporale accoppiata pressione + regressione mesh.
%
% INPUT
%   pre       output di run_predesign
%
%   settings  struct opzionale con campi:
%       .tmax               [s]    default 300
%       .ext_diameter       [m]    default 280.1*2e-3
%       .multiplyer         [-]    default 80 per star
%       .n_cylinder         [-]    default 500
%       .fine_ode_boolean   logical default false
%       .plot_case          "ode45" oppure "refined"; default "ode45"
%       .RelTol             default 1e-8
%       .AbsTol             default 1e-10
%       .n_output           default 8
%       .make_plots         logical default false
%       .make_animation     logical default false
%       .verbose            logical default false
%
% OUTPUT
%   sim       struct con risultati temporali, geometrici e di combustione
%
% NOTA
%   Per l'ottimizzazione conviene:
%       settings.make_plots = false;
%       settings.make_animation = false;
%       settings.verbose = false;

if nargin < 2
    settings = struct();
end

%% ============================================================
%  SETTINGS
% ============================================================

settings = set_default(settings, "tmax", 300);
settings = set_default(settings, "ext_diameter", 280.1*2e-3);
settings = set_default(settings, "multiplyer", 80);
settings = set_default(settings, "n_cylinder", 500);
settings = set_default(settings, "fine_ode_boolean", false);
settings = set_default(settings, "plot_case", "ode45");
settings = set_default(settings, "RelTol", 1e-8);
settings = set_default(settings, "AbsTol", 1e-10);
settings = set_default(settings, "n_output", 8);
settings = set_default(settings, "make_plots", false);
settings = set_default(settings, "make_animation", false);
settings = set_default(settings, "verbose", false);

type = string(pre.type);

%% ============================================================
%  STEP 0: DATA
% ============================================================

vars = struct();

% ---------------- GEOMETRY DATA ----------------

if type == "star"

    n = settings.multiplyer * pre.n_tips;

    vars.geometry.ntips = pre.n_tips;
    vars.geometry.diametro_est = pre.diam_e;
    vars.geometry.diametro_int = pre.diam_i;

elseif type == "cylinder"

    n = settings.n_cylinder;
    vars.geometry.diametro = pre.diam;

else

    error("Tipo non riconosciuto nella simulazione temporale: %s", type);

end

vars.geometry.cyl_tol_rel = 0.03;
vars.geometry.cyl_tol_abs = 0.15e-3;

% ---------------- FUEL DATA ----------------

rho_f = pre.rho_f;
a_rf  = pre.a_rf;
n_rf  = pre.n_rf;

% ---------------- COMBUSTION DATA ----------------

mox = pre.mox;

if isfield(pre, "CEA")
    T_fun_of_p = pre.CEA.T_fun_of_p;
    k_fun_of_p = pre.CEA.k_fun_of_p;
    R_fun_of_p = pre.CEA.R_fun_of_p;
else
    cea = load("CEA_functions.mat");
    T_fun_of_p = cea.T_fun_of_p;
    k_fun_of_p = cea.k_fun_of_p;
    R_fun_of_p = cea.R_fun_of_p;
end

OFmin = 1.1;
OFmax = 19.5;

pmin = 101325;
pmax = 99e5;

tmax = settings.tmax;
pamb = 0; 

% ---------------- ENGINE DATA ----------------

chamber_length = pre.L;
D_camera = settings.ext_diameter;

%% ============================================================
%  STEP 1: INITIAL MESH
% ============================================================

[coordstar_cart, coord_mesh, coord_v, idx_v] = ...
    make_mesh0(type, n, vars, "cartesiane"); %#ok<ASGLU>

r_v = vecnorm(coord_mesh(idx_v,:), 2, 2);

idx_v_interni = [];
idx_v_esterni = [];

if type == "star"

    r_thr = 0.5*(min(r_v) + max(r_v));

    idx_v_interni = idx_v(r_v < r_thr);
    idx_v_esterni = idx_v(r_v >= r_thr);

    n_mesh = size(coord_mesh,1);

    if max(idx_v) > n_mesh || ...
       max(idx_v_interni) > n_mesh || ...
       max(idx_v_esterni) > n_mesh

        error("Indici punte non coerenti con coord_mesh.");

    end
end

%% ============================================================
%  STEP 2: INITIAL GEOMETRY / STEADY STATE
% ============================================================

[perim0, Ap0] = eval_mesh(coord_mesh, "cartesian");

tol = 1e-4;

if abs(pre.Ap - Ap0) > tol
    error("L'area di porto del temporale e del predesign non corrispondono. pre.Ap=%.6e, Ap0=%.6e", ...
        pre.Ap, Ap0);
end

Gox_test = mox / Ap0;
rf_test = a_rf * Gox_test^n_rf;

pc_target_bar = pre.pch_bar;
pc_target = pc_target_bar * 1e5;

% Initial fuel flow from current mesh
Gox0 = mox / Ap0;
rf0 = a_rf * Gox0^n_rf;
Ab0 = perim0 * chamber_length;
mdot_f0 = rho_f * Ab0 * rf0;

mdot_tot = mox + mdot_f0;
OF0_geom = mox / mdot_f0;

Tc0 = T_fun_of_p(OF0_geom, pc_target_bar);
k0  = k_fun_of_p(OF0_geom, pc_target_bar);
R0  = R_fun_of_p(OF0_geom, pc_target_bar);

K2 = k0 * ((2/(k0+1))^((k0+1)/(k0-1)));
cstar0 = sqrt(R0*Tc0/K2);

At = (mdot_tot*cstar0)/pc_target;
throat_diameter = 2*sqrt(At/pi);

%% ============================================================
%  COMPACT DATA IN VARS
% ============================================================

vars.geometry.port_area = Ap0;
vars.geometry.burning_perimeter = perim0;
vars.geometry.grain_length = chamber_length;
vars.geometry.throat_area = At;
vars.geometry.type = type;
vars.geometry.D_camera = D_camera;
vars.geometry.wall_tol = 1e-3*D_camera;

if isfield(settings, "wall_event_tol")
    vars.geometry.wall_event_tol = settings.wall_event_tol;
end

if type == "star"
    vars.geometry.idx_v = idx_v;
    vars.geometry.idx_v_interni = idx_v_interni;
    vars.geometry.idx_v_esterni = idx_v_esterni;
end

vars.fuel.a_rf = a_rf;
vars.fuel.n_rf = n_rf;
vars.fuel.rho_f = rho_f;

vars.combustion.mdot_ox = mox;
vars.combustion.Tc_fun = T_fun_of_p;
vars.combustion.R_fun = R_fun_of_p;
vars.combustion.k_fun = k_fun_of_p;

vars.combustion.OFmin = OFmin;
vars.combustion.OFmax = OFmax;
vars.combustion.pcmin_bar = pmin*1e-5;
vars.combustion.pcmax_bar = pmax*1e-5;

%% ============================================================
%  INITIAL STEADY-STATE PRESSURE
% ============================================================

p0 = fzero(@(pc) Z_chamber_stst(pc, vars), [pmin, pmax]);

[~, properties0] = Z_chamber_stst(p0, vars);

OF0 = properties0.O_F;
gox0 = properties0.Gox;

%% ============================================================
%  STEP 3: COUPLED ODE
% ============================================================

fine_ode_boolean = settings.fine_ode_boolean;
plot_case = string(settings.plot_case);

Y0 = [coord_mesh(:,1); coord_mesh(:,2)];
y0 = [p0; Y0(:)];

options = odeset( ...
    'RelTol', settings.RelTol, ...
    'AbsTol', settings.AbsTol, ...
    'Events', @(t,y) chamber_full_event_local(t, y, vars, plot_case, type));

tspan = linspace(0, tmax, settings.n_output);

odefun = @(t,y) ode_coupled( ...
    t, ...
    y, ...
    vars, ...
    fine_ode_boolean, ...
    plot_case);

[t_ode, y_in_time, t_event, y_event, i_event] = ...
    ode15s(odefun, tspan, y0, options);

p_time = y_in_time(:,1);
Y_ode = y_in_time(:,2:end);

%% ============================================================
%  OPTIONAL POSTERIOR REFINEMENT
% ============================================================

if plot_case == "refined"

    n_steps = length(t_ode);
    Y_ode_refined = zeros(size(Y_ode));

    for k = 1:n_steps

        Yk_raw = Y_ode(k,:)';
        Pk_raw = stato_to_punti_locale(Yk_raw);

        if type == "star"
            Pk_ref = refine_mesh_v3(Pk_raw, idx_v_interni);
        else
            Pk_ref = Pk_raw;
        end

        [Pk_ref, ~] = refine_mesh_camera_local(Pk_ref, D_camera);

        Y_ode_refined(k,:) = punti_to_stato_locale(Pk_ref)';

    end

    Y_plot = Y_ode_refined;
    plot_label = "mesh raffinata a posteriori";

elseif plot_case == "ode45"

    Y_plot = Y_ode;
    plot_label = "mesh grezza ode45";

else

    error('plot_case non valido. Usa "refined" oppure "ode45".');

end

t_plot = t_ode;

n_common = min(size(Y_plot,1), length(t_plot));

Y_plot = Y_plot(1:n_common,:);
t_plot = t_plot(1:n_common);
p_plot = p_time(1:n_common);

%% ============================================================
%  TIME-HISTORY POST-PROCESSING
% ============================================================

n_steps_geom = size(Y_plot,1);

area_t_m2 = NaN(n_steps_geom,1);
area_t_mm2 = NaN(n_steps_geom,1);

perimetro_t_m = NaN(n_steps_geom,1);
perimetro_t_mm = NaN(n_steps_geom,1);

Gox_t = NaN(n_steps_geom,1);
rf_t = NaN(n_steps_geom,1);

Ab_t = NaN(n_steps_geom,1);
mdot_f_t = NaN(n_steps_geom,1);
mdot_in_t = NaN(n_steps_geom,1);
OF_t = NaN(n_steps_geom,1);

L = vars.geometry.grain_length;

for k = 1:n_steps_geom

    Yk = Y_plot(k,:)';
    Pk = stato_to_punti_locale(Yk);

    [perimetro_t_m(k), area_t_m2(k)] = eval_mesh(Pk, "cartesian");

    perimetro_t_mm(k) = perimetro_t_m(k)*1e3;
    area_t_mm2(k) = area_t_m2(k)*1e6;

    if area_t_m2(k) > 0 && perimetro_t_m(k) > 0

        Gox_t(k) = mox/area_t_m2(k);

        rf_t(k) = a_rf*Gox_t(k)^n_rf;

        Ab_t(k) = perimetro_t_m(k)*L;

        mdot_f_t(k) = rho_f*Ab_t(k)*rf_t(k);

        mdot_in_t(k) = mdot_f_t(k) + mox;

        if mdot_f_t(k) > 0
            OF_t(k) = mox/mdot_f_t(k);
        end

    end
end

% Global mesh variation
delta_mesh = zeros(n_steps_geom,1);
Y_initial = Y_plot(1,:)';

for k = 1:n_steps_geom
    delta_mesh(k) = norm(Y_plot(k,:)' - Y_initial);
end

%% ============================================================
%  OUTPUT STRUCT
% ============================================================

sim.pre = pre;
sim.settings = settings;
sim.vars = vars;

sim.type = type;
sim.plot_case = plot_case;
sim.plot_label = plot_label;

sim.coord_mesh0 = coord_mesh;
sim.coord_v = coord_v;
sim.idx_v = idx_v;
sim.idx_v_interni = idx_v_interni;
sim.idx_v_esterni = idx_v_esterni;

sim.perim0 = perim0;
sim.Ap0 = Ap0;
sim.Gox0_geometry = Gox_test;
sim.rf0 = rf_test;

sim.At = At;
sim.throat_diameter = throat_diameter;

sim.p0 = p0;
sim.OF0 = OF0;
sim.Gox0 = gox0;

sim.t = t_plot;
sim.p = p_plot;
sim.Y = Y_plot;
sim.Y_raw = Y_ode;

sim.area_m2 = area_t_m2;
sim.area_mm2 = area_t_mm2;
sim.perimeter_m = perimetro_t_m;
sim.perimeter_mm = perimetro_t_mm;

sim.Gox = Gox_t;
sim.rf = rf_t;
sim.Ab = Ab_t;
sim.mdot_f = mdot_f_t;
sim.mdot_ox = mox*ones(size(t_plot));
sim.mdot_in = mdot_in_t;
sim.OF = OF_t;

sim.delta_mesh = delta_mesh;

sim.t_event = t_event;
sim.y_event = y_event;
sim.i_event = i_event;
sim.event_occurred = ~isempty(t_event);

%% ============================================================
%  VERBOSE OUTPUT
% ============================================================

if settings.verbose

    fprintf("\nGeometria iniziale:\n");
    fprintf("perimeter = %.6e m\n", perim0);
    fprintf("area      = %.6e m^2\n", Ap0);

    fprintf("\nInitial values:\n");
    fprintf("\tpressure = %.3f bar\n", p0*1e-5);
    fprintf("\tGOx      = %.3f kg/m2s\n", gox0);
    fprintf("\tO/F      = %.4f\n", OF0);

    if sim.event_occurred
        fprintf("\nSimulazione arrestata per evento camera.\n");
        fprintf("Tempo evento = %.6f s\n", t_event(end));
    end

    fprintf("\nValori finali:\n");
    fprintf("Pressione finale = %.6f bar\n", p_plot(end)*1e-5);
    fprintf("Area finale = %.6e m^2\n", area_t_m2(end));
    fprintf("Perimetro finale = %.6e m\n", perimetro_t_m(end));
    fprintf("Gox finale = %.6f kg/(m^2 s)\n", Gox_t(end));
    fprintf("O/F finale = %.6f\n", OF_t(end));

end

%% ============================================================
%  OPTIONAL PLOTS
% ============================================================

if settings.make_plots
    plot_temporal_results(sim);
end

if settings.make_animation
    animate_temporal_mesh(sim);
end

end


%% ============================================================
%  LOCAL HELPERS
% ============================================================

function s = set_default(s, fieldname, value)

if ~isfield(s, fieldname)
    s.(fieldname) = value;
end

end


function P = stato_to_punti_locale(Y)

Y = Y(:);

n = length(Y)/2;

x = Y(1:n);
y = Y(n+1:end);

P = [x y];

end


function Y = punti_to_stato_locale(P)

Y = [P(:,1); P(:,2)];

end


function [P_ref, idx_wall] = refine_mesh_camera_local(P, D_camera)

P_ref = P;

R_camera = D_camera/2;

r = hypot(P(:,1), P(:,2));

tol_wall = 1e-12*R_camera;

idx_wall = r >= R_camera - tol_wall;

for i = 1:size(P,1)

    if idx_wall(i) && r(i) > eps
        P_ref(i,:) = R_camera*P(i,:)/r(i);
    end

end

end


function [value, isterminal, direction] = ...
    chamber_full_event_local(~, y, vars, plot_case, type)

Y = y(2:end);
P_raw = stato_to_punti_locale(Y);

D_camera = vars.geometry.D_camera;
R_camera = D_camera/2;

switch plot_case

    case "refined"

        [P_raw, ~] = refine_mesh_camera_local(P_raw, D_camera);

        if type == "star"
            P_check = refine_mesh_v3(P_raw, vars.geometry.idx_v_interni);
        else
            P_check = P_raw;
        end

        [P_check, ~] = refine_mesh_camera_local(P_check, D_camera);

    case "ode45"

        [P_check, ~] = refine_mesh_camera_local(P_raw, D_camera);

    otherwise

        error("plot_case non valido in chamber_full_event.");

end

r = hypot(P_check(:,1), P_check(:,2));

if isfield(vars.geometry, "wall_event_tol")
    wall_event_tol = vars.geometry.wall_event_tol;
else
    % Mantiene il valore numerico del main originale: 1e-3 m = 1 mm.
    wall_event_tol = 1e-3;
end

value = min(r) - (R_camera - wall_event_tol);

isterminal = 1;
direction = 0;

end


function plot_temporal_results(sim)

t = sim.t;

figure;
plot(t, sim.p*1e-5, 'LineWidth', 1.8);
grid on;
xlabel('Tempo [s]');
ylabel('Pressione camera [bar]');
title("Evoluzione pressione camera - " + sim.plot_label);

figure;
yyaxis left
plot(t, sim.area_mm2, 'LineWidth', 1.8);
ylabel('Area [mm^2]');
yyaxis right
plot(t, sim.perimeter_mm, 'LineWidth', 1.8);
ylabel('Perimetro [mm]');
grid on;
xlabel('Tempo [s]');
title("Evoluzione area e perimetro - " + sim.plot_label);

figure;
plot(t, sim.Gox, 'LineWidth', 1.8);
grid on;
xlabel('Tempo [s]');
ylabel('G_{ox} [kg/(m^2 s)]');
title("Evoluzione G_{ox}");

figure;
plot(t, sim.OF, 'LineWidth', 1.8);
grid on;
hold on;
if isfield(sim.vars.combustion, "OFmin")
    yline(sim.vars.combustion.OFmin, 'k--');
end
if isfield(sim.vars.combustion, "OFmax")
    yline(sim.vars.combustion.OFmax, 'k--');
end
xlabel('Tempo [s]');
ylabel('O/F [-]');
title("Evoluzione O/F");

figure;
plot(t, sim.mdot_ox, 'LineWidth', 1.8, ...
    'DisplayName', '\dot{m}_{ox}');
hold on;
plot(t, sim.mdot_f, 'LineWidth', 1.8, ...
    'DisplayName', '\dot{m}_{f}');
plot(t, sim.mdot_in, 'LineWidth', 1.8, ...
    'DisplayName', '\dot{m}_{in}');
grid on;
xlabel('Tempo [s]');
ylabel('Portata [kg/s]');
legend('Location','best');
title("Evoluzione portate");

figure;
plot(t, sim.delta_mesh*1e3, 'LineWidth', 1.8);
grid on;
xlabel('Tempo [s]');
ylabel('Norma variazione mesh [mm]');
title('Variazione globale della mesh');

% Mesh evolution
figure;
hold on;
grid on;
axis equal;

n_steps = size(sim.Y,1);
n_plot = min(10,n_steps);
idx_plot = unique(round(linspace(1,n_steps,n_plot)), 'stable');

for jj = 1:length(idx_plot)

    k = idx_plot(jj);

    Pk = stato_to_punti_locale(sim.Y(k,:)')*1e3;

    plot(Pk(:,1), Pk(:,2), ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('t = %.3f s', t(k)));

end

R_camera_mm = sim.vars.geometry.D_camera*1e3/2;
theta_cam = linspace(0,pi,600);

plot(R_camera_mm*cos(theta_cam), ...
     R_camera_mm*sin(theta_cam), ...
     'k--', ...
     'DisplayName','Camera');

xlabel('x [mm]');
ylabel('y [mm]');
title("Evoluzione mesh - " + sim.plot_label);
legend('Location','bestoutside');

end


function animate_temporal_mesh(sim)

figure;
hold on;
grid on;
axis equal;

xlabel('x [mm]');
ylabel('y [mm]');
title("Animazione regressione mesh - " + sim.plot_label);

n_steps = size(sim.Y,1);
idx_anim = unique(round(linspace(1,n_steps,min(200,n_steps))), 'stable');

% Limits from all selected frames
all_x = [];
all_y = [];

for jj = 1:length(idx_anim)

    Pk = stato_to_punti_locale(sim.Y(idx_anim(jj),:)')*1e3;

    all_x = [all_x; Pk(:,1)]; %#ok<AGROW>
    all_y = [all_y; Pk(:,2)]; %#ok<AGROW>

end

xmin = min(all_x);
xmax = max(all_x);
ymin = min(all_y);
ymax = max(all_y);

pad = 0.20*max(xmax-xmin, ymax-ymin);

if pad == 0
    pad = 1;
end

xlim([xmin-pad, xmax+pad]);
ylim([ymin-pad, ymax+pad]);

R_camera_mm = sim.vars.geometry.D_camera*1e3/2;
theta_cam = linspace(0,pi,600);

plot(R_camera_mm*cos(theta_cam), ...
     R_camera_mm*sin(theta_cam), ...
     'k--');

P_first = stato_to_punti_locale(sim.Y(idx_anim(1),:)')*1e3;

h_mesh = plot(P_first(:,1), P_first(:,2), 'b-', 'LineWidth',1.8);

h_time = text(0.02,0.95,'', ...
    'Units','normalized', ...
    'FontSize',12, ...
    'FontWeight','bold');

for jj = 1:length(idx_anim)

    k = idx_anim(jj);

    Pk = stato_to_punti_locale(sim.Y(k,:)')*1e3;

    set(h_mesh, 'XData',Pk(:,1), 'YData',Pk(:,2));

    set(h_time, 'String', ...
        sprintf('t = %.3f s | p_c = %.2f bar', ...
        sim.t(k), sim.p(k)*1e-5));

    drawnow;
    pause(0.02);

end

end
