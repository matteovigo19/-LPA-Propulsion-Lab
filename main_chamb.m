close all;
clc;


%% STEP 0: data

% GEOMETRY DATA

% % Inner port data:
% type = "cylinder";
% refine_mesh_boolean = false;
% vars.diameter = 0.03;       % [m]
% npoints = 100;

% Inner port data:
type="star";
%refine_mesh_boolean = true;
multiplyer = 100;
vars.geometry.ntips = 8;
vars.geometry.diametro_est = 3;
vars.geometry.diametro_int = 1;
n = multiplyer * vars.geometry.ntips;

% ========== FUEL DATA ===================
% FUEL: HTPB
% regression rate: rf = a*GOx^n
rho_f = 920;                  % [kg/m3]
a_rf = 0.027;                 % [(mm/s)/((kg/m2 s)^n)]
n_rf = 0.75;                  % [-]
a_rf = a_rf*1e-3;             % [m/s]

% ======== COMBUSTION DATA ===============
% OXIDIZER: O2
mox = 0.015;                  % [kg/s]
load("CEA_functions.mat");
% ! Limits
OFmin = 1.1;  % 1
OFmax = 19.5; % 20
pmax = 99e5;  % 100 bar    % [Pa]
pmin = 101325; % 1 bar     % [Pa]
tmax = 20;        % maximum burning time [s]
pamb = 101325;    % [Pa]

fine_ode_boolean = true;

% =========== ENGINE DATA =================
ext_diameter = 0.05;           % [m]
chamber_length = 0.09;         % [m]
throat_diameter = 0.005;       % [m]
eps = 2;
At = 0.25*pi*throat_diameter^2; % [m2]

% =========== end of data =================



% Chiamate alla tua funzione grainmesh
[coordstar_cart, coord_mesh, coord_v, idx_v, idx_v_interni, idx_v_esterni] = ...
    make_mesh0("star", n, vars, "cartesiane");

%% OUTPUT MESH INIZIALE

figure()
plot(coord_mesh(:,1), coord_mesh(:,2), 'b.-')
hold on
plot(coord_mesh(idx_v,1), coord_mesh(idx_v,2), 'ko', ...
     'MarkerSize', 8, 'LineWidth', 1.5)
plot(coord_mesh(idx_v_interni,1), coord_mesh(idx_v_interni,2), 'ro', ...
     'MarkerSize', 10, 'LineWidth', 2)
plot(coord_mesh(idx_v_esterni,1), coord_mesh(idx_v_esterni,2), 'gs', ...
     'MarkerSize', 8, 'LineWidth', 1.5)

axis equal
grid on
legend('mesh', 'tutti i vertici', 'punte interne', 'punte esterne')
title('Controllo indici vertici')

%% STEP 2: Initialize chamber, st. st.
%COMPUTE PERIMETER AND AREA
[perim0, Ap0] = eval_mesh(coordstar_cart, "cartesian");
fprintf("perimeter = %f m\n", perim0)
fprintf("area = %f m^2\n", Ap0)


% geometry
vars.geometry.port_area = Ap0;
vars.geometry.burning_perimeter = perim0;
vars.geometry.grain_length = chamber_length;
vars.geometry.throat_area = At;

% fuel
vars.fuel.a_rf = a_rf;
vars.fuel.n_rf = n_rf;
vars.fuel.rho_f = rho_f;

% combustion
vars.combustion.mdot_ox = mox;
vars.combustion.Tc_fun = T_fun_of_p;
vars.combustion.R_fun = R_fun_of_p;
vars.combustion.k_fun = k_fun_of_p;


% --- fmincon
p0 = fzero(@(pc) Z_chamber_stst(pc, vars), [pmin, pmax]);
[~, properties0] = Z_chamber_stst(p0, vars);
OF0 = properties0.O_F;
gox0 = properties0.Gox;

fprintf("Initial values:\n")
fprintf("\tpressure = %.1f bar\n", p0*1e-5);
fprintf("\tGOx = %.1f kg/m2s\n", gox0);
fprintf("\tO/F = %.2f\n", OF0);

%% STEP 3: integrate the mesh

rf = 1;         % Velocità di regressione [m/s]

Y0 = [coord_mesh(:,1); coord_mesh(:,2)];