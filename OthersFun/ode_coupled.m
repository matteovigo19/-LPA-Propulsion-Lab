function dy = ode_coupled(~, y, vars, fine_ode_boolean, plot_case, type)

% DESCRIPTION:
%   ODE accoppiata pressione-camera + regressione mesh.
%
% Stato:
%   y(1)     = pressione camera pc [Pa]
%   y(2:end) = mesh [x1...xn, y1...yn]
%
% Quando la geometria è ancora stellata:
%   - uso refine_mesh_v3
%   - uso ode_mesh_v4 con correzione delle punte interne
%
% Quando la geometria è quasi cilindrica:
%   - non uso refine_mesh_v3
%   - uso ode_mesh_v4 con sola regressione normale

%% ============================================================
%  1. UNPACK STATO
% ============================================================

pc = y(1);          % [Pa]
Y  = y(2:end);
Y  = Y(:);

%% ============================================================
%  2. UNPACK VARS
% ============================================================

% geometry
% geometry
L  = vars.geometry.grain_length;
At = vars.geometry.throat_area;

D_camera = vars.geometry.D_camera;
R_camera = D_camera / 2;

idx_v         = vars.geometry.idx_v;
idx_v_esterni = vars.geometry.idx_v_esterni;
idx_v_interni = vars.geometry.idx_v_interni;

% fuel
a_rf  = vars.fuel.a_rf;
n_rf  = vars.fuel.n_rf;
rho_f = vars.fuel.rho_f;

% combustion
mdot_ox = vars.combustion.mdot_ox;
Tc_fun  = vars.combustion.Tc_fun;
R_fun   = vars.combustion.R_fun;
k_fun   = vars.combustion.k_fun;

%% ============================================================
%  3. VALUTAZIONE MESH GREZZA
% ============================================================

P_raw = stato_to_punti_locale(Y);

% ============================================================
%  CORREZIONE CAMERA SEMPRE ATTIVA
% ============================================================
% Se qualche punto supera la camera, viene proiettato sulla parete.

[P_raw, ~] = refine_mesh_camera(P_raw, D_camera);
Y = punti_to_stato_locale(P_raw);

[~, Ap_raw] = eval_mesh(P_raw, "cartesian");

if Ap_raw <= 0
    error("Area portante non positiva in ode_coupled.");
end

Gox_raw = mdot_ox / Ap_raw;
rf = a_rf * Gox_raw^n_rf;

%% ============================================================
%  4. CONTROLLO GEOMETRIA QUASI CILINDRICA
% ============================================================
% Se la differenza media tra raggi delle punte esterne e interne
% diventa piccola, la geometria non è più realmente una stella.
%
% In quel caso le vecchie punte interne non devono più essere trattate
% come intersezioni di due lati regrediti.

is_cylindrical = is_quasi_cylindrical_geometry( ...
    P_raw, ...
    idx_v_esterni, ...
    idx_v_interni, ...
    vars);

if is_cylindrical
    mesh_mode = "cylindrical";
else
    mesh_mode = "star";
end

%% ============================================================
%  5. EVENTUALE RAFFINAMENTO MESH CORRENTE
% ============================================================

switch plot_case

    case "refined"

        P_eval = refine_mesh_v3(P_raw, idx_v_interni);

        % Correzione camera sempre attiva anche dopo refine_mesh_v3
        [P_eval, ~] = refine_mesh_camera(P_eval, D_camera);

        Y_eval = punti_to_stato_locale(P_eval);

    case "ode45"

    P_eval = P_raw;

    % Anche se non raffino la stella, rispetto comunque la camera
    [P_eval, ~] = refine_mesh_camera(P_eval, D_camera);

    Y_eval = punti_to_stato_locale(P_eval);

    otherwise

        error(['plot_case non valido in ode_coupled. ', ...
               'Usa "refined" oppure "ode45".']);

end

%% ============================================================
%  6. AREA, PERIMETRO E REGRESSIONE SULLA MESH SCELTA
% ============================================================

[perimB, Ap] = eval_mesh(P_eval, "cartesian");

if Ap <= 0
    error("Area portante valutata non positiva in ode_coupled.");
end

Gox = mdot_ox / Ap;
rf = a_rf * Gox^n_rf;

Ab = perimB * L;

mdot_f  = rho_f * Ab * rf;
mdot_in = mdot_f + mdot_ox;

O_F = mdot_ox / mdot_f;

%% ============================================================
%  7. USCITA UGELLO
% ============================================================

pc_bar = pc*1e-5;

%% ============================================================
%  CHECK DOMINIO CEA
% ============================================================

if isfield(vars.combustion, "OFmin")
    OFmin = vars.combustion.OFmin;
else
    OFmin = -inf;
end

if isfield(vars.combustion, "OFmax")
    OFmax = vars.combustion.OFmax;
else
    OFmax = inf;
end

if isfield(vars.combustion, "pcmin_bar")
    pcmin_bar = vars.combustion.pcmin_bar;
else
    pcmin_bar = -inf;
end

if isfield(vars.combustion, "pcmax_bar")
    pcmax_bar = vars.combustion.pcmax_bar;
else
    pcmax_bar = inf;
end

if ~isfinite(O_F) || O_F <= 0
    error("O_F non valido in ode_coupled: O_F = %.6e", O_F);
end

if ~isfinite(pc_bar) || pc_bar <= 0
    error("pc_bar non valida in ode_coupled: pc_bar = %.6e", pc_bar);
end

if O_F < OFmin || O_F > OFmax
    error("O_F fuori dal dominio CEA: O_F = %.6f, dominio = [%.6f, %.6f]", ...
        O_F, OFmin, OFmax);
end

if pc_bar < pcmin_bar || pc_bar > pcmax_bar
    error("pc_bar fuori dal dominio CEA: pc_bar = %.6f bar, dominio = [%.6f, %.6f] bar", ...
        pc_bar, pcmin_bar, pcmax_bar);
end

%% ============================================================
%  VALUTAZIONE PROPRIETÀ CEA
% ============================================================

Tc = Tc_fun(O_F, pc_bar);
R  = R_fun(O_F, pc_bar);
k  = k_fun(O_F, pc_bar);

if any(~isfinite([Tc, R, k]))
    error("CEA ha restituito NaN/Inf: O_F = %.6f, pc_bar = %.6f, Tc = %.6e, R = %.6e, k = %.6e", ...
        O_F, pc_bar, Tc, R, k);
end

if Tc <= 0 || R <= 0 || k <= 1
    error("Proprietà CEA non fisiche: Tc = %.6e, R = %.6e, k = %.6e", ...
        Tc, R, k);
end

K2 = k * ((2/(k+1))^((k+1)/(k-1)));

if ~isfinite(K2) || K2 <= 0
    error("K2 non valido: K2 = %.6e, k = %.6e", K2, k);
end

cstar = sqrt(R * Tc / K2);

if ~isfinite(cstar) || cstar <= 0
    error("cstar non valida: cstar = %.6e", cstar);
end

mdot_out = pc * At / cstar;

rho_g = pc / (R * Tc);

if any(~isfinite([mdot_out, rho_g]))
    error("mdot_out o rho_g NaN: mdot_out = %.6e, rho_g = %.6e", ...
        mdot_out, rho_g);
end

%% ============================================================
%  8. DINAMICA MESH
% ============================================================
% Se mesh_mode == "star":
%   ode_mesh_v4 corregge le punte interne.
%
% Se mesh_mode == "cylindrical":
%   ode_mesh_v4 usa solo regressione normale per tutti i punti.

dY = ode_mesh_v4(0, Y, idx_v_interni, rf, mesh_mode);

% I punti sulla camera restano fissi
dY = freeze_points_on_camera(Y, dY, D_camera);

%% ============================================================
%  9. DERIVATE AREA/PERIMETRO PER EQUAZIONE DI PRESSIONE
% ============================================================

dY_eval = ode_mesh_v4(0, Y_eval, idx_v_interni, rf, mesh_mode);

% Coerenza anche per il calcolo di dA/dt e dP/dt
dY_eval = freeze_points_on_camera(Y_eval, dY_eval, D_camera);

[d_perim_perim, d_area_area, dperimeter, darea] = ...
    ode_eval_mesh(Y_eval, dY_eval);

%% ============================================================
%  10. EQUAZIONE PRESSIONE CAMERA
% ============================================================

vars.geometry.grain_length = L;
vars.geometry.port_area = Ap;
vars.geometry.d_area_area = d_area_area;

if fine_ode_boolean

    dRTdOF_fun_OF_p = vars.combustion.dRTdOF_fun_OF_p;
    dRTdp_fun_OF_p  = vars.combustion.dRTdp_fun_OF_p;

    dRTdOF = dRTdOF_fun_OF_p(O_F, pc_bar);
    dRTdp  = dRTdp_fun_OF_p(O_F, pc_bar);

    deltapc = pc_bar * dRTdp  / (R * Tc);
    deltaOF = O_F    * dRTdOF / (R * Tc);

    fine_ode.deltapc = deltapc;
    fine_ode.deltaOF = deltaOF;
    fine_ode.n_rf = n_rf;
    fine_ode.d_perim_perim = d_perim_perim;

    dp_p = ode_chamber_pressure( ...
        mdot_in, ...
        mdot_out, ...
        rho_g, ...
        vars, ...
        fine_ode);

else

    dp_p = ode_chamber_pressure( ...
        mdot_in, ...
        mdot_out, ...
        rho_g, ...
        vars);

end

dp = dp_p * pc;

%% ============================================================
%  11. ASSEMBLAGGIO OUTPUT
% ============================================================

dy = [dp; dY(:)];

end


%% ============================================================
%  FUNZIONE LOCALE: CONTROLLO GEOMETRIA QUASI CILINDRICA
% ============================================================

function is_cyl = is_quasi_cylindrical_geometry(P, idx_v_esterni, idx_v_interni, vars)

idx_v_esterni = idx_v_esterni(:);
idx_v_interni = idx_v_interni(:);

r = vecnorm(P, 2, 2);

r_ext = mean(r(idx_v_esterni));
r_int = mean(r(idx_v_interni));

r_mean = mean(r);

% Profondità media delle gole della stella
star_depth = r_ext - r_int;

% Soglia relativa.
% Quando la differenza tra raggi esterni e interni è inferiore,
% per esempio, al 3% del raggio medio, considero la geometria quasi cilindrica.
if isfield(vars.geometry, "cyl_tol_rel")
    cyl_tol_rel = vars.geometry.cyl_tol_rel;
else
    cyl_tol_rel = 0.03;
end

% Soglia assoluta legata alla risoluzione media della mesh.
% Serve per evitare switch instabili dovuti a rumore numerico.
ds = vecnorm(diff(P,1,1), 2, 2);
h_ref = median(ds);

if isfield(vars.geometry, "cyl_tol_abs")
    cyl_tol_abs = vars.geometry.cyl_tol_abs;
else
    cyl_tol_abs = 3*h_ref;
end

depth_tol = max(cyl_tol_abs, cyl_tol_rel*r_mean);

is_cyl = star_depth < depth_tol;

end


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