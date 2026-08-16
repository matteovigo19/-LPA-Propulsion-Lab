function pre = run_predesign(design)
% RUN_PREDESIGN
% Esegue il predesign dell'endoreattore e restituisce tutti i risultati
% necessari alla successiva simulazione temporale.
%
% INPUT
%   design.thrust        [N]
%   design.eps           [-]
%   design.pch_bar       [bar]
%   design.O_F           [-]
%   design.GOX           [kg/(m^2 s)]
%   design.type          "cylinder", "star" oppure "RAT"
%
% Parametri necessari solo per alcune geometrie:
%   design.n_tips        numero punte, per "star"
%   design.radius_factor fattore rispetto a rad_cyl_ref, per "star"/"RAT"
%
% Parametri fuel opzionali:
%   design.rho_f         [kg/m^3]    default 920
%   design.a_rf          [(mm/s)/(kg/m^2 s)^n] default 0.027
%   design.n_rf          [-]         default 0.75
%
% OUTPUT
%   pre                  struct contenente input, geometria e risultati
%
% NOTA
%   Questa funzione NON salva prevars.mat. I risultati vengono restituiti
%   direttamente in 'pre' e possono essere passati a run_temporal_simulation.

%% ============================================================
%  INPUT
% ============================================================

thrust  = design.thrust;
eps_noz = design.eps;
pch_bar = design.pch_bar;
pch     = pch_bar * 1e5;
O_F     = design.O_F;
GOX     = design.GOX;
type    = string(design.type);

if isfield(design, "rho_f")
    rho_f = design.rho_f;
else
    rho_f = 920;
end

if isfield(design, "a_rf")
    a_rf_mm = design.a_rf;
else
    a_rf_mm = 0.027;
end

if isfield(design, "n_rf")
    n_rf = design.n_rf;
else
    n_rf = 0.75;
end

a_rf = a_rf_mm * 1e-3;   % SI: [m/s]/[(kg/m^2 s)^n]

if isfield(design, "radius_factor")
    radius_factor = design.radius_factor;
else
    radius_factor = 0.7;
end

%% ============================================================
%  CEA DATA
% ============================================================

cea = load("CEA_functions.mat");

T_fun_of_p = cea.T_fun_of_p;
k_fun_of_p = cea.k_fun_of_p;
R_fun_of_p = cea.R_fun_of_p;

%% ============================================================
%  THERMODYNAMIC / NOZZLE PREDESIGN
% ============================================================

Tch = T_fun_of_p(O_F, pch_bar);
k   = k_fun_of_p(O_F, pch_bar);
R   = R_fun_of_p(O_F, pch_bar);

% Vandenkerckhove coefficient squared
K2 = k * ((2/(k+1))^((k+1)/(k-1)));

% Characteristic velocity
cstar = sqrt(R*Tch/K2);

% Exit Mach number from nozzle expansion ratio
eps_fun = @(Me, kk) ...
    1/Me * sqrt(((1 + 0.5*(kk-1)*Me.^2)/(0.5*(kk+1)))^((kk+1)/(kk-1)));

Me = fzero(@(Me) eps_noz - eps_fun(Me, k), 2);

% Exit pressure
pe = pch / ((1 + 0.5*(k-1)*Me^2)^(k/(k-1)));

% Exit temperature
Te = Tch / (1 + 0.5*(k-1)*Me^2);

% Exit speed of sound and velocity
ae = sqrt(k*R*Te);
ue = Me*ae;

% Total mass flow
mdot = thrust / (ue + pe*eps_noz*(cstar/pch));

% Throat area and diameter
At = cstar*mdot/pch;
dt = 2*sqrt(At/pi);

% Oxidizer/fuel mass flow split
fox = mdot/(O_F + 1);
mox = O_F*fox;

% Initial port area
Ap = mox/GOX;

% Equivalent cylindrical port
diam_cyl_ref = 2*sqrt(Ap/pi);
rad_cyl_ref  = 0.5*diam_cyl_ref;

%% ============================================================
%  GEOMETRY
% ============================================================

switch lower(type)

case "star"

        n_tips = design.n_tips;
        radius = radius_factor * rad_cyl_ref;
        
        % Angolo tra una punta interna e una punta esterna consecutiva
        alpha = pi/n_tips;

        %% FUNZIONI GEOMETRICHE

        % Area del poligono formato collegando le punte interne.
        % ri è il raggio circoscritto, non l'apotema.
        Area_inner_polygon = @(ri) ...
            0.5*n_tips*ri.^2*sin(2*alpha);

        % Lato del poligono che collega due punte interne consecutive
        side_inner_polygon = @(ri) ...
            2*ri.*sin(alpha);

        % Area di uno dei triangoli esterni.
        %
        % La base è il lato del poligono interno.
        % L'altezza è:
        % re - ri*cos(alpha)
        Area_external_triangle = @(re,ri) ...
            0.5.*side_inner_polygon(ri).* ...
            (re - ri.*cos(alpha));

        % Area complessiva dei triangoli esterni
        Area_external_triangles = @(re,ri) ...
            n_tips.*Area_external_triangle(re,ri);

        % Area complessiva della stella
        Ap_fun = @(re,ri) ...
            Area_inner_polygon(ri) + ...
            Area_external_triangles(re,ri);

        % Formula equivalente semplificata:
        Ap_fun_direct = @(re,ri) ...
            n_tips.*re.*ri.*sin(alpha);

        % Lunghezza del segmento tra una punta interna e una esterna
        star_side_fun = @(re,ri) ...
            sqrt(re.^2 + ri.^2 - ...
            2.*re.*ri.*cos(alpha));

        % Perimetro complessivo della stella
        Pb_fun = @(re,ri) ...
            2*n_tips.*star_side_fun(re,ri);


        %% DETERMINAZIONE DI ri E re

        if radius < rad_cyl_ref

            % radius rappresenta il raggio delle punte interne
            ri = radius;

            % Dall'equazione:
            % Ap = n_tips*ri*re*sin(alpha)
            re = Ap/(n_tips*ri*sin(alpha));

        elseif radius > rad_cyl_ref

            % radius rappresenta il raggio delle punte esterne
            re = radius;

            % Dall'equazione:
            % Ap = n_tips*ri*re*sin(alpha)
            ri = Ap/(n_tips*re*sin(alpha));

        else

            % Caso degenere: tutte le punte sono sulla stessa
            % circonferenza.
            ri = radius;
            re = radius;

        end


        %% CONTROLLI DI COERENZA

        if ri <= 0 || re <= 0
            error('I raggi ri e re devono essere positivi.');
        end

        if re < ri
            error(['Geometria non valida: il raggio delle punte esterne ' ...
                   're risulta minore del raggio interno ri.']);
        end


        %% AREA E PERIMETRO RISULTANTI

        Ap_geometry = Ap_fun(re,ri);
        Pb = Pb_fun(re,ri);

        % Controllo numerico dell'area
        area_relative_error = abs(Ap_geometry - Ap)/max(abs(Ap),eps);

        if area_relative_error > 1e-10
            warning(['L''area geometrica non coincide con Ap. ' ...
                     'Errore relativo: %.3e'], ...
                     area_relative_error);
        end


        %% LUNGHEZZA DEL GRANO

        L = fox/(rho_f*Pb*a_rf*(mox^n_rf)/Ap^n_rf);


        %% DIAMETRI CARATTERISTICI

        diam_i = 2*ri;
        diam_e = 2*re;

        geom.n_tips = n_tips;
        geom.ri = ri;
        geom.re = re;
        geom.diam_i = diam_i;
        geom.diam_e = diam_e;
        geom.radius = radius;

    case "rat"

        if isfield(design, "radius")
            radius = design.radius;
        else
            radius = rad_cyl_ref * radius_factor;
        end

        area_in = @(ri) pi*ri.^2;
        area_e  = @(re) pi*re.^2;
        Ap_fun  = @(re,ri) area_e(re) - area_in(ri);

        if radius <= rad_cyl_ref

            ri = radius;
            re = fzero(@(re) Ap - Ap_fun(re,ri), rad_cyl_ref);

        else

            re = radius;
            ri = fzero(@(ri) Ap - Ap_fun(re,ri), rad_cyl_ref);

        end

        Pb = 2*pi*(re + ri);

        L = fox / ...
            (rho_f*Pb*a_rf*(mox^n_rf)/Ap^n_rf);

        diam_in = 2*ri;
        diam_e  = 2*re;

        geom.ri = ri;
        geom.re = re;
        geom.diam_in = diam_in;
        geom.diam_e = diam_e;
        geom.radius = radius;

    case "cylinder"

        diam = diam_cyl_ref;
        Pb = pi*diam;

        L = fox / ...
            (rho_f*Pb*a_rf*(mox^n_rf)/Ap^n_rf);

        geom.diam = diam;
        geom.ri = diam/2;
        geom.re = diam/2;

    otherwise

        error("Geometry case not defined: %s", type);

end

%% ============================================================
%  OUTPUT STRUCT
% ============================================================

% Input/design variables
pre.design = design;

pre.thrust = thrust;
pre.eps = eps_noz;
pre.pch_bar = pch_bar;
pre.pch = pch;
pre.O_F = O_F;
pre.GOX = GOX;
pre.type = type;

% Fuel
pre.rho_f = rho_f;
pre.a_rf = a_rf;
pre.a_rf_mm = a_rf_mm;
pre.n_rf = n_rf;

% Thermodynamics / nozzle
pre.Tch = Tch;
pre.k = k;
pre.R = R;
pre.K2 = K2;
pre.cstar = cstar;
pre.Me = Me;
pre.pe = pe;
pre.Te = Te;
pre.ae = ae;
pre.ue = ue;

pre.mdot = mdot;
pre.fox = fox;
pre.mox = mox;

pre.At = At;
pre.dt = dt;

% Port/grain geometry
pre.Ap = Ap;
pre.Pb = Pb;
pre.L = L;

pre.diam_cyl_ref = diam_cyl_ref;
pre.rad_cyl_ref = rad_cyl_ref;

pre.geometry = geom;

% Replicate convenient flat fields used by the old temporal main
switch lower(type)
    case "star"
        pre.n_tips = geom.n_tips;
        pre.diam_i = geom.diam_i;
        pre.diam_e = geom.diam_e;

    case "rat"
        pre.diam_in = geom.diam_in;
        pre.diam_e = geom.diam_e;

    case "cylinder"
        pre.diam = geom.diam;
end

% CEA handles can be reused directly by the temporal simulation
pre.CEA.T_fun_of_p = T_fun_of_p;
pre.CEA.k_fun_of_p = k_fun_of_p;
pre.CEA.R_fun_of_p = R_fun_of_p;

end
