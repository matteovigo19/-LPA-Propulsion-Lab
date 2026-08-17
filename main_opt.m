close all
clear
clc

%% ============================================================
%  MAIN OPT
%  Predesign + simulazione temporale + opt
% ============================================================


%% ============================================================
%  1. DATI DI PREDESIGN
% ============================================================

design.thrust = 50000;       % [N]

% Rapporto di espansione ugello
design.eps = 200;            % [-]

% Pressione iniziale desiderata in camera
design.pch_bar = 20;         % [bar]

% Rapporto ossidante / combustibile iniziale
design.O_F = 2;              % [-]

% Flusso massico di ossidante iniziale
design.GOX = 500;            % [kg/(m^2 s)]


%% ============================================================
%  2. COMBUSTIBILE
% ============================================================

% HTPB

design.rho_f = 920;          % [kg/m^3]

% Legge di regressione:
%
%   rf = a * GOX^n
%

design.a_rf = 0.027;         % [(mm/s)/(kg/(m^2 s))^n]
design.n_rf = 0.75;          % [-]


%% ============================================================
%  3. GEOMETRIA
% ============================================================

% Scegliere una delle geometrie disponibili:
%
%   "cylinder"
%   "star"
%   "RAT"
%

design.type = "star";


switch lower(design.type)

    case "star"

        % Numero di punte della stella
        design.n_tips = 5;

        % Raggio scelto rispetto al raggio del cilindro equivalente:
        %
        % radius = radius_factor * rad_cyl_ref
        %
        % Se radius_factor < 1:
        % il raggio assegnato viene interpretato come raggio interno.

        design.radius_factor = 0.7;


    case "rat"

        design.radius_factor = 0.7;


    case "cylinder"

        % Nessun parametro geometrico aggiuntivo necessario


    otherwise

        error("Geometria non riconosciuta");

end


%% ============================================================
%  4. ESECUZIONE PREDESIGN
% ============================================================

pre = run_predesign(design);


%% ============================================================
%  5. STAMPA RISULTATI PREDESIGN
% ============================================================

fprintf("\n");
fprintf("=============================================\n");
fprintf("              RISULTATI PREDESIGN\n");
fprintf("=============================================\n");

fprintf("Geometria:              %s\n", pre.type);

fprintf("\nPARAMETRI OPERATIVI\n");

fprintf("Spinta iniziale:        %.3f N\n", pre.thrust);
fprintf("Pressione iniziale:     %.3f bar\n", pre.pch_bar);
fprintf("O/F iniziale:           %.4f\n", pre.O_F);
fprintf("GOX iniziale:           %.3f kg/(m^2 s)\n", pre.GOX);

fprintf("\nPORTATE\n");

fprintf("Portata totale:         %.6f kg/s\n", pre.mdot);
fprintf("Portata ossidante:      %.6f kg/s\n", pre.mox);
fprintf("Portata combustibile:   %.6f kg/s\n", pre.fox);

fprintf("\nGEOMETRIA PORTO\n");

fprintf("Area porto Ap:          %.8e m^2\n", pre.Ap);
fprintf("Perimetro Pb:           %.8e m\n", pre.Pb);
fprintf("Lunghezza grano L:      %.8e m\n", pre.L);

fprintf("\nUGELLO\n");

fprintf("Area di gola:           %.8e m^2\n", pre.At);
fprintf("Diametro di gola:       %.6f m\n", pre.dt);
fprintf("Mach in uscita:         %.4f\n", pre.Me);
fprintf("Pressione in uscita:    %.3f Pa\n", pre.pe);

fprintf("\n");


%% ============================================================
%  6. IMPOSTAZIONI SIMULAZIONE TEMPORALE
% ============================================================

settings = struct();


%% ------------------------------------------------------------
% Tempo massimo simulazione
% -------------------------------------------------------------

settings.tmax = 300;         % [s]


%% ------------------------------------------------------------
% Dimensione camera di combustione
% -------------------------------------------------------------

settings.ext_diameter = 280.1*2e-3;    % [m]


%% ------------------------------------------------------------
% Discretizzazione geometrica
% -------------------------------------------------------------

% Per la stella:
%
% n = multiplyer * n_tips

settings.multiplyer = 80;


% Per il cilindro

settings.n_cylinder = 500;


%% ------------------------------------------------------------
% ODE
% -------------------------------------------------------------

settings.RelTol = 1e-8;
settings.AbsTol = 1e-10;

settings.fine_ode_boolean = false;


%% ------------------------------------------------------------
% Tipo di mesh utilizzata
% -------------------------------------------------------------

% "ode45"
%       usa la mesh grezza integrata dall'ODE
%
% "refined"
%       applica refine_mesh_v3 nel post-processing

settings.plot_case = "refined";


%% ------------------------------------------------------------
% Numero di output temporali richiesti
% -------------------------------------------------------------

settings.n_output = 100;


%% ------------------------------------------------------------
% Output grafici
% -------------------------------------------------------------

settings.make_plots = true;

settings.make_animation = true;

settings.verbose = true;


%% ============================================================
%  7. ESECUZIONE SIMULAZIONE TEMPORALE
% ============================================================

sim = run_temporal_simulation(pre, settings);


%% ============================================================
%  8. RISULTATI FINALI
% ============================================================

fprintf("\n");
fprintf("=============================================\n");
fprintf("          RISULTATI SIMULAZIONE\n");
fprintf("=============================================\n");

fprintf("Tempo finale:            %.6f s\n", ...
    sim.t(end));

fprintf("Pressione iniziale:      %.6f bar\n", ...
    sim.p(1)*1e-5);

fprintf("Pressione finale:        %.6f bar\n", ...
    sim.p(end)*1e-5);

fprintf("\n");

fprintf("GOX iniziale:            %.6f kg/(m^2 s)\n", ...
    sim.Gox(1));

fprintf("GOX finale:              %.6f kg/(m^2 s)\n", ...
    sim.Gox(end));

fprintf("\n");

fprintf("O/F iniziale:            %.6f\n", ...
    sim.OF(1));

fprintf("O/F finale:              %.6f\n", ...
    sim.OF(end));

fprintf("\n");

fprintf("Area iniziale:           %.8e m^2\n", ...
    sim.area_m2(1));

fprintf("Area finale:             %.8e m^2\n", ...
    sim.area_m2(end));

fprintf("\n");

fprintf("Perimetro iniziale:      %.8e m\n", ...
    sim.perimeter_m(1));

fprintf("Perimetro finale:        %.8e m\n", ...
    sim.perimeter_m(end));


%% ============================================================
%  9. EVENTUALE EVENTO DI FINE COMBUSTIONE
% ============================================================

if sim.event_occurred

    fprintf("\n");
    fprintf("Camera raggiunta dalla regressione.\n");
    fprintf("Tempo evento: %.6f s\n", sim.t_event(end));

else

    fprintf("\n");
    fprintf("Nessun evento di raggiungimento camera.\n");

end