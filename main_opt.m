close all
clear
clc

%% ============================================================
%  MAIN OPT
%
%  Per ora viene considerata solamente la geometria STAR.
%
%  Variabile discreta:
%       - numero di punte
%
%  Variabili ottimizzate da fmincon:
%       x(1) = O/F iniziale
%       x(2) = GOX iniziale
%       x(3) = thrust iniziale
%       x(4) = radius factor
%
%  Funzione obiettivo:
%       massimizzazione impulso totale
%
%  Vincolo:
%       thrust medio >= thrust medio minimo
% ============================================================

tic;

%% ============================================================
%  1. DATI GENERALI DEL DESIGN
% ============================================================

% Rapporto di espansione ugello
design.eps = 200;                 % [-]

% Pressione iniziale desiderata in camera
design.pch_bar = 20;              % [bar]


%% ============================================================
%  2. COMBUSTIBILE
% ============================================================

% HTPB

design.rho_f = 920;               % [kg/m^3]

% Legge di regressione:
%
%       rf = a * GOX^n

design.a_rf = 0.027;              % [(mm/s)/(kg/(m^2 s))^n]

design.n_rf = 0.75;               % [-]


%% ============================================================
%  3. GEOMETRIA
% ============================================================

% Per il momento viene ottimizzata solamente la geometria STAR

design.type = "star";


%% ============================================================
%  4. LIMITI DELLE VARIABILI DI OTTIMIZZAZIONE
% ============================================================

%% ------------------------------------------------------------
% O/F iniziale
% -------------------------------------------------------------

OF_min = 1.2;
OF_max_global = 4.0;

% Limite superiore dell'O/F iniziale in funzione del numero di punte.
%
% Modifica qui i valori quando avrai determinato limiti più accurati.
%
% Esempio attuale:
%   10 punte -> O/F max = 2.9
%   12 punte -> O/F max = 2.6
%
% Per i numeri di punte non presenti nella mappa viene usato
% OF_max_global.

OF_max_by_tips = containers.Map( ...
    [10, 12, 14], ...
    [2.7, 2.4, 2.1]);


%% ------------------------------------------------------------
% GOX iniziale
% -------------------------------------------------------------

GOX_min = 300;                    % [kg/(m^2 s)]
GOX_max = 800;                    % [kg/(m^2 s)]


%% ------------------------------------------------------------
% Thrust iniziale
% -------------------------------------------------------------

thrust_min = 45000;               % [N]
thrust_max = 55000;               % [N]


%% ------------------------------------------------------------
% Radius factor
%
% radius = radius_factor * rad_cyl_ref
% -------------------------------------------------------------

radius_factor_min = 0.50;
radius_factor_max = 0.95;


%% ------------------------------------------------------------
% Numero di punte
%
% Non viene passato direttamente a fmincon perché è una
% variabile discreta.
% -------------------------------------------------------------

n_tips_min = 10;
n_tips_max = 12;

n_tips_values = n_tips_min:2:n_tips_max;


%% ============================================================
%  5. BOUNDS DELLE VARIABILI DI OTTIMIZZAZIONE
%
%  I bounds vengono costruiti localmente nel ciclo su n_tips,
%  perché il limite superiore dell'O/F iniziale dipende
%  dal numero di punte.
% ============================================================


%% ============================================================
%  6. IMPOSTAZIONI SIMULAZIONE TEMPORALE
% ============================================================

settings = struct();


%% ------------------------------------------------------------
% Tempo massimo simulazione
% -------------------------------------------------------------

settings.tmax = 300;               % [s]


%% ------------------------------------------------------------
% Camera
% -------------------------------------------------------------

settings.ext_diameter = 1000.1*2e-3;     % [m]


%% ------------------------------------------------------------
% Discretizzazione mesh
%
% Per la stella:
%
%       n = multiplyer * n_tips
% -------------------------------------------------------------

settings.multiplyer = 10;


%% ------------------------------------------------------------
% ODE
% -------------------------------------------------------------

settings.RelTol = 1e-8;

settings.AbsTol = 1e-10;

settings.fine_ode_boolean = false;


%% ------------------------------------------------------------
% Tipo di trattamento mesh
% -------------------------------------------------------------

settings.plot_case = "refined";


%% ------------------------------------------------------------
% Numero di output temporali
% -------------------------------------------------------------

settings.n_output = 10;


%% ------------------------------------------------------------
% Durante l'ottimizzazione:
%
% NON voglio figure, animazioni o stampe del temporale
% -------------------------------------------------------------

settings.make_plots = false;

settings.make_animation = false;

settings.verbose = false;


%% ============================================================
%  VINCOLI SULLA SIMULAZIONE
% ============================================================

% Spinta media minima
thrust_mean_min = 50000;       % [N]

% O/F ammesso durante TUTTA la simulazione
OF_sim_min = 1.15;             % [-]
OF_sim_max = 19.45;            % [-]

% Pressione di camera ammessa durante TUTTA la simulazione
pc_sim_min_bar = 101325 * 1e-5;           % [bar]
pc_sim_max_bar = 99;           % [bar]

% Conversione in Pa
pc_sim_min = pc_sim_min_bar * 1e5;
pc_sim_max = pc_sim_max_bar * 1e5;

%% ============================================================
%  8. OPZIONI FMINCON
% ============================================================

options = optimoptions( ...
    "fmincon", ...
    "Algorithm", "sqp", ...
    "Display", "iter", ...
    "MaxIterations", 100, ...
    "MaxFunctionEvaluations", 1000, ...
    "ConstraintTolerance", 1e-4, ...
    "StepTolerance", 1e-8);


%% ============================================================
%  9. INIZIALIZZAZIONE OTTIMIZZAZIONE
% ============================================================

rng("shuffle");

n_configurations = length(n_tips_values);

% Numero di guess casuali per ogni numero di punte
n_guess = 2;

results = struct([]);


%% ============================================================
%  10. CICLO SUL NUMERO DI PUNTE
% ============================================================

for ii = 1:n_configurations

    n_tips = n_tips_values(ii);

    %% ========================================================
    %  BOUND O/F LOCALE DIPENDENTE DAL NUMERO DI PUNTE
    % ========================================================

    if isKey(OF_max_by_tips, n_tips)

        OF_max_local = OF_max_by_tips(n_tips);

    else

        OF_max_local = OF_max_global;

    end

    %% ========================================================
    %  BOUNDS LOCALI FMINCON
    % ========================================================

    lb_local = [
        OF_min
        GOX_min
        thrust_min
        radius_factor_min
    ];

    ub_local = [
        OF_max_local
        GOX_max
        thrust_max
        radius_factor_max
    ];


    fprintf("\n");
    fprintf("=====================================================\n");
    fprintf("       OTTIMIZZAZIONE STAR - %d PUNTE\n", n_tips);
    fprintf("=====================================================\n");


    %% ========================================================
    %  INIZIALIZZAZIONE MIGLIOR RISULTATO PER QUESTO N_TIPS
    % ========================================================

    best_fval_local = Inf;

    best_x_opt_local = [];
    best_x0_local = [];

    best_exitflag_local = [];
    best_output_local = [];

    best_pre_local = [];
    best_sim_local = [];


    %% ========================================================
    %  CICLO SUI GUESS CASUALI
    % ========================================================

    for jj = 1:n_guess

        fprintf("\n");
        fprintf("---------------------------------------------\n");
        fprintf("Guess %d / %d - %d punte\n", ...
            jj, n_guess, n_tips);
        fprintf("---------------------------------------------\n");


        %% ====================================================
        %  10.1 GENERAZIONE PUNTO INIZIALE CASUALE
        % =====================================================

        x0 = lb_local + rand(size(lb_local)).*(ub_local-lb_local);

        fprintf("Range O/F       = [%.4f, %.4f]\n", OF_min, OF_max_local);
        fprintf("O/F             = %.4f\n", x0(1));
        fprintf("GOX             = %.2f kg/(m^2 s)\n", x0(2));
        fprintf("Thrust          = %.2f N\n", x0(3));
        fprintf("Radius factor   = %.4f\n", x0(4));


        %% ====================================================
        %  10.2 FUNZIONE OBIETTIVO
        % =====================================================

        objective = @(x) objective_function( ...
            x, ...
            n_tips, ...
            design, ...
            settings);


        %% ====================================================
        %  10.3 VINCOLI NON LINEARI
        % =====================================================

        nonlcon = @(x) optimization_constraints( ...
            x, ...
            n_tips, ...
            design, ...
            settings, ...
            thrust_mean_min, ...
            OF_sim_min, ...
            OF_sim_max, ...
            pc_sim_min, ...
            pc_sim_max);

        %% ====================================================
        %  10.4 FMINCON
        % =====================================================

        [x_opt, fval, exitflag, output] = fmincon( ...
            objective, ...
            x0, ...
            [], ...
            [], ...
            [], ...
            [], ...
            lb_local, ...
            ub_local, ...
            nonlcon, ...
            options);


        %% ====================================================
        %  10.5 SIMULAZIONE DEL RISULTATO
        % =====================================================

        design_opt = design;

        design_opt.O_F = x_opt(1);
        design_opt.GOX = x_opt(2);
        design_opt.thrust = x_opt(3);
        design_opt.radius_factor = x_opt(4);

        design_opt.n_tips = n_tips;


        pre_opt = run_predesign(design_opt);

        sim_opt = run_temporal_simulation( ...
            pre_opt, ...
            settings);


        %% ====================================================
        %  10.6 STAMPA RISULTATO DEL GUESS
        % =====================================================

        fprintf("\nRisultato guess %d:\n", jj);

        fprintf("O/F             = %.6f\n", ...
            x_opt(1));

        fprintf("GOX             = %.6f kg/(m^2 s)\n", ...
            x_opt(2));

        fprintf("Thrust iniziale = %.6f N\n", ...
            x_opt(3));

        fprintf("Radius factor   = %.6f\n", ...
            x_opt(4));

        fprintf("Thrust medio    = %.6f N\n", ...
            sim_opt.thrust_mean);

        fprintf("Impulso totale  = %.6e N s\n", ...
            sim_opt.total_impulse);

        fprintf("fval            = %.6e\n", ...
            fval);


        %% ====================================================
        %  10.7 CONTROLLO SE È IL MIGLIORE
        %
        %  Siccome:
        %
        %       J = -I_tot
        %
        %  il risultato migliore è quello con fval MINORE.
        % =====================================================

        if fval < best_fval_local

            best_fval_local = fval;

            best_x_opt_local = x_opt;

            best_x0_local = x0;

            best_exitflag_local = exitflag;

            best_output_local = output;

            best_pre_local = pre_opt;

            best_sim_local = sim_opt;

            fprintf("\n>>> Nuovo miglior risultato per %d punte <<<\n", ...
                n_tips);

        end

    end


    %% ========================================================
    %  10.8 SALVATAGGIO MIGLIOR RISULTATO PER QUESTO N_TIPS
    % ========================================================

    results(ii).n_tips = n_tips;

    results(ii).n_guess = n_guess;


    %% Punto iniziale che ha portato alla soluzione migliore

    results(ii).x0 = best_x0_local;


    %% Parametri ottimizzati

    results(ii).x_opt = best_x_opt_local;

    results(ii).OF_opt = best_x_opt_local(1);

    results(ii).GOX_opt = best_x_opt_local(2);

    results(ii).thrust0_opt = best_x_opt_local(3);

    results(ii).radius_factor_opt = ...
        best_x_opt_local(4);


    %% Prestazioni

    results(ii).thrust_mean = ...
        best_sim_local.thrust_mean;

    results(ii).total_impulse = ...
        best_sim_local.total_impulse;


    %% Risultati fmincon

    results(ii).fval = ...
        best_fval_local;

    results(ii).exitflag = ...
        best_exitflag_local;

    results(ii).output = ...
        best_output_local;


    %% Simulazione completa

    results(ii).pre = ...
        best_pre_local;

    results(ii).sim = ...
        best_sim_local;


    %% ========================================================
    %  10.9 STAMPA MIGLIOR RISULTATO PER QUESTO NUMERO DI PUNTE
    % ========================================================

    fprintf("\n");
    fprintf("=====================================================\n");
    fprintf(" MIGLIOR RISULTATO %d PUNTE SU %d GUESS\n", ...
        n_tips, n_guess);
    fprintf("=====================================================\n");

    fprintf("O/F                   = %.6f\n", ...
        results(ii).OF_opt);

    fprintf("GOX                   = %.6f kg/(m^2 s)\n", ...
        results(ii).GOX_opt);

    fprintf("Thrust iniziale       = %.6f N\n", ...
        results(ii).thrust0_opt);

    fprintf("Radius factor         = %.6f\n", ...
        results(ii).radius_factor_opt);

    fprintf("\n");

    fprintf("Thrust medio          = %.6f N\n", ...
        results(ii).thrust_mean);

    fprintf("Thrust medio minimo   = %.6f N\n", ...
        thrust_mean_min);

    fprintf("Errore thrust medio   = %.6f N\n", ...
        results(ii).thrust_mean - thrust_mean_min);

    fprintf("\n");

    fprintf("Impulso totale        = %.6e N s\n", ...
        results(ii).total_impulse);

end


%% ============================================================
%  11. CONFRONTO DELLE CONFIGURAZIONI
% ============================================================

fprintf("\n");
fprintf("=====================================================\n");
fprintf("             CONFRONTO CONFIGURAZIONI\n");
fprintf("=====================================================\n");


for ii = 1:n_configurations

    fprintf("\n%d punte\n", ...
        results(ii).n_tips);

    fprintf("O/F             = %.6f\n", ...
        results(ii).OF_opt);

    fprintf("GOX             = %.6f kg/(m^2 s)\n", ...
        results(ii).GOX_opt);

    fprintf("Thrust iniziale = %.6f N\n", ...
        results(ii).thrust0_opt);

    fprintf("Radius factor   = %.6f\n", ...
        results(ii).radius_factor_opt);

    fprintf("Thrust medio    = %.6f N\n", ...
        results(ii).thrust_mean);

    fprintf("Impulso totale  = %.6e N s\n", ...
        results(ii).total_impulse);

end


%% ============================================================
%  12. SCELTA DELLA CONFIGURAZIONE MIGLIORE
%
%  La configurazione migliore è quella con impulso totale
%  massimo.
% ============================================================

impulses = [results.total_impulse];


[best_impulse, idx_best] = max(impulses);


best = results(idx_best);


%% ============================================================
%  13. RISULTATI OTTIMIZZAZIONE
% ============================================================

fprintf("\n");
fprintf("=====================================================\n");
fprintf("             CONFIGURAZIONE MIGLIORE\n");
fprintf("=====================================================\n");

fprintf("\n");

fprintf("Geometria             = %s\n", ...
    design.type);

fprintf("Numero punte          = %d\n", ...
    best.n_tips);

fprintf("\n");

fprintf("O/F iniziale          = %.6f\n", ...
    best.OF_opt);

fprintf("GOX iniziale          = %.6f kg/(m^2 s)\n", ...
    best.GOX_opt);

fprintf("Thrust iniziale       = %.6f N\n", ...
    best.thrust0_opt);

fprintf("Radius factor         = %.6f\n", ...
    best.radius_factor_opt);

fprintf("\n");

fprintf("Thrust medio          = %.6f N\n", ...
    best.thrust_mean);



fprintf("\n");

fprintf("IMPULSO TOTALE        = %.6e N s\n", ...
    best_impulse);


%% ============================================================
%  14. SIMULAZIONE FINALE DELLA CONFIGURAZIONE OTTIMA
%
%  Riattivo plot, animazione e stampe.
% ============================================================

design_final = design;

design_final.O_F = best.OF_opt;

design_final.GOX = best.GOX_opt;

design_final.thrust = best.thrust0_opt;

design_final.radius_factor = best.radius_factor_opt;

design_final.n_tips = best.n_tips;


%% ------------------------------------------------------------
% Predesign finale
% -------------------------------------------------------------

pre_final = run_predesign( ...
    design_final);


%% ------------------------------------------------------------
% Settings finali
% -------------------------------------------------------------

settings_final = settings;

settings_final.make_plots = true;

settings_final.make_animation = true;

settings_final.verbose = true;


%% ------------------------------------------------------------
% Simulazione finale
% -------------------------------------------------------------

sim_final = run_temporal_simulation( ...
    pre_final, ...
    settings_final);


%% ============================================================
%  15. RISULTATI FINALI
% ============================================================

fprintf("\n");
fprintf("=============================================\n");
fprintf("          RISULTATI SIMULAZIONE\n");
fprintf("=============================================\n");

fprintf("Tempo finale:            %.6f s\n", ...
    sim_final.t(end));

fprintf("\n");

fprintf("Pressione iniziale:      %.6f bar\n", ...
    sim_final.p(1)*1e-5);

fprintf("Pressione finale:        %.6f bar\n", ...
    sim_final.p(end)*1e-5);

fprintf("\n");

fprintf("GOX iniziale:            %.6f kg/(m^2 s)\n", ...
    sim_final.Gox(1));

fprintf("GOX finale:              %.6f kg/(m^2 s)\n", ...
    sim_final.Gox(end));

fprintf("\n");

fprintf("O/F iniziale:            %.6f\n", ...
    sim_final.OF(1));

fprintf("O/F finale:              %.6f\n", ...
    sim_final.OF(end));

fprintf("\n");

fprintf("Spinta media:            %.6f N\n", ...
    sim_final.thrust_mean);

fprintf("Impulso totale:          %.6e N s\n", ...
    sim_final.total_impulse);


%% ============================================================
%  16. PLOT SPINTA
% ============================================================

figure;

plot( ...
    sim_final.t, ...
    sim_final.thrust, ...
    "LineWidth", ...
    1.8);

hold on
grid on


xlabel("Tempo [s]");

ylabel("Thrust [N]");

title( ...
    sprintf( ...
    "Thrust - configurazione ottima (%d punte)", ...
    best.n_tips));


%% ============================================================
%  17. PLOT IMPULSO VS NUMERO DI PUNTE
% ============================================================

figure;

plot( ...
    n_tips_values, ...
    impulses, ...
    "o-", ...
    "LineWidth", ...
    1.8);

grid on

xlabel("Numero di punte");

ylabel("Impulso totale [N s]");

title("Impulso totale ottimizzato");

toc

%% ============================================================
%  FUNZIONI LOCALI
% ============================================================


%% ============================================================
%  VALUTAZIONE DESIGN CON CACHE
% ============================================================

function sim = evaluate_design_cached( ...
    x, ...
    n_tips, ...
    design_base, ...
    settings)

    persistent last_x
    persistent last_n_tips
    persistent last_sim

    %% ========================================================
    %  CACHE
    % ========================================================

    if ~isempty(last_x) && ...
       ~isempty(last_n_tips) && ...
       ~isempty(last_sim)

        if isequal(x, last_x) && ...
           isequal(n_tips, last_n_tips)

            sim = last_sim;
            return

        end

    end

    %% ========================================================
    %  COSTRUZIONE DESIGN
    % ========================================================

    design = design_base;

    design.O_F = x(1);
    design.GOX = x(2);
    design.thrust = x(3);
    design.radius_factor = x(4);
    design.n_tips = n_tips;

    %% ========================================================
    %  PREDESIGN
    % ========================================================

    pre = run_predesign(design);

    %% ========================================================
    %  SIMULAZIONE TEMPORALE
    %
    %  Nessun try/catch.
    % ========================================================

    sim = run_temporal_simulation( ...
        pre, ...
        settings);

    %% ========================================================
    %  AGGIORNAMENTO CACHE
    % ========================================================

    last_x = x;
    last_n_tips = n_tips;
    last_sim = sim;

end


%% ============================================================
%  FUNZIONE OBIETTIVO
% ============================================================

function J = objective_function( ...
    x, ...
    n_tips, ...
    design_base, ...
    settings)

    sim = evaluate_design_cached( ...
        x, ...
        n_tips, ...
        design_base, ...
        settings);

    %% ========================================================
    %  MASSIMIZZAZIONE IMPULSO TOTALE
    % ========================================================

    J = -sim.total_impulse;

end


%% ============================================================
%  VINCOLI NON LINEARI
% ============================================================

function [c, ceq] = optimization_constraints( ...
    x, ...
    n_tips, ...
    design_base, ...
    settings, ...
    thrust_mean_min, ...
    OF_sim_min, ...
    OF_sim_max, ...
    pc_sim_min, ...
    pc_sim_max)


    sim = evaluate_design_cached( ...
        x, ...
        n_tips, ...
        design_base, ...
        settings);

    OF_valid = sim.OF(isfinite(sim.OF));

    p_valid = sim.p(isfinite(sim.p));


    if isempty(OF_valid) || isempty(p_valid)

        c = 1e6 * ones(5,1);

        ceq = [];

        return

    end


    %% Estremi simulazione

    OF_min_actual = min(OF_valid);
    OF_max_actual = max(OF_valid);

    pc_min_actual = min(p_valid);
    pc_max_actual = max(p_valid);


    %% ========================================================
    % VINCOLI
    %
    % fmincon:
    %
    %       c <= 0
    % ========================================================


    % Fmean >= Fmean_min
    c_thrust = ...
        thrust_mean_min ...
        - sim.thrust_mean;


    % OF <= OFmax
    c_OF_max = ...
        OF_max_actual ...
        - OF_sim_max;


    % OF >= OFmin
    c_OF_min = ...
        OF_sim_min ...
        - OF_min_actual;


    % pc <= pcmax
    c_pc_max = ...
        pc_max_actual ...
        - pc_sim_max;


    % pc >= pcmin
    c_pc_min = ...
        pc_sim_min ...
        - pc_min_actual;


    c = [
        c_thrust
        c_OF_max
        c_OF_min
        c_pc_max
        c_pc_min
    ];


    ceq = [];

end