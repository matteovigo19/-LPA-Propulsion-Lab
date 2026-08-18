close all
clear
clc

%% ============================================================
%  MAIN OPT
%
%  Geometrie supportate:
%       design.type = "star"
%       design.type = "cylinder"
%
%  STAR:
%       x(1) = O/F iniziale
%       x(2) = GOX iniziale
%       x(3) = thrust iniziale
%       x(4) = radius factor
%       + ciclo esterno su n_tips
%
%  CYLINDER:
%       x(1) = O/F iniziale
%       x(2) = GOX iniziale
%       x(3) = thrust iniziale
%
%  Il diametro cilindrico iniziale viene determinato nel predesign.
% ============================================================

tic;

%% ============================================================
%  1. DATI GENERALI DEL DESIGN
% ============================================================

design.eps = 200;
design.pch_bar = 20;

%% ============================================================
%  2. COMBUSTIBILE
% ============================================================

design.rho_f = 920;
design.a_rf = 0.027;
design.n_rf = 0.75;

%% ============================================================
%  3. SCELTA GEOMETRIA
% ============================================================

% "star" oppure "cylinder"
design.type = "cylinder";

% ============================================================
% SCELTA FUNZIONE OBIETTIVO
%
%   "impulse" -> massimizza l'impulso totale
%   "shift"   -> minimizza gli shift relativi di O/F e GOX
% ============================================================

design.opt = "shift";

% ============================================================
% OUTPUT E PLOT
%
% Questi booleani sono indipendenti da design.opt.
% Quindi è possibile ottimizzare una metrica e stampare/plottare
% anche l'altra.
% ============================================================

show_impulse_output = true;
show_impulse_plots  = true;

show_shift_output = true;
show_shift_plots  = true;

%% ============================================================
%  4. LIMITI VARIABILI DI OTTIMIZZAZIONE
% ============================================================

OF_min = 3;
OF_max_global = 4.5;

% SOLO STAR: hn-òpopmassimo O/F iniziale dipendente dal numero di punte.
OF_max_by_tips = containers.Map( ...
    [7, 8, 9, 10, 12, 14], ...
    [3.7, 3.4, 3.1, 2.7, 2.4, 2.1]);

GOX_min = 500;
GOX_max = 700;

thrust_min = 45000;
thrust_max = 55000;

% pesi shift
w_OF  = 0.5;
w_GOX = 0.5;

% SOLO STAR
radius_factor_min = 0.5;
radius_factor_max = 0.8;

% SOLO STAR
n_tips_min = 4;
n_tips_max = 4;
n_tips_values = n_tips_min:n_tips_max;

%% ============================================================
%  5. SETTINGS SIMULAZIONE TEMPORALE
% ============================================================

settings = struct();

settings.tmax = 300;
settings.ext_diameter = 1000.1*2e-3;

% STAR
settings.multiplyer = 10;

% CYLINDER
settings.n_cylinder = 70;
  
settings.RelTol = 1e-8;
settings.AbsTol = 1e-10;
settings.fine_ode_boolean = false;

settings.plot_case = "ode45";
settings.n_output = 10;

settings.make_plots = false;
settings.make_animation = false;
settings.verbose = false;

%% ============================================================
%  6. VINCOLI
% ============================================================

thrust_mean_max = 50000;

OF_sim_min = 1.15;
OF_sim_max = 19.45;

pc_sim_min_bar = 101325 * 1e-5;
pc_sim_max_bar = 99;

pc_sim_min = pc_sim_min_bar * 1e5;
pc_sim_max = pc_sim_max_bar * 1e5;

%% ============================================================
%  7. OPZIONI FMINCON
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
%  8. INIZIALIZZAZIONE
% ============================================================

rng("shuffle");

n_guess = 20;
results = struct([]);

%% ============================================================
%  9. OTTIMIZZAZIONE
% ============================================================

switch lower(string(design.type))

    %% ========================================================
    %  STAR
    % ========================================================
    case "star"

        n_configurations = length(n_tips_values);

        for ii = 1:n_configurations

            n_tips = n_tips_values(ii);

            % O/F max locale
            if isKey(OF_max_by_tips, n_tips)
                OF_max_local = OF_max_by_tips(n_tips);
            else
                OF_max_local = OF_max_global;
            end

            % x = [OF, GOX, thrust, radius_factor]
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
            fprintf("       OBIETTIVO: %s\n", upper(string(design.opt)));
            fprintf("=====================================================\n");

            best_fval_local = Inf;
            best_x_opt_local = [];
            best_x0_local = [];
            best_exitflag_local = [];
            best_output_local = [];
            best_pre_local = [];
            best_sim_local = [];

            for jj = 1:n_guess

                fprintf("\n");
                fprintf("---------------------------------------------\n");
                fprintf("Guess %d / %d - STAR %d punte\n", ...
                    jj, n_guess, n_tips);
                fprintf("---------------------------------------------\n");

                x0 = lb_local + ...
                    rand(size(lb_local)).*(ub_local-lb_local);

                fprintf("Range O/F       = [%.4f, %.4f]\n", ...
                    OF_min, OF_max_local);
                fprintf("O/F             = %.4f\n", x0(1));
                fprintf("GOX             = %.2f kg/(m^2 s)\n", x0(2));
                fprintf("Thrust          = %.2f N\n", x0(3));
                fprintf("Radius factor   = %.4f\n", x0(4));

                geometry_id = n_tips;

                objective = make_objective( ...
                    design.opt, ...
                    geometry_id, ...
                    design, ...
                    settings, ...
                    w_OF, ...
                    w_GOX);

                nonlcon = @(x) optimization_constraints( ...
                    x, geometry_id, design, settings, ...
                    thrust_mean_max, ...
                    OF_sim_min, OF_sim_max, ...
                    pc_sim_min, pc_sim_max);

                [x_opt, fval, exitflag, output] = fmincon( ...
                    objective, ...
                    x0, ...
                    [], [], [], [], ...
                    lb_local, ...
                    ub_local, ...
                    nonlcon, ...
                    options);

                design_opt = build_design_from_x( ...
                    x_opt, geometry_id, design);

                pre_opt = run_predesign(design_opt);

                sim_opt = run_temporal_simulation( ...
                    pre_opt, settings);

                fprintf("\nRisultato guess %d:\n", jj);
                fprintf("O/F             = %.6f\n", x_opt(1));
                fprintf("GOX             = %.6f kg/(m^2 s)\n", x_opt(2));
                fprintf("Thrust iniziale = %.6f N\n", x_opt(3));
                fprintf("Radius factor   = %.6f\n", x_opt(4));
                fprintf("Thrust medio    = %.6f N\n", sim_opt.thrust_mean);
                fprintf("Impulso totale  = %.6e N s\n", sim_opt.total_impulse);
                fprintf("J (%s)          = %.6e\n", upper(string(design.opt)), fval);

                if fval < best_fval_local

                    best_fval_local = fval;
                    best_x_opt_local = x_opt;
                    best_x0_local = x0;
                    best_exitflag_local = exitflag;
                    best_output_local = output;
                    best_pre_local = pre_opt;
                    best_sim_local = sim_opt;

                    fprintf("\n>>> Nuovo miglior risultato STAR %d punte <<<\n", ...
                        n_tips);

                end

            end

            results(ii).type = "star";
            results(ii).n_tips = n_tips;
            results(ii).n_guess = n_guess;
            results(ii).x0 = best_x0_local;
            results(ii).x_opt = best_x_opt_local;

            results(ii).OF_opt = best_x_opt_local(1);
            results(ii).GOX_opt = best_x_opt_local(2);
            results(ii).thrust0_opt = best_x_opt_local(3);
            results(ii).radius_factor_opt = best_x_opt_local(4);

            results(ii).thrust_mean = best_sim_local.thrust_mean;
            results(ii).total_impulse = best_sim_local.total_impulse;

            results(ii).fval = best_fval_local;
            results(ii).exitflag = best_exitflag_local;
            results(ii).output = best_output_local;

            results(ii).pre = best_pre_local;
            results(ii).sim = best_sim_local;

            fprintf("\n");
            fprintf("=====================================================\n");
            fprintf(" MIGLIOR RISULTATO STAR %d PUNTE SU %d GUESS\n", ...
                n_tips, n_guess);
            fprintf("=====================================================\n");

            fprintf("O/F                   = %.6f\n", results(ii).OF_opt);
            fprintf("GOX                   = %.6f kg/(m^2 s)\n", results(ii).GOX_opt);
            fprintf("Thrust iniziale       = %.6f N\n", results(ii).thrust0_opt);
            fprintf("Radius factor         = %.6f\n", results(ii).radius_factor_opt);
            fprintf("Thrust medio          = %.6f N\n", results(ii).thrust_mean);
            fprintf("Thrust medio massimo  = %.6f N\n", thrust_mean_max);
            fprintf("Impulso totale        = %.6e N s\n", results(ii).total_impulse);

        end


    %% ========================================================
    %  CYLINDER
    %
    %  Nessun n_tips.
    %  Nessun radius_factor.
    %  x = [OF, GOX, thrust]
    % ========================================================
    case "cylinder"

        n_configurations = 1;

        lb_local = [
            OF_min
            GOX_min
            thrust_min
        ];

        ub_local = [
            OF_max_global
            GOX_max
            thrust_max
        ];

        fprintf("\n");
        fprintf("=====================================================\n");
        fprintf("          OTTIMIZZAZIONE CYLINDER\n");
        fprintf("          OBIETTIVO: %s\n", upper(string(design.opt)));
        fprintf("=====================================================\n");

        best_fval_local = Inf;
        best_x_opt_local = [];
        best_x0_local = [];
        best_exitflag_local = [];
        best_output_local = [];
        best_pre_local = [];
        best_sim_local = [];

        % =====================================================
        %  SALVATAGGIO RISULTATI DI TUTTI I GUESS - CYLINDER
        %
        %  Servono per confrontare l'andamento dei parametri
        %  ottimizzati e dell'impulso totale tra i vari guess.
        % =====================================================

        cylinder_guess_OF = NaN(n_guess,1);
        cylinder_guess_GOX = NaN(n_guess,1);
        cylinder_guess_impulse = NaN(n_guess,1);

        for jj = 1:n_guess

            fprintf("\n");
            fprintf("---------------------------------------------\n");
            fprintf("Guess %d / %d - CYLINDER\n", jj, n_guess);
            fprintf("---------------------------------------------\n");

            x0 = lb_local + ...
                rand(size(lb_local)).*(ub_local-lb_local);

            fprintf("Range O/F       = [%.4f, %.4f]\n", ...
                OF_min, OF_max_global);
            fprintf("O/F             = %.4f\n", x0(1));
            fprintf("GOX             = %.2f kg/(m^2 s)\n", x0(2));
            fprintf("Thrust          = %.2f N\n", x0(3));

            geometry_id = [];

            objective = make_objective( ...
                design.opt, ...
                geometry_id, ...
                design, ...
                settings, ...
                w_OF, ...
                w_GOX);

            nonlcon = @(x) optimization_constraints( ...
                x, geometry_id, design, settings, ...
                thrust_mean_max, ...
                OF_sim_min, OF_sim_max, ...
                pc_sim_min, pc_sim_max);

            [x_opt, fval, exitflag, output] = fmincon( ...
                objective, ... 
                x0, ...
                [], [], [], [], ...
                lb_local, ...
                ub_local, ...
                nonlcon, ...
                options);

            design_opt = build_design_from_x( ...
                x_opt, geometry_id, design);

            pre_opt = run_predesign(design_opt);

            sim_opt = run_temporal_simulation( ...
                pre_opt, settings);

            % =================================================
            % Salvataggio del risultato di questo guess
            % =================================================

            cylinder_guess_OF(jj) = x_opt(1);
            cylinder_guess_GOX(jj) = x_opt(2);
            cylinder_guess_impulse(jj) = sim_opt.total_impulse;

            fprintf("\nRisultato guess %d:\n", jj);
            fprintf("O/F             = %.6f\n", x_opt(1));
            fprintf("GOX             = %.6f kg/(m^2 s)\n", x_opt(2));
            fprintf("Thrust iniziale = %.6f N\n", x_opt(3));
            fprintf("Thrust medio    = %.6f N\n", sim_opt.thrust_mean);
            fprintf("Impulso totale  = %.6e N s\n", sim_opt.total_impulse);
            fprintf("J (%s)          = %.6e\n", upper(string(design.opt)), fval);

            if fval < best_fval_local

                best_fval_local = fval;
                best_x_opt_local = x_opt;
                best_x0_local = x0;
                best_exitflag_local = exitflag;
                best_output_local = output;
                best_pre_local = pre_opt;
                best_sim_local = sim_opt;

                fprintf("\n>>> Nuovo miglior risultato CYLINDER <<<\n");

            end

        end

        results(1).type = "cylinder";
        results(1).n_guess = n_guess;
        results(1).x0 = best_x0_local;
        results(1).x_opt = best_x_opt_local;

        results(1).OF_opt = best_x_opt_local(1);
        results(1).GOX_opt = best_x_opt_local(2);
        results(1).thrust0_opt = best_x_opt_local(3);

        results(1).thrust_mean = best_sim_local.thrust_mean;
        results(1).total_impulse = best_sim_local.total_impulse;

        results(1).fval = best_fval_local;
        results(1).exitflag = best_exitflag_local;
        results(1).output = best_output_local;

        results(1).pre = best_pre_local;
        results(1).sim = best_sim_local;

        % Risultati di tutti i guess del cilindro
        results(1).guess_OF = cylinder_guess_OF;
        results(1).guess_GOX = cylinder_guess_GOX;
        results(1).guess_impulse = cylinder_guess_impulse;

        fprintf("\n");
        fprintf("=====================================================\n");
        fprintf(" MIGLIOR RISULTATO CYLINDER SU %d GUESS\n", n_guess);
        fprintf("=====================================================\n");

        fprintf("O/F                   = %.6f\n", results(1).OF_opt);
        fprintf("GOX                   = %.6f kg/(m^2 s)\n", results(1).GOX_opt);
        fprintf("Thrust iniziale       = %.6f N\n", results(1).thrust0_opt);
        fprintf("Thrust medio          = %.6f N\n", results(1).thrust_mean);
        fprintf("Thrust medio massimo  = %.6f N\n", thrust_mean_max);
        fprintf("Impulso totale        = %.6e N s\n", results(1).total_impulse);


    otherwise

        error( ...
            "Tipo di geometria non riconosciuto: %s. Usare 'star' o 'cylinder'.", ...
            design.type);

end


%% ============================================================
%  10. CONFIGURAZIONE MIGLIORE
%
%  La configurazione migliore è quella che minimizza la funzione
%  obiettivo selezionata in design.opt.
% ============================================================

objective_values = [results.fval];

[best_objective, idx_best] = min(objective_values);

best = results(idx_best);

best_impulse = best.total_impulse;


%% ============================================================
%  11. RISULTATI OTTIMIZZAZIONE
% ============================================================

fprintf("\n");
fprintf("=====================================================\n");
fprintf("             CONFIGURAZIONE MIGLIORE\n");
fprintf("=====================================================\n");

fprintf("Geometria             = %s\n", design.type);
fprintf("Funzione obiettivo    = %s\n", design.opt);
fprintf("Valore obiettivo      = %.8e\n", best_objective);

if lower(string(design.type)) == "star"
    fprintf("Numero punte          = %d\n", best.n_tips);
end

fprintf("O/F iniziale          = %.6f\n", best.OF_opt);
fprintf("GOX iniziale          = %.6f kg/(m^2 s)\n", best.GOX_opt);
fprintf("Thrust iniziale       = %.6f N\n", best.thrust0_opt);

if lower(string(design.type)) == "star"
    fprintf("Radius factor         = %.6f\n", best.radius_factor_opt);
end


%% ============================================================
%  12. SIMULAZIONE FINALE
% ============================================================

design_final = design;

design_final.O_F = best.OF_opt;
design_final.GOX = best.GOX_opt;
design_final.thrust = best.thrust0_opt;

if lower(string(design.type)) == "star"
    design_final.radius_factor = best.radius_factor_opt;
    design_final.n_tips = best.n_tips;
end

pre_final = run_predesign(design_final);

settings_final = settings;

% I plot specifici vengono gestiti separatamente dal main.
settings_final.make_plots = false;
settings_final.make_animation = false;
settings_final.verbose = false;

sim_final = run_temporal_simulation( ...
    pre_final, ...
    settings_final);


%% ============================================================
%  13. CALCOLO METRICHE FINALI
%
%  Entrambe vengono sempre calcolate, indipendentemente da
%  design.opt.
% ============================================================

impulse_metrics = evaluate_impulse_metrics(sim_final);

shift_metrics = evaluate_shift_metrics( ...
    sim_final, ...
    w_OF, ...
    w_GOX);


%% ============================================================
%  14. OUTPUT IMPULSO
% ============================================================

if show_impulse_output

    fprintf("\n");
    fprintf("=====================================================\n");
    fprintf("              OUTPUT IMPULSO\n");
    fprintf("=====================================================\n");

    fprintf("Tempo finale           = %.6f s\n", ...
        impulse_metrics.duration);

    fprintf("Thrust medio           = %.6f N\n", ...
        impulse_metrics.thrust_mean);

    fprintf("Thrust iniziale sim.   = %.6f N\n", ...
        impulse_metrics.thrust_initial);

    fprintf("Thrust finale sim.     = %.6f N\n", ...
        impulse_metrics.thrust_final);

    fprintf("Impulso totale         = %.6e N s\n", ...
        impulse_metrics.total_impulse);

end


%% ============================================================
%  15. OUTPUT SHIFT
% ============================================================

if show_shift_output

    fprintf("\n");
    fprintf("=====================================================\n");
    fprintf("               OUTPUT SHIFT\n");
    fprintf("=====================================================\n");

    fprintf("O/F iniziale           = %.6f\n", ...
        shift_metrics.OF0);

    fprintf("O/F finale             = %.6f\n", ...
        shift_metrics.OF_final);

    fprintf("GOX iniziale           = %.6f kg/(m^2 s)\n", ...
        shift_metrics.GOX0);

    fprintf("GOX finale             = %.6f kg/(m^2 s)\n", ...
        shift_metrics.GOX_final);

    fprintf("RMS shift O/F          = %.6f %%\n", ...
        100*shift_metrics.RMS_OF);

    fprintf("RMS shift GOX          = %.6f %%\n", ...
        100*shift_metrics.RMS_GOX);

    fprintf("MAX shift O/F          = %.6f %%\n", ...
        100*shift_metrics.MAX_OF);

    fprintf("MAX shift GOX          = %.6f %%\n", ...
        100*shift_metrics.MAX_GOX);

    fprintf("J shift O/F            = %.8e\n", ...
        shift_metrics.J_OF);

    fprintf("J shift GOX            = %.8e\n", ...
        shift_metrics.J_GOX);

    fprintf("J shift complessivo    = %.8e\n", ...
        shift_metrics.J_total);

end


%% ============================================================
%  16. PLOT IMPULSO
% ============================================================

if show_impulse_plots

    figure;

    plot( ...
        sim_final.t, ...
        sim_final.thrust, ...
        "LineWidth", ...
        1.8);

    grid on

    xlabel("Tempo [s]");
    ylabel("Thrust [N]");

    title(sprintf( ...
        "Thrust nel tempo - %s - opt: %s", ...
        upper(string(design.type)), ...
        upper(string(design.opt))));


    figure;

    plot( ...
        sim_final.t, ...
        impulse_metrics.cumulative_impulse, ...
        "LineWidth", ...
        1.8);

    grid on

    xlabel("Tempo [s]");
    ylabel("Impulso cumulativo [N s]");

    title(sprintf( ...
        "Impulso cumulativo - %s - opt: %s", ...
        upper(string(design.type)), ...
        upper(string(design.opt))));

end


%% ============================================================
%  17. PLOT SHIFT
% ============================================================

if show_shift_plots

    figure;

    plot( ...
        sim_final.t, ...
        shift_metrics.OF_relative, ...
        "LineWidth", ...
        1.8);

    hold on
    yline(1, "--");

    grid on

    xlabel("Tempo [s]");
    ylabel("(O/F)/(O/F)_0 [-]");

    title(sprintf( ...
        "Evoluzione relativa O/F - %s - opt: %s", ...
        upper(string(design.type)), ...
        upper(string(design.opt))));


    figure;

    plot( ...
        sim_final.t, ...
        shift_metrics.GOX_relative, ...
        "LineWidth", ...
        1.8);

    hold on
    yline(1, "--");

    grid on

    xlabel("Tempo [s]");
    ylabel("GOX/GOX_0 [-]");

    title(sprintf( ...
        "Evoluzione relativa GOX - %s - opt: %s", ...
        upper(string(design.type)), ...
        upper(string(design.opt))));


    figure;

    plot( ...
        sim_final.t, ...
        100*shift_metrics.shift_OF, ...
        "LineWidth", ...
        1.8);

    hold on

    plot( ...
        sim_final.t, ...
        100*shift_metrics.shift_GOX, ...
        "LineWidth", ...
        1.8);

    yline(0, "--");

    grid on

    xlabel("Tempo [s]");
    ylabel("Shift relativo [%]");

    title(sprintf( ...
        "Shift O/F e GOX - %s - opt: %s", ...
        upper(string(design.type)), ...
        upper(string(design.opt))));

    legend( ...
        "O/F", ...
        "GOX", ...
        "Location", ...
        "best");

end


%% ============================================================
%  18. CONFRONTO CONFIGURAZIONI STAR
% ============================================================

if lower(string(design.type)) == "star" && length(results) > 1

    figure;

    plot( ...
        [results.n_tips], ...
        [results.fval], ...
        "o-", ...
        "LineWidth", ...
        1.8);

    grid on

    xlabel("Numero di punte");
    ylabel("Valore funzione obiettivo");

    title(sprintf( ...
        "Obiettivo ottimizzato vs numero di punte - %s", ...
        upper(string(design.opt))));

end

toc


%% ============================================================
%  FUNZIONI LOCALI
% ============================================================

function design = build_design_from_x( ...
    x, ...
    geometry_id, ...
    design_base)

    design = design_base;

    switch lower(string(design.type))

        case "star"

            design.O_F = x(1);
            design.GOX = x(2);
            design.thrust = x(3);
            design.radius_factor = x(4);

            design.n_tips = geometry_id;


        case "cylinder"

            design.O_F = x(1);
            design.GOX = x(2);
            design.thrust = x(3);

            % Nessun radius_factor.
            % Nessun n_tips.
            % Il diametro viene determinato dal predesign.


        otherwise

            error( ...
                "Tipo di geometria non riconosciuto: %s", ...
                design.type);

    end

end


function sim = evaluate_design_cached( ...
    x, ...
    geometry_id, ...
    design_base, ...
    settings)

    persistent last_x
    persistent last_geometry_id
    persistent last_type
    persistent last_sim

    current_type = lower(string(design_base.type));

    if ~isempty(last_x) && ...
       ~isempty(last_type) && ...
       ~isempty(last_sim)

        same_x = isequal(x, last_x);
        same_type = isequal(current_type, last_type);
        same_geometry = isequal(geometry_id, last_geometry_id);

        if same_x && same_type && same_geometry

            sim = last_sim;
            return

        end

    end

    design = build_design_from_x( ...
        x, geometry_id, design_base);

    pre = run_predesign(design);

    sim = run_temporal_simulation( ...
        pre, settings);

    last_x = x;
    last_geometry_id = geometry_id;
    last_type = current_type;
    last_sim = sim;

end


function objective = make_objective( ...
    opt_type, ...
    geometry_id, ...
    design_base, ...
    settings, ...
    w_OF, ...
    w_GOX)

    switch lower(string(opt_type))

        case "impulse"

            objective = @(x) objective_impulse( ...
                x, ...
                geometry_id, ...
                design_base, ...
                settings);


        case "shift"

            objective = @(x) objective_shift( ...
                x, ...
                geometry_id, ...
                design_base, ...
                settings, ...
                w_OF, ...
                w_GOX);


        otherwise

            error( ...
                "Funzione obiettivo non riconosciuta: %s. Usare 'impulse' o 'shift'.", ...
                opt_type);

    end

end


function J = objective_impulse( ...
    x, ...
    geometry_id, ...
    design_base, ...
    settings)

    sim = evaluate_design_cached( ...
        x, ...
        geometry_id, ...
        design_base, ...
        settings);

    J = -sim.total_impulse / 15e6;

end


function J = objective_shift( ...
    x, ...
    geometry_id, ...
    design_base, ...
    settings, ...
    w_OF, ...
    w_GOX)

    sim = evaluate_design_cached( ...
        x, ...
        geometry_id, ...
        design_base, ...
        settings);

    metrics = evaluate_shift_metrics( ...
        sim, ...
        w_OF, ...
        w_GOX);

    J = metrics.J_total;

end


function [c, ceq] = optimization_constraints( ...
    x, ...
    geometry_id, ...
    design_base, ...
    settings, ...
    thrust_mean_max, ...
    OF_sim_min, ...
    OF_sim_max, ...
    pc_sim_min, ...
    pc_sim_max)

    sim = evaluate_design_cached( ...
        x, geometry_id, design_base, settings);

    OF_valid = sim.OF(isfinite(sim.OF));
    p_valid = sim.p(isfinite(sim.p));

    if isempty(OF_valid) || isempty(p_valid)

        c = 1e6 * ones(5,1);
        ceq = [];
        return

    end

    OF_min_actual = min(OF_valid);
    OF_max_actual = max(OF_valid);

    pc_min_actual = min(p_valid);
    pc_max_actual = max(p_valid);

    % fmincon richiede c <= 0

    % Fmean <= Fmean_max
    c_thrust = ...
        sim.thrust_mean ...
        - thrust_mean_max;

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

%% ============================================================
%  METRICHE IMPULSO
% ============================================================

function metrics = evaluate_impulse_metrics(sim)

    t = sim.t(:);
    thrust = sim.thrust(:);

    metrics = struct();

    metrics.duration = t(end) - t(1);
    metrics.thrust_initial = thrust(1);
    metrics.thrust_final = thrust(end);
    metrics.thrust_mean = sim.thrust_mean;
    metrics.total_impulse = sim.total_impulse;

    metrics.cumulative_impulse = ...
        cumtrapz(t, thrust);

end


%% ============================================================
%  METRICHE SHIFT
% ============================================================

function metrics = evaluate_shift_metrics( ...
    sim, ...
    w_OF, ...
    w_GOX)

    t = sim.t(:);

    OF = sim.OF(:);
    GOX = sim.Gox(:);

    OF0 = OF(1);
    GOX0 = GOX(1);

    shift_OF = ...
        (OF - OF0) / OF0;

    shift_GOX = ...
        (GOX - GOX0) / GOX0;

    T = t(end) - t(1);

    if T <= 0

        RMS_OF = Inf;
        RMS_GOX = Inf;

    else

        RMS_OF = sqrt( ...
            trapz(t, shift_OF.^2) / T);

        RMS_GOX = sqrt( ...
            trapz(t, shift_GOX.^2) / T);

    end

    MAX_OF = max(abs(shift_OF));
    MAX_GOX = max(abs(shift_GOX));

    J_OF = RMS_OF;
    J_GOX = RMS_GOX;

    J_total = ...
        w_OF * J_OF ...
        + ...
        w_GOX * J_GOX;

    metrics = struct();

    metrics.OF0 = OF0;
    metrics.OF_final = OF(end);

    metrics.GOX0 = GOX0;
    metrics.GOX_final = GOX(end);

    metrics.shift_OF = shift_OF;
    metrics.shift_GOX = shift_GOX;

    metrics.OF_relative = OF / OF0;
    metrics.GOX_relative = GOX / GOX0;

    metrics.RMS_OF = RMS_OF;
    metrics.RMS_GOX = RMS_GOX;

    metrics.MAX_OF = MAX_OF;
    metrics.MAX_GOX = MAX_GOX;

    metrics.J_OF = J_OF;
    metrics.J_GOX = J_GOX;
    metrics.J_total = J_total;

end


%% ============================================================
%  NORMALIZZAZIONE PER PLOT
%
%  Porta un vettore nell'intervallo [0,1].
%  Se tutti i valori sono uguali restituisce un vettore di 0.5.
% ============================================================

function y_norm = normalize_for_plot(y)

    y = y(:);

    y_min = min(y, [], "omitnan");
    y_max = max(y, [], "omitnan");

    if ~isfinite(y_min) || ~isfinite(y_max)

        y_norm = NaN(size(y));
        return

    end

    if abs(y_max - y_min) < eps(max(abs([y_min y_max 1])))

        y_norm = 0.5 * ones(size(y));

    else

        y_norm = (y - y_min) / (y_max - y_min);

    end

end

