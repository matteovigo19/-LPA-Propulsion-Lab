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

%% ============================================================
%  4. LIMITI VARIABILI DI OTTIMIZZAZIONE
% ============================================================

OF_min = 1.1;
OF_max_global = 2.0;

% SOLO STAR: massimo O/F iniziale dipendente dal numero di punte.
OF_max_by_tips = containers.Map( ...
    [7, 8, 9, 10, 12, 14], ...
    [3.7, 3.4, 3.1, 2.7, 2.4, 2.1]);

GOX_min = 300;
GOX_max = 800;

thrust_min = 45000;
thrust_max = 55000;

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

                objective = @(x) objective_function( ...
                    x, geometry_id, design, settings);

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
                fprintf("fval            = %.6e\n", fval);

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

            objective = @(x) objective_function( ...
                x, geometry_id, design, settings);

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
            fprintf("Thrust medio    = %.6f N\n", sim_opt.thrust_mean);
            fprintf("Impulso totale  = %.6e N s\n", sim_opt.total_impulse);
            fprintf("fval            = %.6e\n", fval);

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
% ============================================================

impulses = [results.total_impulse];

[best_impulse, idx_best] = max(impulses);

best = results(idx_best);


%% ============================================================
%  11. RISULTATI OTTIMIZZAZIONE
% ============================================================

fprintf("\n");
fprintf("=====================================================\n");
fprintf("             CONFIGURAZIONE MIGLIORE\n");
fprintf("=====================================================\n");

fprintf("Geometria             = %s\n", design.type);

if lower(string(design.type)) == "star"
    fprintf("Numero punte          = %d\n", best.n_tips);
end

fprintf("O/F iniziale          = %.6f\n", best.OF_opt);
fprintf("GOX iniziale          = %.6f kg/(m^2 s)\n", best.GOX_opt);
fprintf("Thrust iniziale       = %.6f N\n", best.thrust0_opt);

if lower(string(design.type)) == "star"
    fprintf("Radius factor         = %.6f\n", best.radius_factor_opt);
end

fprintf("Thrust medio          = %.6f N\n", best.thrust_mean);
fprintf("IMPULSO TOTALE        = %.6e N s\n", best_impulse);


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
settings_final.make_plots = true;
settings_final.make_animation = true;
settings_final.verbose = true;

sim_final = run_temporal_simulation( ...
    pre_final, settings_final);


%% ============================================================
%  13. RISULTATI FINALI
% ============================================================

fprintf("\n");
fprintf("=============================================\n");
fprintf("          RISULTATI SIMULAZIONE\n");
fprintf("=============================================\n");

fprintf("Tempo finale:            %.6f s\n", sim_final.t(end));
fprintf("Pressione iniziale:      %.6f bar\n", sim_final.p(1)*1e-5);
fprintf("Pressione finale:        %.6f bar\n", sim_final.p(end)*1e-5);
fprintf("GOX iniziale:            %.6f kg/(m^2 s)\n", sim_final.Gox(1));
fprintf("GOX finale:              %.6f kg/(m^2 s)\n", sim_final.Gox(end));
fprintf("O/F iniziale:            %.6f\n", sim_final.OF(1));
fprintf("O/F finale:              %.6f\n", sim_final.OF(end));
fprintf("Spinta media:            %.6f N\n", sim_final.thrust_mean);
fprintf("Impulso totale:          %.6e N s\n", sim_final.total_impulse);


%% ============================================================
%  14. PLOT SPINTA
% ============================================================

figure;

plot( ...
    sim_final.t, ...
    sim_final.thrust, ...
    "LineWidth", 1.8);

grid on

xlabel("Tempo [s]");
ylabel("Thrust [N]");

if lower(string(design.type)) == "star"

    title(sprintf( ...
        "Thrust - configurazione ottima STAR (%d punte)", ...
        best.n_tips));

else

    title("Thrust - configurazione ottima CYLINDER");

end


%% ============================================================
%  15. PLOT IMPULSO VS NUMERO DI PUNTE - SOLO STAR
% ============================================================

if lower(string(design.type)) == "star" && length(results) > 1

    figure;

    plot( ...
        [results.n_tips], ...
        impulses, ...
        "o-", ...
        "LineWidth", 1.8);

    grid on

    xlabel("Numero di punte");
    ylabel("Impulso totale [N s]");
    title("Impulso totale ottimizzato - STAR");

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


function J = objective_function( ...
    x, ...
    geometry_id, ...
    design_base, ...
    settings)

    sim = evaluate_design_cached( ...
        x, geometry_id, design_base, settings);

    J = -sim.total_impulse / 15e6;

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
