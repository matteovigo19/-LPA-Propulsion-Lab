%% ============================================================
%  CONFRONTO OSSIDANTI - HTPB
%
%  Sam Kobliska CEA - 02/07/2017
%
%  Output:
%       Isp = f(O/F)
%
%  Ossidanti:
%       LOX
%       GOX
%       H2O2
%       N2O
%       N2O4
%
%  Ogni ossidante ha il proprio vettore O/F
% =============================================================

clear
close all
clc

%% ============================================================
% DATI DEL PROBLEMA
% ============================================================

% Pressione di camera
pc = 20;                    % bar

% Rapporto di espansione dell'ugello
epsilon = 200;              % Ae/At

% Temperatura iniziale combustibile
T_fuel = 300;               % K

%% ============================================================
% DEFINIZIONE DEGLI OSSIDANTI
% ============================================================

% Nome visualizzato nel grafico
nome_ossidante = { ...
    'LOX', ...
    'GOX', ...
    'H_2O_2', ...
    'N_2O', ...
    'N_2O_4'};

% Nome della specie nel database CEA
ossidante_CEA = { ...
    'O2', ...
    'O2', ...
    'H2O2', ...
    'N2O', ...
    'N2O4'};

% Temperatura iniziale dell'ossidante
T_iniziale_ox = [ ...
    161, ...                 % LOX
    300, ...                 % GOX
    300, ...                 % H2O2
    298, ...                 % N2O
    300];                    % N2O4

% Numero di ossidanti
N_ossidanti = length(nome_ossidante);

%% ============================================================
% VETTORI O/F
% ============================================================

OF_LOX  = 0.5:0.1:9;
OF_GOX  = 0.5:0.1:9;
OF_H2O2 = 0.5:0.1:9;
OF_N2O  = 2.5:0.1:9;
OF_N2O4 = 0.5:0.1:9;

% Raccolta dei vettori
OF_all = { ...
    OF_LOX, ...
    OF_GOX, ...
    OF_H2O2, ...
    OF_N2O, ...
    OF_N2O4};

%% ============================================================
% PREALLOCAMENTO
% ============================================================

Isp = cell(N_ossidanti,1);

%% ============================================================
% CICLO CEA
% ============================================================

for i = 1:N_ossidanti

    % Vettore O/F dell'ossidante corrente
    OF_current = OF_all{i};

    % Preallocamento Isp
    Isp_current = zeros(size(OF_current));

    fprintf('\n');
    fprintf('============================================\n');
    fprintf(' Ossidante: %s\n',nome_ossidante{i});
    fprintf(' Temperatura: %.1f K\n',T_iniziale_ox(i));
    fprintf('============================================\n');

    for j = 1:length(OF_current)

        of = OF_current(j);

        fprintf('O/F = %.2f\n',of);

        %% -----------------------------------------------------
        % CEA
        % ------------------------------------------------------

        CEA_output = CEA( ...
            'problem','o/f',of, ...
            'rocket','eq', ...
            'p,bar',pc, ...
            'sup,ae/at',epsilon, ...
            'reactants', ...
                'oxid',ossidante_CEA{i}, ...
                'wt%',100, ...
                't,k',T_iniziale_ox(i), ...
                'fuel','HTPB', ...
                'wt%',100, ...
                't,k',T_fuel, ...
                'h,kj/mol',-58, ...
                'C',7.075, ...
                'H',10.650, ...
                'O',0.223, ...
                'N',0.063, ...
            'output', ...
                'short', ...
            'end');

        %% -----------------------------------------------------
        % DIAGNOSTICA CEA
        %
        % Viene eseguita solo per il primo punto:
        % LOX, O/F = 0.5
        %
        % Serve per verificare la struttura restituita
        % dalla versione di CEA utilizzata.
        % ------------------------------------------------------

       if i == 1 && j == 20

            fprintf('\n');
            fprintf('============================================================\n');
            fprintf('              DIAGNOSTICA CEA\n');
            fprintf('============================================================\n');

            fprintf('O/F = %.2f\n',of);
            fprintf('pc = %.2f bar\n',pc);
            fprintf('epsilon = %.2f\n',epsilon);

            fprintf('\n--- Campi disponibili in output.eql ---\n');
            disp(fieldnames(CEA_output.output.eql));

            fprintf('\n--- Struttura completa output.eql ---\n');
            disp(CEA_output.output.eql);

            fprintf('\n--- Vettore Isp ---\n');
            disp(CEA_output.output.eql.isp);

            fprintf('\n--- Pressione [bar] ---\n');
            disp(CEA_output.output.eql.pressure);
            
            fprintf('\n--- Temperatura [K] ---\n');
            disp(CEA_output.output.eql.temperature);
            
            fprintf('\n--- Mach [-] ---\n');
            disp(CEA_output.output.eql.mach);
            
            fprintf('\n--- Rapporto Ae/At [-] ---\n');
            disp(CEA_output.output.eql.aeat);
            
            fprintf('\n--- Isp [s] ---\n');
            disp(CEA_output.output.eql.isp);
            
            fprintf('\n--- Cstar [m/s] ---\n');
            disp(CEA_output.output.eql.cstar);


            fprintf('============================================================\n');
            fprintf('          FINE DIAGNOSTICA CEA\n');
            fprintf('============================================================\n');
            fprintf('\n');

        end

        %% -----------------------------------------------------
        % ESTRAZIONE Isp
        %
        % Il valore finale del vettore isp viene utilizzato
        % come Isp alla stazione finale dell'ugello.
        % ------------------------------------------------------

        Isp_current(j) = ...
            CEA_output.output.eql.isp(end);

    end

    % Salvataggio risultati dell'ossidante
    Isp{i} = Isp_current;

end

%% ============================================================
% PLOT Isp - O/F
% ============================================================

figure

hold on
grid on
box on

for i = 1:N_ossidanti

    plot( ...
        OF_all{i}, ...
        Isp{i}, ...
        'LineWidth',2, ...
        'DisplayName',nome_ossidante{i});

end

xlabel('O/F [-]')
ylabel('I_{sp} [s]')

title(sprintf( ...
    'HTPB - I_{sp} vs O/F - p_c = %.1f bar - \epsilon = %.0f', ...
    pc,epsilon))

legend( ...
    'show', ...
    'Location','best')

hold off

%% ============================================================
% MASSIMO Isp
% ============================================================

fprintf('\n\n');
fprintf('============================================================\n');
fprintf('                  RISULTATI Isp MASSIMO\n');
fprintf('============================================================\n');

fprintf('\n');
fprintf('%-10s  %-12s  %-15s\n', ...
    'Ossidante','O/F opt','Isp max [s]');

fprintf('------------------------------------------------------------\n');

for i = 1:N_ossidanti

    [Isp_max,index_max] = max(Isp{i});

    OF_opt = OF_all{i}(index_max);

    fprintf('%-10s  %-12.2f  %-15.2f\n', ...
        nome_ossidante{i}, ...
        OF_opt, ...
        Isp_max);

end

fprintf('============================================================\n');