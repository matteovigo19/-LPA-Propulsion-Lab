%% ============================================================
% DEFINIZIONE DEGLI OSSIDANTI
% ============================================================

nome_ossidante = { ...
    'LOX', ...
    'GOX'};

ossidante_CEA = { ...
    'O2(L)', ...
    'O2'};

T_iniziale_ox = [ ...
    298, ...        % LOX
    298];             % GOX

N_ossidanti = length(nome_ossidante);
%% ============================================================
% CICLO CEA
% ============================================================

for i = 1:N_ossidanti

    % Vettore O/F relativo all'ossidante corrente
    OF_current = OF_all{i};

    % Preallocamento Isp
    Isp_current = NaN(size(OF_current));

    fprintf('\n');
    fprintf('============================================\n');
    fprintf(' Ossidante: %s\n',nome_ossidante{i});
    fprintf(' Temperatura: %.2f K\n',T_iniziale_ox(i));
    fprintf('============================================\n');

    for j = 1:length(OF_current)

        of = OF_current(j);

        fprintf('O/F = %.2f\n',of);

        try

            %% -------------------------------------------------
            % CEA
            % --------------------------------------------------

            CEA_output = CEA( ...
                'problem','o/f',of, ...
                'rocket','eq', ...
                'p,psia',pc, ...
                'pi/p',pi_p, ...
                'reactants', ...
                    'oxid',ossidante_CEA{i}, ...
                    'wt%',100, ...
                    't,k',T_iniziale_ox(i), ...
                    'fuel','HTPB', ...
                    'wt%',100, ...
                    't,k',T_fuel, ...
                    'h,cal',-2970, ...
                    'C',7.337, ...
                    'H',10.982, ...
                    'O',0.058, ...
                'output', ...
                    'short', ...
                'end');

            %% ---------------------------------------------
            % ESTRAZIONE Isp
            % ---------------------------------------------

            Isp_current(j) = ...
                CEA_output.output.eql.isp(end);

        catch

            fprintf('   --> CEA non converge: punto saltato\n');

            Isp_current(j) = NaN;

        end

    end

    % Salvataggio risultati
    Isp{i} = Isp_current;

end
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

    [Isp_max,index_max] = max(Isp{i},[],'omitnan');

    OF_opt = OF_all{i}(index_max);

    fprintf('%-10s  %-12.2f  %-15.2f\n', ...
        nome_ossidante{i}, ...
        OF_opt, ...
        Isp_max);

end

fprintf('============================================================\n');
%% ============================================================
% CONTROLLO DEL MASSIMO
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('           DETTAGLIO CURVA ATTORNO AL MASSIMO\n');
fprintf('============================================================\n');

for i = 1:N_ossidanti

    fprintf('\n%s\n',nome_ossidante{i});
    fprintf('%-10s %-15s\n','O/F','Isp [s]');
    fprintf('----------------------------\n');

    for j = 1:length(OF_all{i})

        if OF_all{i}(j) >= 1.5 && OF_all{i}(j) <= 3.5

            fprintf('%-10.2f %-15.3f\n', ...
                OF_all{i}(j), ...
                Isp{i}(j));

        end

    end

end


%% ============================================================
% ZOOM SULLA ZONA DEL MASSIMO
% ============================================================

figure

hold on
grid on
box on

for i = 1:N_ossidanti

    plot( ...
        OF_all{i}, ...
        Isp{i}, ...
        'o-', ...
        'LineWidth',2, ...
        'DisplayName',nome_ossidante{i});

end

xlim([1.5 3.5])

xlabel('O/F [-]')
ylabel('I_{sp} [s]')

title('Dettaglio del massimo - HTPB')

legend('show','Location','best')

hold off