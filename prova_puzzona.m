%% Verifica interpolazione Pb_s(t)

Pb_s_mm = Pb_s * 1e3;  % conversione m -> mm, coerente con perimetro_t_mm

% Controllo di monotonia richiesto da interp1
assert(all(diff(t_vals(valid)) > 0), 'ATTENZIONE: t_vals non è strettamente crescente');

Pb_s_interp_mm = interp1(t_vals(valid), Pb_s_mm(valid), t_ode, 'linear', 'extrap');

figure('Name', 'Verifica interpolazione perimetro', 'Position', [100 100 900 600]);

% 1. Curva "sorgente" (level-set, campionata su s) + punti di interpolazione richiesti
subplot(2,1,1); hold on;
plot(t_vals(valid), Pb_s_mm(valid), 'b.-', 'MarkerSize', 8, 'DisplayName', 'Pb\_s originale (nodi level-set)');
plot(t_ode, Pb_s_interp_mm, 'ro', 'MarkerSize', 5, 'DisplayName', 'punti interpolati a t\_ode');
plot(t_ode, perimetro_t_mm, 'g.', 'MarkerSize', 10, 'DisplayName', 'perimetro Lagrangiano (ODE)');
xlabel('t [s]'); ylabel('Perimetro [mm]');
legend('Location','best');
title('Confronto diretto sulle stesse ascisse temporali');
grid on;

% 2. Densità dei nodi level-set nel tempo (per capire dove l'interpolazione "tira" di più)
subplot(2,1,2);
plot(t_vals(valid), 'b.-');
xlabel('indice campione s'); ylabel('t\_vals [s]');
title('Non-uniformità della griglia temporale del level-set (dt/ds = 1/r_f)');
grid on;

% Errore relativo, ora correttamente allineato nel tempo
err = abs((perimetro_t_mm - Pb_s_interp_mm) ./ Pb_s_interp_mm);
figure();
plot(t_ode, err, 'LineWidth', 1.2);
xlabel('t [s]'); ylabel('errore relativo perimetro');
grid on;