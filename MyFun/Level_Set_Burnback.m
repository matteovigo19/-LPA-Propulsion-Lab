% Level Set / Fast Marching Method per il burnback di un grano ibrido a stella.

% clear; clc; close all;

%% 1. GEOMETRIA E DATI DI COMBUSTIONE (da prevars.mat se presente)

prevars_path = 'prevars.mat';

if isfile(prevars_path)
    d = load(prevars_path);
    n_tips  = double(d.n_tips);
    diam_i  = double(d.diam_i);
    diam_e  = double(d.diam_e);
    Ap0_ref = double(d.Ap);
    Pb0_ref = double(d.Pb);
    L_grain = double(d.L);
    mox     = double(d.mox);
    fprintf('Parametri caricati da ''%s''.\n', prevars_path);
else
    % Valori di fallback (solo per test standalone senza MATLAB a monte)
    % n_tips = 5; diam_i = 0.0500; diam_e = 0.0850;
    % Ap0_ref = NaN; Pb0_ref = NaN; L_grain = 0.30;
    % mox = 0.9;
    warning('''%s'' non trovato: uso valori di default per test standalone.\n', prevars_path);
end

n_tips = 9;

ri = diam_i / 2.0;
re = diam_e / 2.0;

% fuel: HTPB (stessi valori di main1_predesign.m / main2.m)
rho_f = 920.0;
a_rf  = 0.027e-3;    % (m/s)/(kg/m2 s)^n
n_rf  = 0.75;

% raggio della camera (parete esterna del grano)
%ext_diameter = 4*diam_e;
ext_diameter = 300*2e-3;
R_case = ext_diameter/2;

%% 2. GRIGLIA E GEOMETRIA INIZIALE DELLA STELLA

%creo la griglia per l'analisi del burnback

Ngrid = 500;
Lgrid = 1.15 * R_case;
x = linspace(-Lgrid, Lgrid, Ngrid);
y = linspace(-Lgrid, Lgrid, Ngrid);
dx = x(2) - x(1);
[X, Y] = meshgrid(x, y);  % griglia contenuta in vettori

% inserisco i dati geometrici delle stelle
alpha = pi ./ n_tips;
angles = (0:(2*n_tips - 1)) * alpha;
radii = zeros(1, 2*n_tips);
radii(1:2:end) = ri; % Inner radii
radii(2:2:end) = re; % Outer radii

star_x = radii .* cos(angles);
star_y = radii .* sin(angles);

% Chiudiamo il poligono della stella per il plot e la validazione
star_x(end+1) = star_x(1);
star_y(end+1) = star_y(1);

% Equivalente di star_path.contains_points()
inside_star = inpolygon(X, Y, star_x, star_y); % port a t=0

R_grid = hypot(X, Y);
in_case = R_grid <= R_case;            % dominio fisico (fuel + port), fino alla parete
fuel0 = in_case & (~inside_star);      % propellente solido a t=0

cell_area = dx^2;
Ap0_geom = sum(inside_star(in_case)) * cell_area;

fprintf('\nArea porta iniziale (da griglia FMM):  %.4f mm^2\n', Ap0_geom * 1e6);
if ~isnan(Ap0_ref)
    err = abs(Ap0_geom - Ap0_ref) / Ap0_ref;
    fprintf('Area porta iniziale (da predesign):    %.4f mm^2  (errore relativo di discretizzazione: %.2f%%)\n', Ap0_ref * 1e6, err * 100);
end

%% 3. FAST MARCHING: DISTANZA (a velocita' unitaria) DAL BORDO STELLA
% L'equivalente del FMM con v=1 è la Euclidean Distance Transform (bwdist).
% Calcola la distanza minima da ogni punto fino all'area 'inside_star'

T = double(bwdist(inside_star)) * dx; 
T_filled = T;
T_filled(~in_case) = NaN; % fuori camera = fuori dominio, metto NaN

%% 4. TABELLA Ap(s), Pb(s) IN FUNZIONE DELLA DISTANZA DI REGRESSIONE s

s_max = max(T_filled(fuel0)); % distanza massima raggiungibile prima di toccare la parete
n_s = 300;
s_vals = linspace(0, s_max, n_s);

Ap_s = NaN(1, n_s);
Pb_s = NaN(1, n_s);

% Copia di T per i contour in cui rimuoviamo i NaN rimpiazzandoli con un valore irraggiungibile
T_for_contour = T_filled;
T_for_contour(isnan(T_filled)) = 1e9;

for idx = 1:n_s
    s = s_vals(idx);
    
    burned_mask = in_case & (inside_star | (T_filled <= s));
    Ap_s(idx) = sum(burned_mask(:)) * cell_area;

    if s < s_max - 1e-9
        % Trova le curve di livello equivalenti a measure.find_contours
        C = contourc(x, y, T_for_contour, [s s]);
        perim = 0.0;
        
        % Parsing della matrice generata da contourc
        c_idx = 1;
        while c_idx <= size(C, 2)
            num_pts = C(2, c_idx);
            x_pts = C(1, c_idx+1 : c_idx+num_pts);
            y_pts = C(2, c_idx+1 : c_idx+num_pts);
            
            % Calcola la distanza tra i punti del segmento
            seg = sqrt(diff(x_pts).^2 + diff(y_pts).^2);
            perim = perim + sum(seg);
            
            % Passa al contour successivo
            c_idx = c_idx + num_pts + 1;
        end
        Pb_s(idx) = perim;
    else
        Pb_s(idx) = 0.0; % fronte annegato nella parete
    end
end

%% 5. INTEGRAZIONE TEMPORALE: t(s) = INTEGRALE DI ds / rf(s)

Gox_s = mox ./ Ap_s;
rf_s = a_rf .* Gox_s .^ n_rf;

% Maschera per i valori finiti validi
valid = isfinite(rf_s) & (rf_s > 0);

t_vals = zeros(1, n_s);
% cumtrapz(x, y) = integrale di y(x) in dx
t_vals(valid) = cumtrapz(s_vals(valid), 1.0 ./ rf_s(valid));

%% 6. OUTPUT NUMERICO (da confrontare con l'output finale di main2.m)

fprintf('\n--- Risultati Level-Set / FMM ---\n');
fprintf('s_max (regressione fino a parete camera): %.4f mm\n', s_max * 1e3);
fprintf('tempo stimato di web-out:                  %.4f s\n', t_vals(end));
fprintf('Area finale:                               %.4f mm^2\n', Ap_s(end) * 1e6);
fprintf('Perimetro finale:                          %.4f mm\n', Pb_s(end) * 1e3);
fprintf('Gox finale:                                %.4f kg/(m^2 s)\n', Gox_s(end));

%% 7. PLOT

figure('Name', 'FMM Burnback', 'Position', [100, 100, 1000, 800]);
theta = linspace(0, 2*pi, 200);

% 7a. Campo T(x,y)
subplot(2, 2, 1);
hold on;
contourf(X * 1e3, Y * 1e3, T_filled, 30, 'LineStyle', 'none');
colormap(gca);
plot(R_case * cos(theta) * 1e3, R_case * sin(theta) * 1e3, 'r--', 'DisplayName', 'parete camera');
plot(star_x * 1e3, star_y * 1e3, 'w-', 'LineWidth', 1.2, 'DisplayName', 'bordo iniziale');
title('Campo T(x,y): tempo di arrivo a v=1 (FMM)');
xlabel('x [mm]'); ylabel('y [mm]');
axis equal;
legend('Location', 'best', 'FontSize', 8);
colorbar;
hold off;

% 7b. Fronti a diversi s
subplot(2, 2, 2);
hold on;
fracs = linspace ( 0.0, 0.9, 11);
colors = lines(length(fracs)); % palette di colori di default
for i = 1:length(fracs)
    frac = fracs(i);
    s_show = frac * s_max;
    contour(X * 1e3, Y * 1e3, T_filled, [s_show s_show], 'Color', colors(i,:), 'LineWidth', 1.5);
end
plot(R_case * cos(theta) * 1e3, R_case * sin(theta) * 1e3, 'k--');
title('Fronte di combustione a diversi s');
xlabel('x [mm]'); ylabel('y [mm]');
axis equal;
hold off;

% 7c. Area vs tempo
subplot(2, 2, 3);
plot(t_vals(valid), Ap_s(valid) * 1e6, 'b-', 'LineWidth', 1.5);
xlabel('t [s]'); ylabel('Area porta [mm^2]');
title('Area vs tempo (Level Set / FMM)');
grid on;

% 7d. Perimetro vs tempo
subplot(2, 2, 4);
plot(t_vals(valid), Pb_s(valid) * 1e3, 'r-', 'LineWidth', 1.5);
xlabel('t [s]'); ylabel('Perimetro bruciante [mm]');
title('Perimetro vs tempo (Level Set / FMM)');
grid on;

