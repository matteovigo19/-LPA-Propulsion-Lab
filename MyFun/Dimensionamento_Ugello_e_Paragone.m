clc;
close all;
clear;

%% Variabili da Predesign
load("prevars.mat");  % Estrazione At e k da predesign
eps = 200;            % Rapporto di espansione desiderato
pch_bar = 20;         % Pressione in camera [bar]
Rt = sqrt(At/pi);
l_percent = 80;       % Percentuale di lunghezza standard (es. 80%)

%% Calcolo dimensioni e plot


bell_nozzle_main(At, eps, k, l_percent); % funzione principale 
[angles, contour] = bell_nozzle_func(k, eps, Rt, 80); 
compare_pressure_rao_vs_cone(contour, Rt, eps, pch_bar, k)


%% =========================================================================
%%  FUNZIONE MAIN per costruzione ugello con RAO
% =========================================================================
function bell_nozzle_main(At, eps, k, l_percent)
%
% GIUSTO DARE I CREDITI, non sono un puzzone
%
%  Implemented from the following technical notes 
%  The thrust optimised parabolic nozzle
%  http://www.aspirespace.org.uk/downloads/Thrust%20optimised%20parabolic%20nozzle.pdf
%  .................................................................
% 
% The radius of the nozzle exit: 
% Re = √ε * Rt							[Eqn. 2]
% and nozzle length 
% LN = 0.8 ((√∈−1) * Rt )/ tan(15)		[Eqn. 3]
% .................................................................
% For the throat entrant section:
% x = 1.5 Rt cosθ
% y = 1.5 Rt sinθ + 1.5 Rt + Rt			[Eqn. 4]
% where: −135 ≤ θ ≤ −90
% (The initial angle isn't defined and is up to the
% combustion chamber designer, -135 degrees is typical.)
% .................................................................
% For the throat exit section:
% x = 0.382 Rt cosθ
% y = 0.382 Rt sinθ + 0.382 Rt + Rt		[Eqn. 5]
% where: −90 ≤ θ ≤ (θn − 90)
% .................................................................
% The bell is a quadratic Bézier curve, which has equations:
% x(t) = (1 − t)^2 * Nx + 2(1 − t)t * Qx + t^2 * Ex, 0≤t≤1
% y(t) = (1 − t)^2 * Ny + 2(1 − t)t * Qy + t^2 * Ey, 0≤t≤1 [Eqn. 6]
% .................................................................
% Selecting equally spaced divisions between 0 and 1 produces 
% the points described earlier in the graphical method, 
% for example 0.25, 0.5, and 0.75.
% .................................................................
% Equations 6 are defined by points N, Q, and E (see the graphical method 
% earlier for the locations of these points).
% 
% Point N is defined by equations 5 setting the angle to (θn – 90).
% Nx = 0.382 Rt cos(θn – 90)
% Ny = 0.382 Rt sin(θn – 90) + 0.382 Rt + Rt
% .................................................................
% Coordinate Ex is defined by equation 3, and coordinate Ey is defined by equation 2.
% Ex = 0.8*(((√ε−1)-1)*Rt)/(tan(15)) # degrees in rad
% Ey = √ε * Rt
% .................................................................
% Point Q is the intersection of the lines: ⃗⃗⃗⃗⃗⃗
% NQ = m1 x + C1 and: ⃗⃗⃗⃗⃗
% QE = m2 x + C2 			[Eqn. 7]
% 
% where: gradient 
% m1 = tan(θn ) , m2 = tan(θe )	[Eqn. 8]
% 
% and: intercept 
% C1 = Ny − m1 Nx
% C2 = Ey − m2 Ex		[Eqn. 9]
% .................................................................
% The intersection of these two lines (at point Q) is given by:
% Qx = (C2 − C1 ) /(m1 − m2 )
% Qy = (m1 C2 − m2 C1 ) / (m1 − m2 ) [Eqn. 10]
% .................................................................
% 

    aratio = eps;     % Ae / At = 200
    throat_radius = sqrt(At/pi);   % [m]
    
    [angles, contour] = bell_nozzle_func(k, aratio, throat_radius, l_percent);
    title_str = sprintf('Bell Nozzle \n [Area Ratio = %.1f, Throat Area = %.3f m, Rt = %.3f m]', aratio, At, throat_radius);
    plot_nozzle_and_3d(title_str, throat_radius, angles, contour);
end

% =========================================================================
%  FUNZIONE DI CALCOLO GEOMETRIA UGELLO
% =========================================================================
function [angles, contour] = bell_nozzle_func(k, aratio, Rt, l_percent)
    entrant_angle = -135;
    ea_radian = deg2rad(entrant_angle);
    
    if l_percent == 60
        Lnp = 0.6;
    elseif l_percent == 80
        Lnp = 0.8;
    elseif l_percent == 90
        Lnp = 0.9;
    else
        Lnp = 0.8;
    end
    
    % Trova gli angoli di parete con gestione estrapolazione sicura per eps=200
    [nozzle_length, theta_n, theta_e] = find_wall_angles(aratio, Rt, l_percent);
    angles = [nozzle_length, theta_n, theta_e];
    data_interval = 150;
    
    % 1. Sezione Entrante (Throat Entrant)
    angle_list = linspace(ea_radian, -pi/2, data_interval);
    xe = 1.5 * Rt * cos(angle_list);
    ye = 1.5 * Rt * sin(angle_list) + 2.5 * Rt;
    
    % 2. Sezione di Uscita dalla Gola (Throat Exit)
    angle_list = linspace(-pi/2, theta_n - pi/2, data_interval);
    xe2 = 0.382 * Rt * cos(angle_list);
    ye2 = 0.382 * Rt * sin(angle_list) + 1.382 * Rt;
    
    % 3. Sezione a Campana (Curva di Bézier Quadratica)
    Nx = 0.382 * Rt * cos(theta_n - pi/2);
    Ny = 0.382 * Rt * sin(theta_n - pi/2) + 1.382 * Rt;
    
    Ex = Lnp * ((sqrt(aratio) - 1) * Rt) / tan(deg2rad(15));
    Ey = sqrt(aratio) * Rt;
    m1 = tan(theta_n);
    m2 = tan(theta_e);
    
    C1 = Ny - m1 * Nx;
    C2 = Ey - m2 * Ex;
    
    Qx = (C2 - C1) / (m1 - m2);
    Qy = (m1 * C2 - m2 * C1) / (m1 - m2);
    
    int_list = linspace(0, 1, data_interval);
    xbell = zeros(size(int_list));
    ybell = zeros(size(int_list));
    
    for idx = 1:length(int_list)
        t = int_list(idx);
        xbell(idx) = ((1-t)^2)*Nx + 2*(1-t)*t*Qx + (t^2)*Ex;
        ybell(idx) = ((1-t)^2)*Ny + 2*(1-t)*t*Qy + (t^2)*Ey;
    end
    
    % Specchiatura simmetrica
    nye  = -ye;
    nye2 = -ye2;
    nybell = -ybell;
    contour = {xe, ye, nye, xe2, ye2, nye2, xbell, ybell, nybell};
end

% =========================================================================
%  RICERCA ANGOLI DI PARETE ED INTERPOLAZIONE / ESTRAPOLAZIONE SICURA
% =========================================================================
function [Ln, tn_rad, te_rad] = find_wall_angles(ar, Rt, l_percent)
    aratio      = [4,      5,    10,    20,    30,    40,    50,   100];
    theta_n_60  = [26.5, 28.0,  32.0,  35.0,  36.2,  37.1,  35.0,  40.0];
    theta_n_80  = [21.5, 23.0,  26.3,  28.8,  30.0,  31.0,  31.5,  33.5];
    theta_n_90  = [20.0, 21.0,  24.0,  27.0,  28.5,  29.5,  30.2,  32.0];
    theta_e_60  = [20.5, 20.5,  16.0,  14.5,  14.0,  13.5,  13.0,  11.2];
    theta_e_80  = [14.0, 13.0,  11.0,   9.0,   8.5,   8.0,   7.5,   7.0];
    theta_e_90  = [11.5, 10.5,   8.0,   7.0,   6.5,   6.0,   6.0,   6.0];
    
    f1 = ((sqrt(ar) - 1) * Rt) / tan(deg2rad(15));
    
    if l_percent == 60
        theta_n = theta_n_60; theta_e = theta_e_60; Ln = 0.6 * f1;
    elseif l_percent == 80
        theta_n = theta_n_80; theta_e = theta_e_80; Ln = 0.8 * f1;
    elseif l_percent == 90
        theta_n = theta_n_90; theta_e = theta_e_90; Ln = 0.9 * f1;
    else
        theta_n = theta_n_80; theta_e = theta_e_80; Ln = 0.8 * f1;
    end
    
    if ar <= max(aratio)
        [x_index, ~] = find_nearest(aratio, ar);
        if round(aratio(x_index), 1) == round(ar, 1)
            tn_val = theta_n(x_index);
            te_val = theta_e(x_index);
        else
            if (x_index > 2) && (x_index < length(aratio) - 1)
                idx_range = x_index-1 : x_index+2;
            elseif x_index <= 2
                idx_range = 1:4;
            else
                idx_range = length(aratio)-3 : length(aratio);
            end
            tn_val = interp1(aratio(idx_range), theta_n(idx_range), ar, 'linear');
            te_val = interp1(aratio(idx_range), theta_e(idx_range), ar, 'linear');
        end
    else
        % Estrapolazione per ar > 100 (es. ar = 200) con clamp di sicurezza per te
        tn_val = interp1(aratio(end-1:end), theta_n(end-1:end), ar, 'linear', 'extrap');
        te_val = interp1(aratio(end-1:end), theta_e(end-1:end), ar, 'linear', 'extrap');
        te_val = max(2.5, te_val); % Clamp di sicurezza per evitare angoli negativi
    end
    
    tn_rad = deg2rad(tn_val);
    te_rad = deg2rad(te_val);
end

function [idx, val] = find_nearest(array, value)
    [~, idx] = min(abs(array - value));
    val = array(idx);
end

% =========================================================================
%  FUNZIONI DI PLOT CON TUTTE LE QUOTE RICHIESTE
% =========================================================================
function plot_nozzle_and_3d(title_str, Rt, angles, contour)
    figure('Position', [50, 50, 1500, 700]);
    
    % Subplot 1: Profilo 2D dettagliato con tutte le quote
    ax1 = subplot(1, 2, 1);
    plot_nozzle_2d(ax1, title_str, Rt, angles, contour);
    
    % Subplot 2: Modello 3D
    ax2 = subplot(1, 2, 2, 'Projection', 'perspective');
    plot3D(ax2, contour);
    
    sgtitle(title_str, 'FontSize', 12, 'FontWeight', 'bold');
end

function plot_nozzle_2d(ax, title_str, Rt, angles, contour)
    nozzle_length = angles(1); theta_n = angles(2); theta_e = angles(3);
    
    xe = contour{1};   ye = contour{2};   nye = contour{3};
    xe2 = contour{4};  ye2 = contour{5};  nye2 = contour{6};
    xbell = contour{7}; ybell = contour{8}; nybell = contour{9};
    
    hold(ax, 'on');
    ax.DataAspectRatio = [1 1 1];
    
    % Plot sezioni ugello
    plot(ax, xe, ye, 'LineWidth', 2.0, 'Color', 'g');
    plot(ax, xe, nye, 'LineWidth', 2.0, 'Color', 'g');
    plot(ax, xe2, ye2, 'LineWidth', 2.0, 'Color', 'r');
    plot(ax, xe2, nye2, 'LineWidth', 2.0, 'Color', 'r');
    plot(ax, xbell, ybell, 'LineWidth', 2.0, 'Color', 'b');
    plot(ax, xbell, nybell, 'LineWidth', 2.0, 'Color', 'b');
    
    % Assi e griglia
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    grid(ax, 'on');
    grid(ax, 'minor');
    
    % --- AGGIUNTA DI TUTTE LE QUOTE E ANNOTAZIONI RICHIESTE ---
    
    % 1. Raggio di gola (Rt)
    plot(ax, [0, xe(end)], [0, ye(end)], 'k--', 'LineWidth', 0.8);
    text(ax, xe(end)/2, ye(end)/2 + Rt*0.1, sprintf(' R_t = %.2f m', Rt), 'FontSize', 9, 'Color', 'r', 'FontWeight', 'bold');
    
    % 2. Raggio in ingresso (Ri) e Lunghezza in ingresso (Li)
    plot(ax, [xe(1), xe(1)], [0, nye(1)], 'k:', 'LineWidth', 0.8);
    text(ax, xe(1)*1.05, nye(1)/2, sprintf(' R_i = %.2f m', abs(nye(1))), 'FontSize', 9, 'FontWeight', 'bold');
    text(ax, xe(1)/2, -Rt*0.15, sprintf(' L_i = %.2f m', abs(xe(1))), 'FontSize', 9, 'FontWeight', 'bold');
    
    % 3. Raccordi geometrici (1.5 * Rt e 0.382 * Rt)
    text(ax, xe(round(end/2)), ye(round(end/2)) + Rt*0.2, sprintf(' 1.5 \\cdot R_t = %.2f', 1.5*Rt), 'FontSize', 8, 'Color', [0 0.5 0]);
    text(ax, xe2(round(end/2)), ye2(round(end/2)) + Rt*0.2, sprintf(' 0.382 \\cdot R_t = %.2f', 0.382*Rt), 'FontSize', 8, 'Color', 'r');
    
    % 4. Raggio di uscita (Re) e Lunghezza totale (Ln)
    plot(ax, [xbell(end), xbell(end)], [0, ybell(end)], 'k:', 'LineWidth', 0.8);
    text(ax, xbell(end)*1.02, ybell(end)/2, sprintf(' R_e = %.2f m', ybell(end)), 'FontSize', 9, 'FontWeight', 'bold');
    text(ax, xbell(end)/2, -Rt*0.15, sprintf(' L_n = %.2f m', nozzle_length), 'FontSize', 9, 'FontWeight', 'bold');
    
    % 5. Angoli theta_n e theta_e
    text(ax, xe2(end)*1.1, nye2(end)*1.1, sprintf(' \\theta_n = %.1f°', rad2deg(theta_n)), 'FontSize', 9, 'FontWeight', 'bold', 'Color', 'm');
    text(ax, xbell(end)*0.9, ybell(end)*0.9, sprintf(' \\theta_e = %.1f°', rad2deg(theta_e)), 'FontSize', 9, 'FontWeight', 'bold', 'Color', 'm');
    
    title(ax, title_str, 'FontSize', 10);
    xlabel(ax, 'Asse Assiale x [m]');
    ylabel(ax, 'Asse Radiale y [m]');
    hold(ax, 'off');
end

function plot3D(ax, contour)
    xe = contour{1};   ye = contour{2};   nye = contour{3};
    xe2 = contour{4};  ye2 = contour{5};  nye2 = contour{6};
    xbell = contour{7}; ybell = contour{8}; nybell = contour{9};
    
    x_all = [xe(:); xe2(:); xbell(:)];
    y_all = [ye(:); ye2(:); ybell(:)];
    
    thick = 5 * (x_all(2) - x_all(1));
    
    hold(ax, 'on');
    for i = 1:length(y_all)
        [X, Y, Z] = ring(y_all(i), thick, x_all(i));
        surf(ax, X, Y, Z, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', '1');
    end
    
    axis(ax, 'equal');
    grid(ax, 'on');
    view(ax, -170, -15);
    colormap(ax, winter);
    title(ax, 'Modello 3D Estruso', 'FontSize', 10);
    hold(ax, 'off');
end

function [X, Y, Z] = ring(r, h, a)
    n_theta = 36;
    n_height = 5;
    theta = linspace(0, 2*pi, n_theta);
    v = linspace(a, a + h, n_height);
    [theta, v] = meshgrid(theta, v);
    X = r * cos(theta);
    Y = r * sin(theta);
    Z = v;
end

% =========================================================================
%% Funzione di paragone
% =========================================================================
function compare_pressure_rao_vs_cone(contour, Rt, eps, pch_bar, k)

    % COMPARE_PRESSURE_RAO_VS_CONE Confronta geometria e caduta di pressione 
    % tra l'ugello a campana (Rao) e un ugello a cono retto (15°).
    %
    % INPUT:
    % contour   -> Cell array con le coordinate geometriche della campana
    % Rt        -> Raggio di gola [m]
    % eps       -> Rapporto di espansione (Ae / At) [-]
    % pch_bar   -> Pressione in camera [bar]
    % k         -> Coefficiente di espansione isentropica [-]

    pch = pch_bar * 1e5; % [Pa]
    At = pi * Rt^2;
    R_e = sqrt(eps) * Rt;

    %% 1. ESTRAZIONE PROFILO SUPERSONICO RAO DALLA CONTOUR
    xe2 = contour{4};  ye2 = contour{5};   % Sezione d'uscita gola
    xbell = contour{7}; ybell = contour{8}; % Sezione campana di Bézier
    
    % Uniamo le coordinate per avere l'intero tratto supersonico dalla gola all'uscita
    x_rao = [xe2(:)', xbell(:)'];
    r_rao = [ye2(:)', ybell(:)'];
    
    % Ordiniamo i punti lungo l'asse x per sicurezza
    [x_rao, sort_idx] = sort(x_rao);
    r_rao = r_rao(sort_idx);

    %% 2. GENERAZIONE UGELLO A CONO RETTO (15 gradi)
    alpha_cone = deg2rad(15);
    L_cone = (R_e - Rt) / tan(alpha_cone);
    
    % Creiamo un vettore x con lo stesso numero di punti per un confronto pulito
    x_cone = linspace(0, L_cone, length(x_rao));
    r_cone = Rt + (R_e - Rt) * (x_cone / L_cone);

    %% 3. CALCOLO PRESSIONE STATICA ISENTROPICA (Rao vs Conico)
    % Funzione implicita per legare il rapporto d'area locale al numero di Mach
    eps_fun = @(M) (1./M) .* (((2./(k+1)) .* (1 + 0.5*(k-1).*(M.^2))) .^ ((k+1)/(2*(k-1))));  

    % Calcolo per la campana di Rao
    p_rao = zeros(size(x_rao));
    Me_vect_Rao = zeros(size(x_rao));
    for i = 1:length(x_rao)
        eps_local = (pi * r_rao(i)^2) / At;
        if abs(eps_local - 1.0) < 1e-4
            Me = 1.0;
            Me_vect_Rao(i) = Me;
        else
            Me = fzero(@(M) eps_fun(M) - eps_local, [1.0, 30]);
            Me_vect_Rao(i) = Me;
        end
        p_rao(i) = pch / ((1 + 0.5*(k-1)*(Me^2))^(k/(k-1)));
    end

    % Calcolo per il cono retto
    p_cone = zeros(size(x_cone));
    Me_vect_cone = zeros(size(x_cone));
    for i = 1:length(x_cone)
        eps_local = (pi * r_cone(i)^2) / At;
        if abs(eps_local - 1.0) < 1e-4
            Me = 1.0;
            Me_vect_cone(i) = Me;
        else
            Me = fzero(@(M) eps_fun(M) - eps_local, [1.0, 30]);
            Me_vect_cone(i) = Me;
        end
        p_cone(i) = pch / ((1 + 0.5*(k-1)*(Me^2))^(k/(k-1)));
    end

    %% 4. PLOT DEI CONFRONTI (Geometria e Pressione)
    figure('Name', 'Confronto Rao vs Conico: Geometria, Pressioni e Mach', 'Position', [100, 100, 600, 9000]);
    movegui(gcf, 'center');

    % Subplot 1: Profili Geometrici
    ax1 = subplot(3,1,1);
    hold on; grid on; axis equal;
    plot(x_rao*1e3, r_rao*1e3, 'b-', 'LineWidth', 2, 'DisplayName', 'Campana Rao [L_n = ' + string(round(x_rao(end)*1e3,1)) + ' mm]');
    plot(x_cone*1e3, r_cone*1e3, 'r--', 'LineWidth', 2, 'DisplayName', 'Cono Retto (15°) [L_c = ' + string(round(L_cone*1e3,1)) + ' mm]');
    xlabel('Coordinata Assiale x [mm]');
    ylabel('Raggio r [mm]');
    title('Confronto Ingegneristico dei Profili');
    legend('Location', 'best');

    % Subplot 2: Gradiente di Pressione Statica
    ax2 = subplot(3,1,2);
    hold on; grid on;
    plot(x_rao*1e3, p_rao*1e-5, 'b-', 'LineWidth', 2, 'DisplayName', 'Pressione Rao');
    plot(x_cone*1e3, p_cone*1e-5, 'r--', 'LineWidth', 2, 'DisplayName', 'Pressione Conico');
    xlabel('Coordinata Assiale x [mm]');
    ylabel('Pressione Statica p [bar]');
    title('Confronto Caduta di Pressione Interna');
    legend('Location', 'best');

    ax3 = subplot(3,1,3);
    hold on; grid on;
    plot(x_rao*1e3, Me_vect_Rao,'b-', 'LineWidth', 2,'DisplayName', 'Campana Rao')
    plot(x_cone*1e3, Me_vect_cone,'r--','LineWidth', 2,'DisplayName', 'Cono Retto')
    xlabel('Coordinata Assiale x [mm]')
    ylabel('Mach')
    legend('Location', 'best');

    linkaxes([ax1, ax2, ax3], 'x');
end



% function analyze_bell_angles(At, eps, pch_bar, k)
%     % ANALYZE_BELL_ANGLES Confronta diversi angoli di un ugello a campana
%     %
%     % INPUT:
%     % At      -> Area di gola [m^2]
%     % eps     -> Rapporto di espansione (Ae / At) [-]
%     % pch_bar -> Pressione in camera [bar]
%     % k       -> Coefficiente di espansione isentropica [-]
% 
%     pch = pch_bar * 1e5; % [Pa]
%     r_th = sqrt(At / pi);
%     r_e = sqrt(eps * At / pi);
%     n_points = 150;
% 
%     % Lunghezza di riferimento basata su un cono standard a 15 gradi
%     alpha_cone = deg2rad(15);
%     L_cone = (r_e - r_th) / tan(alpha_cone);
%     L_bell = 0.8 * L_cone; % La campana è generalmente più compatta
% 
%     % Matrice dei casi da testare: [theta_initial (deg), theta_exit (deg)]
%     angle_cases = [
%         30, 5;   % Caso 1: Standard / Ottimizzato tipico (Rao)
%         45, 2;   % Caso 2: Molto aggressivo in gola, parallelo in uscita
%         20, 12   % Caso 3: Molto dolce in gola, ma con alta divergenza in uscita
%     ];
% 
%     figure('Name', 'Analisi Parametrica Angoli Ugello a Campana', 'Position', [100, 100, 1000, 650]);
%     colors = lines(size(angle_cases, 1));
%     eps_fun = @(Me) 1/Me * sqrt(((1 + 0.5*(k-1)*(Me^2)) / (0.5*(k+1)))^((k+1)/(k-1)));
% 
%     for idx = 1:size(angle_cases, 1)
%         theta_i = deg2rad(angle_cases(idx, 1));
%         theta_e = deg2rad(angle_cases(idx, 2));
% 
%         x_bell = linspace(0, L_bell, n_points);
% 
%         % Condizioni al contorno per l'interpolazione cubica della campana
%         c0 = r_th;
%         c1 = tan(theta_i);
% 
%         L = L_bell;
%         A_mat = [L^3, L^2; 3*L^2, 2*L];
%         B_vec = [r_e - r_th - c1*L; tan(theta_e) - c1];
%         coeffs = A_mat \ B_vec;
%         c3 = coeffs(1);
%         c2 = coeffs(2);
% 
%         % Profilo del raggio della campana
%         r_bell = c3 * x_bell.^3 + c2 * x_bell.^2 + c1 * x_bell + c0;
% 
%         % Calcolo della pressione statica interna via isentropica
%         p_bell = zeros(size(x_bell));
%         for i = 1:length(x_bell)
%             eps_local = (pi * r_bell(i)^2) / At;
%             if i == 1
%                 Me = 1.0;
%             else
%                 Me = fzero(@(M) eps_fun(M) - eps_local, [1.0, 15]);
%             end
%             p_bell(i) = pch / ((1 + 0.5*(k-1)*(Me^2))^(k/(k-1)));
%         end
% 
%         %% Subplot 1: Geometria della Campana
%         subplot(2,1,1);
%         hold on; grid on; axis equal;
%         plot(x_bell*1e3, r_bell*1e3, 'LineWidth', 2, 'Color', colors(idx,:), ...
%             'DisplayName', sprintf('\\theta_i = %d°, \\theta_e = %d°', angle_cases(idx,1), angle_cases(idx,2)));
%         plot(x_bell*1e3, -r_bell*1e3, 'LineWidth', 2, 'Color', colors(idx,:), 'HandleVisibility', 'off');
% 
%         %% Subplot 2: Gradiente di Pressione
%         subplot(2,1,2);
%         hold on; grid on;
%         plot(x_bell/L_bell * 100, p_bell*1e-5, 'LineWidth', 2, 'Color', colors(idx,:), ...
%             'DisplayName', sprintf('\\theta_i = %d°, \\theta_e = %d°', angle_cases(idx,1), angle_cases(idx,2)));
%     end
% 
%     % Rifinitura Grafico 1 (Geometria)
%     subplot(2,1,1);
%     xline(0, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
%     xlabel('Ascissa assiale x [mm]');
%     ylabel('Raggio r [mm]');
%     title('Confronto dei Profili Geometrici al variare degli angoli');
%     legend('Location', 'best');
% 
%     % Rifinitura Grafico 2 (Pressione)
%     subplot(2,1,2);
%     xlabel('Posizione assiale normalizzata [% della lunghezza ugello]');
%     ylabel('Pressione Statica p_c [bar]');
%     title('Confronto dei Gradienti di Pressione Interna');
%     legend('Location', 'best');
% end