%  CONVENZIONE LEVEL SET:
%    phi < 0  -->  PORTA  (gas, gia' bruciato)
%    phi > 0  -->  SOLIDO (propellente non bruciato)
%
%  EQUAZIONE RISOLTA (Osher-Sethian, fronte in espansione):
%    d(phi)/dt + F*|grad(phi)| = 0 ,   F = r_F > 0
%
%  IPOTESI CHIAVE (da dichiarare nell'elaborato):
%   1) r_F dipende solo da G_ox "bulk" (m_ox / A_port totale), NON da un
%      valore locale del flusso -> F e' spazialmente uniforme ad ogni
%      istante. Questo e' cio' che rende trattabile il problema con
%      un semplice level-set a velocita' costante (in spazio, non in tempo).
%   2) Combustione 2D estrusa: si trascura il burning sulle facce di
%      testa del grano (end-burning). A_burn = perimetro(phi=0) * Lp.
%   3) m_ox costante nel tempo (come richiesto dal testo, throttling escluso).
%   4) Il merging di porte multiple e' gestito automaticamente dal level
%      set (nessuna logica ad hoc necessaria).
%   5) Nessun accoppiamento con la chimica (CEA) in questo script: qui si
%      calcola solo la GEOMETRIA (A_port(t), A_burn(t)) e da questa G_ox(t)
%      e O/F(t). L'accoppiamento con Tc, c*, Isp va fatto a valle,
%      leggendo Gox(t) e O/F(t) da qui e interpolando i risultati CEA.
%
%  REQUISITI: Image Processing Toolbox (bwdist) per la reinizializzazione.
%             Se non disponibile, vedi nota alla funzione reinit_phi().
% =========================================================================

clear; clc; close all;

%% ---------------- DATI DAL TESTO D'ESAME ----------------
a_reg   = 0.027;      % (mm/s)/(kg/m2/s)^n  - costante regressione HTPB
n_reg   = 0.75;       % [-]                 - esponente regressione
rho_f   = 920;        % kg/m3               - densita' HTPB

t_burn  = 300;        % s  - durata combustione (fissa, da testo)

%% ---------------- OUTPUT DEL TUO DIMENSIONAMENTO (PUNTO i) ----------------
%  !!! SOSTITUISCI QUESTI VALORI CON QUELLI DEL TUO PREDIMENSIONAMENTO !!!
%  Qui sono solo valori PLACEHOLDER plausibili per far girare lo script.
mox     = 10;          % kg/s  - portata ossidante (COSTANTE, tolleranza al grammo)
Lp      = 1.20;         % m     - lunghezza del grano
R_case  = 0.10;         % m     - raggio interno del casing/liner (limite fisico di burnout)

geometry = 'multi';     % 'single' | 'multi' | 'star'

switch geometry
    case 'single'
        % Porta singola centrale (caso "standard" del punto i)
        R0 = sqrt(mox/(pi*500));   % dimensionata da Gox0 = 500 kg/m2s (esempio da testo)
        ports = struct('xc',0,'yc',0,'R',R0);

    case 'multi'
        % N porte circolari disposte su una circonferenza (caso punto iii)
        Nports  = 7;
        R0      = 0.014;    % m - raggio iniziale di ciascuna porta
        Rpitch  = 0.045;    % m - raggio di disposizione dei centri
        ang = linspace(0, 2*pi, Nports+1); ang(end) = [];
        ports = struct('xc', num2cell(Rpitch*cos(ang)), ...
                        'yc', num2cell(Rpitch*sin(ang)), ...
                        'R',  R0);

    case 'star'
        % Gestito direttamente nella sezione INITIAL LEVEL SET sotto
end

%% ---------------- GRIGLIA DI CALCOLO ----------------
Ng     = 300;                 % risoluzione griglia (aumenta per accuratezza)
Lgrid  = 1.05*R_case;
xv = linspace(-Lgrid, Lgrid, Ng);
yv = linspace(-Lgrid, Lgrid, Ng);
dx = xv(2)-xv(1); dy = yv(2)-yv(1);
[X,Y] = meshgrid(xv, yv);
A_cell = dx*dy;

%% ---------------- LEVEL SET INIZIALE (signed distance) ----------------
switch geometry
    case {'single','multi'}
        phi = inf(size(X));
        for k = 1:numel(ports)
            d = sqrt((X-ports(k).xc).^2 + (Y-ports(k).yc).^2) - ports(k).R;
            phi = min(phi, d);      % unione di cerchi = minimo delle distanze
        end

    case 'star'
        Npunte = 7; Rin = 0.020; Rout = 0.050;
        th = atan2(Y,X); rr = sqrt(X.^2 + Y.^2);
        starR = Rin + (Rout-Rin)*0.5*(1+cos(Npunte*th));   % profilo a stella illustrativo
        phi = rr - starR;
end

%% ---------------- TIME MARCHING ----------------
t = 0; it = 0;
reinit_every = 15;     % iterazioni tra una reinizializzazione e l'altra

t_hist = []; Aport_hist = []; Ab_hist = []; Gox_hist = []; rf_hist = []; OF_hist = [];

A_case_tot = pi*R_case^2;

while t < t_burn

    % ---- misure geometriche allo stato corrente ----
    portMask = phi < 0;
    A_port = sum(portMask(:)) * A_cell;

    % perimetro di combustione via marching-squares del contorno phi=0
    Cmat = contourc(xv, yv, phi, [0 0]);
    Ab_perimeter = 0; idxp = 1;
    while idxp < size(Cmat,2)
        npts = Cmat(2,idxp);
        pts  = Cmat(:, idxp+1 : idxp+npts);
        Ab_perimeter = Ab_perimeter + sum(sqrt(sum(diff(pts,1,2).^2,1)));
        idxp = idxp + npts + 1;
    end
    A_burn = Ab_perimeter * Lp;    % area di combustione (assunzione: estrusione uniforme)

    % ---- accoppiamento balistico ----
    Gox   = mox / A_port;                  % kg/(m2 s)
    rf_mm = a_reg * Gox^n_reg;             % mm/s
    rf    = rf_mm * 1e-3;                  % m/s
    mf    = rho_f * rf * A_burn;           % kg/s
    OF    = mox / mf;

    % ---- salvataggio storico ----
    t_hist(end+1)      = t;      %#ok<SAGROW>
    Aport_hist(end+1)  = A_port; %#ok<SAGROW>
    Ab_hist(end+1)     = A_burn; %#ok<SAGROW>
    Gox_hist(end+1)    = Gox;    %#ok<SAGROW>
    rf_hist(end+1)     = rf;     %#ok<SAGROW>
    OF_hist(end+1)     = OF;     %#ok<SAGROW>

    % ---- controllo di burn-through (porta arriva al casing) ----
    if A_port > 0.97*A_case_tot
        warning('Grano quasi completamente consumato a t = %.2f s. Arresto anticipato.', t);
        break;
    end

    % ---- passo temporale (condizione CFL) ----
    dt = 0.4 * min(dx,dy) / rf;
    dt = min(dt, t_burn - t);
    if dt <= 0, break; end

    % ---- differenze upwind con padding "replicate" (no wrap-around) ----
    phiP = [phi(:,1), phi, phi(:,end)];
    phiP = [phiP(1,:); phiP; phiP(end,:)];

    Dxm = (phiP(2:end-1,2:end-1) - phiP(2:end-1,1:end-2)) / dx;
    Dxp = (phiP(2:end-1,3:end)   - phiP(2:end-1,2:end-1)) / dx;
    Dym = (phiP(2:end-1,2:end-1) - phiP(1:end-2,2:end-1)) / dy;
    Dyp = (phiP(3:end,2:end-1)   - phiP(2:end-1,2:end-1)) / dy;

    gradmag2 = max(Dxm,0).^2 + min(Dxp,0).^2 + max(Dym,0).^2 + min(Dyp,0).^2;

    % ---- avanzamento del fronte: phi_t + F|grad phi| = 0 ----
    phi = phi + dt * rf .* sqrt(gradmag2);

    % ---- reinizializzazione periodica (phi torna a essere signed distance) ----
    it = it + 1;
    if mod(it, reinit_every) == 0
        phi = reinit_phi(phi, dx);
    end

    t = t + dt;
end

%% ---------------- POST-PROCESSING ----------------
figure('Name','Grain Burnback - Level Set');

subplot(2,2,1);
plot(t_hist, Aport_hist*1e4, 'LineWidth', 1.5); hold on;
plot(t_hist, Ab_hist*1e4/max(Ab_hist)*max(Aport_hist*1e4)*0 + Ab_hist*1e4,'--','LineWidth',1.5); %#ok
plot(t_hist, Ab_hist*1e4, '--', 'LineWidth', 1.5);
legend('A_{port} [cm^2]','A_{burn} [cm^2]','Location','best');
xlabel('t [s]'); ylabel('Area [cm^2]'); grid on; title('Geometria');

subplot(2,2,2);
plot(t_hist, Gox_hist, 'LineWidth', 1.5);
xlabel('t [s]'); ylabel('G_{ox} [kg/m^2 s]'); grid on; title('Flusso di massa ossidante');

subplot(2,2,3);
plot(t_hist, rf_hist*1e3, 'LineWidth', 1.5);
xlabel('t [s]'); ylabel('r_F [mm/s]'); grid on; title('Tasso di regressione');

subplot(2,2,4);
plot(t_hist, OF_hist, 'LineWidth', 1.5);
xlabel('t [s]'); ylabel('O/F [-]'); grid on; title('O/F shift');

% visualizzazione finale della geometria di porta
figure('Name','Sezione finale del grano');
contourf(xv, yv, phi, [0 0], 'LineWidth', 1.5); axis equal tight;
hold on;
theta = linspace(0,2*pi,200);
plot(R_case*cos(theta), R_case*sin(theta), 'k--', 'LineWidth', 1.2);
title(sprintf('Sezione a t = %.1f s (tratteggio = casing)', t_hist(end)));
xlabel('x [m]'); ylabel('y [m]');

fprintf('\nImpulso totale approssimato (F costante, solo per riferimento geometrico):\n');
fprintf('Durata simulata: %.1f s su %.0f s richiesti\n', t_hist(end), t_burn);


%% ========================================================================
%  FUNZIONE DI REINIZIALIZZAZIONE
% =========================================================================
function phi_r = reinit_phi(phi, dx)
% Ricostruisce phi come distanza euclidea con segno rispetto al fronte
% phi=0, usando bwdist (Image Processing Toolbox).
%
% ALTERNATIVA SENZA TOOLBOX: risolvere l'equazione di reinizializzazione
% di Sussman (d(phi)/d(tau) = sign(phi0)*(1-|grad phi|)) con qualche decina
% di sotto-iterazioni upwind; piu' lento da implementare ma non richiede
% Image Processing Toolbox. Utile se lavori su un PC senza licenza completa.

    mask = phi < 0;                     % true = porta
    dist_out = bwdist(mask)  * dx;       % distanza dal solido alla porta piu' vicina
    dist_in  = bwdist(~mask) * dx;       % distanza dalla porta al solido piu' vicino
    phi_r = dist_out - dist_in;          % phi<0 dentro la porta, per convenzione
end