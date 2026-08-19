function [dY, T] = ode_mesh_v4(~, Y, idx_v_interni, rf, mesh_mode)

if nargin < 5 || isempty(mesh_mode)
    mesh_mode = "star";
end

mesh_mode = string(mesh_mode);

Y = Y(:);
n = length(Y)/2;

x = Y(1:n);
y = Y(n+1:end);

P = [x y];

idx_v_interni = idx_v_interni(:);

%% ============================================================
%  0. CHECK BASE
% ============================================================

if n < 2
    dY = zeros(size(Y));
    T = zeros(n,2);
    return;
end

%% ============================================================
%  1. VELOCITÀ NORMALE PER TUTTI I PUNTI
% ============================================================
% Tutti i punti, comprese le punte esterne, regrediscono normalmente.
%
% Per i bordi uso punti ghost speculari rispetto all'asse x.
% Questo è importante perché primo e ultimo punto della mezza geometria
% devono restare sull'asse x.

x_ext = [x(2); x; x(end-1)];
y_ext = [-y(2); y; -y(end-1)];

tx = x_ext(3:end) - x_ext(1:end-2);
ty = y_ext(3:end) - y_ext(1:end-2);

mod_t = hypot(tx, ty);
mod_t(mod_t == 0) = eps;

tx = tx ./ mod_t;
ty = ty ./ mod_t;

T = [tx ty];

% Normale candidata
nx = ty;
ny = -tx;

%% ============================================================
%  2. SCELTA DEL VERSO DELLA NORMALE
% ============================================================

r_all = hypot(x, y);

erx = zeros(n,1);
ery = zeros(n,1);

valid = r_all > eps;

erx(valid) = x(valid) ./ r_all(valid);
ery(valid) = y(valid) ./ r_all(valid);

verso = nx .* erx + ny .* ery;

flip = verso < 0;

nx(flip) = -nx(flip);
ny(flip) = -ny(flip);

% Velocità normale base per tutti i punti
dx = rf * nx;
dy = rf * ny;

% Sicurezza: primo e ultimo punto restano sull'asse x
dy(1) = 0;
dy(n) = 0;

%% ============================================================
%  3. EVENTUALE CORREZIONE SOLO DELLE PUNTE INTERNE
% ============================================================
% Questa correzione viene usata solo quando la geometria è ancora stellata.
%
% Differenza rispetto alla versione iniziale:
% se una punta interna coincide con i = 1 oppure i = n,
% non viene saltata. Viene trattata usando un ghost point speculare.

if mesh_mode == "star"

    for k = 1:length(idx_v_interni)

        i = idx_v_interni(k);

        if i < 1 || i > n
            continue;
        end

        P_tip = P(i,:);

        %% ------------------------------------------------------------
        % Costruzione dei punti adiacenti alla punta
        %% ------------------------------------------------------------

        if i == 1

            % Prima punta sul bordo di simmetria.
            % Il lato mancante viene ricostruito con il ghost speculare.
            P_left  = [P(2,1), -P(2,2)];
            P_right = P(2,:);

        elseif i == n

            % Ultima punta sul bordo di simmetria.
            % Il lato mancante viene ricostruito con il ghost speculare.
            P_left  = P(n-1,:);
            P_right = [P(n-1,1), -P(n-1,2)];

        else

            % Punta interna standard
            i_left  = i - 1;
            i_right = i + 1;

            if i_left < 1 || i_right > n
                continue;
            end

            P_left  = P(i_left,:);
            P_right = P(i_right,:);

        end

        %% ------------------------------------------------------------
        % Tangenti dei due lati che formano la punta interna
        %% ------------------------------------------------------------

        t_left  = P_tip - P_left;
        t_right = P_right - P_tip;

        nt_left  = norm(t_left);
        nt_right = norm(t_right);

        if nt_left < eps || nt_right < eps
            continue;
        end

        t_left  = t_left  / nt_left;
        t_right = t_right / nt_right;

        %% ------------------------------------------------------------
        % Normali ai due lati
        %% ------------------------------------------------------------

        n_left  = [ t_left(2),  -t_left(1)];
        n_right = [ t_right(2), -t_right(1)];

        %% ------------------------------------------------------------
        % Verso corretto delle normali
        %% ------------------------------------------------------------

        r_tip = norm(P_tip);

        if r_tip < eps
            continue;
        end

        e_r = P_tip / r_tip;

        if dot(n_left, e_r) < 0
            n_left = -n_left;
        end

        if dot(n_right, e_r) < 0
            n_right = -n_right;
        end

        %% ------------------------------------------------------------
        % Velocità della punta interna come intersezione dei lati regrediti
        %% ------------------------------------------------------------

        A = [n_left;
             n_right];

        b = rf * [1; 1];

        if abs(det(A)) > 1e-10

            v_tip = A \ b;

        else

            n_mean = n_left + n_right;

            if norm(n_mean) > eps
                v_tip = rf * n_mean / norm(n_mean);
            else
                v_tip = rf * e_r;
            end

        end

        % Sovrascrivo SOLO la punta interna
        dx(i) = v_tip(1);
        dy(i) = v_tip(2);

        % Se la punta corretta è sul bordo della mezza geometria,
        % deve rimanere sull'asse x.
        if i == 1 || i == n
            dy(i) = 0;
        end

    end

end

%% ============================================================
%  4. SICUREZZA FINALE SUI PUNTI DI SIMMETRIA
% ============================================================

dy(1) = 0;
dy(n) = 0;

%% ============================================================
%  5. OUTPUT
% ============================================================

dY = [dx; dy];

end