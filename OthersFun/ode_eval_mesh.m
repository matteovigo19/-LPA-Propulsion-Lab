function [d_perim_perim, d_area_area, dperimeter, darea] = ode_eval_mesh(Y, dY)

%ODE_EVAL_MESH Calcola le derivate normalizzate di perimetro e area.
%
% INPUT:
%   Y  -> mesh corrente:
%         può essere [x y] oppure vettore [x1 ... xn y1 ... yn]
%
%   dY -> velocità della mesh:
%         stesso formato di Y
%
% OUTPUT:
%   d_perim_perim -> (dP/dt)/P    [1/s]
%   d_area_area   -> (dA/dt)/A    [1/s]
%   dperimeter    -> dP/dt
%   darea         -> dA/dt

%% ============================================================
%  CONVERSIONE INPUT IN FORMATO [x y]
% ============================================================

P = to_points(Y);
V = to_points(dY);

if size(P,1) ~= size(V,1)
    error("Y e dY devono avere lo stesso numero di punti.");
end

n = size(P,1);

if n < 3
    error("Servono almeno 3 punti per calcolare area e perimetro.");
end

x  = P(:,1);
y  = P(:,2);

vx = V(:,1);
vy = V(:,2);

%% ============================================================
%  CHIUSURA GEOMETRIA
% ============================================================

% Chiudo la curva aggiungendo il primo punto in fondo.
x_closed  = [x;  x(1)];
y_closed  = [y;  y(1)];

vx_closed = [vx; vx(1)];
vy_closed = [vy; vy(1)];

%% ============================================================
%  PERIMETRO E DERIVATA DEL PERIMETRO
% ============================================================

dx = diff(x_closed);
dy = diff(y_closed);

dvx = diff(vx_closed);
dvy = diff(vy_closed);

edge_length = hypot(dx, dy);

perimeter = sum(edge_length);

dperimeter = 0;

for i = 1:length(edge_length)

    if edge_length(i) > 1e-14

        % Derivata della lunghezza del lato:
        % d/dt ||e|| = e/||e|| dot de/dt
        dperimeter = dperimeter + ...
            (dx(i)*dvx(i) + dy(i)*dvy(i)) / edge_length(i);

    end

end

%% ============================================================
%  AREA E DERIVATA DELL'AREA
% ============================================================

% Area signed con formula di shoelace
area_signed = 0;

for i = 1:n

    ip1 = i + 1;

    area_signed = area_signed + ...
        x_closed(i)*y_closed(ip1) - ...
        y_closed(i)*x_closed(ip1);

end

area_signed = 0.5 * area_signed;

area = abs(area_signed);

% Derivata dell'area signed:
%
% A = 1/2 sum_i cross(P_i, P_{i+1})
%
% dA/dt = 1/2 sum_i [
%           cross(V_i, P_{i+1}) + cross(P_i, V_{i+1})
%        ]

darea_signed = 0;

for i = 1:n

    ip1 = i + 1;

    term1 = vx_closed(i)*y_closed(ip1) - ...
            vy_closed(i)*x_closed(ip1);

    term2 = x_closed(i)*vy_closed(ip1) - ...
            y_closed(i)*vx_closed(ip1);

    darea_signed = darea_signed + term1 + term2;

end

darea_signed = 0.5 * darea_signed;

% Siccome eval_mesh restituisce area positiva, uso la derivata di abs(area_signed)
if area_signed >= 0
    darea = darea_signed;
else
    darea = -darea_signed;
end

%% ============================================================
%  DERIVATE NORMALIZZATE
% ============================================================

if perimeter > 1e-14
    d_perim_perim = dperimeter / perimeter;
else
    d_perim_perim = 0;
end

if area > 1e-14
    d_area_area = darea / area;
else
    d_area_area = 0;
end

end


%% ============================================================
%  FUNZIONE LOCALE: CONVERSIONE STATO -> PUNTI
% ============================================================

function P = to_points(Y)

Y = Y(:);

% Se Y era già una matrice [n x 2], questa conversione con Y(:)
% la rovina. Quindi controllo prima la forma originale.
if ismatrix(Y) && size(Y,2) == 2
    P = Y;
    return
end

n = length(Y)/2;

if abs(n - round(n)) > 0
    error("Il vettore Y deve avere lunghezza pari.");
end

n = round(n);

x = Y(1:n);
y = Y(n+1:end);

P = [x y];

end