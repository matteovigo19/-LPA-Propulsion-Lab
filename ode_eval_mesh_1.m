function [d_perimeter_perimeter, d_area_area] = ode_eval_mesh_1(Y, dY)

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% NAME: ode_eval_mesh(Y, dY)
% DESCRIPTION: calculate differential normalized perimeter and area
%
% INPUT:
%
%   Y --> cartesian! mesh coordinates as line vector
%         Y = [x1 ... xn, y1 ... yn] (2nx1) [m]
%
%   dY --> current velocities [dx1...dxn, dy1...dyn]
%
% OUTPUT:
%
%   d_perimeter_perimeter --> normalized differential perimeter [1/s]
%   d_area_area            --> normalized differential area [1/s]
%
% AUTHOR: Valerio Santolini
% LAST UPDATE: 25/02/2026
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

%%

Y = Y(:);
dY = dY(:);

n = 0.5 * length(Y);

% Extract coordinates
x = Y(1:n);
y = Y(n+1:end);

% Extract velocities (dY/dt)
vx = dY(1:n);
vy = dY(n+1:end);

% Compute differences between adjacent vertices
dx_diff = diff(x);
dy_diff = diff(y);

% Compute relative velocities between adjacent vertices
dvx_diff = diff(vx);
dvy_diff = diff(vy);

% Compute segment lengths (needed for the denominator)
lengths = sqrt(dx_diff.^2 + dy_diff.^2);

current_semiperimeter = sum(lengths);
current_perimeter = 2 * current_semiperimeter;

% Apply the chain rule to each segment
d_lengths = (dx_diff .* dvx_diff + dy_diff .* dvy_diff) ./ lengths;

% Sum and account for symmetry
d_semiperimeter = sum(d_lengths);
d_perimeter = 2 * d_semiperimeter;

%%

% Terms for the trapezoidal rule derivative:
% Sum of heights (y_i + y_{i+1}) and their rates
y_sum = y(1:end-1) + y(2:end);
vy_sum = vy(1:end-1) + vy(2:end);

% Width of segments (x_{i+1} - x_i) and their rates
dx_diff = diff(x);
dvx_diff = diff(vx);

% Apply the product rule to the trapezoidal formula
% d/dt [ 0.5 * y_sum * dx_diff ]
d_semiarea_segments = 0.5 * ...
    (vy_sum .* dx_diff + y_sum .* dvx_diff);

% Sum the segments
d_semiarea = sum(abs(d_semiarea_segments));

% Handle the 'abs' from your original function
% Derivative of abs(u) is sign(u) * du/dt
current_semiarea = abs(trapz(x, y));
current_area = 2 * current_semiarea;

d_semiarea = sign(current_semiarea) * d_semiarea;

% Account for symmetry
d_area = 2 * d_semiarea;

% OUTPUT
d_perimeter_perimeter = d_perimeter / current_perimeter;
d_area_area = d_area / current_area;

end