function [perimeter, area] = eval_mesh(Y, type)

% eval_mesh computes perimeter and area of a geometry mesh
%
% INPUT:
% Y     -> mesh coordinates
% type  -> "cartesian" or "polar"
%
% OUTPUT:
% perimeter -> geometry perimeter
% area      -> geometry area

%% Convert coordinates if necessary
switch lower(type)

    case "cartesian"

        x = Y(:,1);
        y = Y(:,2);

    case "polar"

        r = Y(:,1);
        theta = Y(:,2);

        [x, y] = pol2cart(theta, r);

    otherwise

        error("Unknown coordinate type")

end

%% Close the geometry if needed
if x(1) ~= x(end) || y(1) ~= y(end)

    x = [x; x(1)];
    y = [y; y(1)];

end

%% PERIMETER

dx = diff(x);
dy = diff(y);

perimeter = 2*(sum(hypot(dx, dy)) - abs(x(1)-x(end-1)));


%% AREA (shoelace formula)

area_sum = 0;

for i = 1:length(x)-1
    area_sum = area_sum + ...
        (x(i)*y(i+1) - y(i)*x(i+1));
end

area = abs(area_sum);

end