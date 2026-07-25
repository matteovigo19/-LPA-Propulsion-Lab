function dY_out = freeze_points_on_camera(Y, dY, D_camera)

Y  = Y(:);
dY = dY(:);

n = length(Y)/2;

x = Y(1:n);
y = Y(n+1:end);

dx = dY(1:n);
dy = dY(n+1:end);

R_camera = D_camera / 2;

r = hypot(x, y);

tol_wall = 1e-12 * R_camera;

idx_wall = r >= R_camera - tol_wall;

dx(idx_wall) = 0;
dy(idx_wall) = 0;

dY_out = [dx; dy];

end