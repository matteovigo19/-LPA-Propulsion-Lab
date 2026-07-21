function [P_ref, idx_wall] = refine_mesh_camera(P, D_camera)

% REFINE_MESH_CAMERA
%
% Proietta sulla parete della camera tutti i punti che superano
% il raggio della camera.
%
% INPUT:
%   P        = mesh [n x 2]
%   D_camera = diametro camera [m]
%
% OUTPUT:
%   P_ref    = mesh corretta
%   idx_wall = punti che hanno raggiunto/superato la camera

P_ref = P;

R_camera = D_camera / 2;

x = P(:,1);
y = P(:,2);

r = hypot(x, y);

tol_wall = 1e-12 * R_camera;

idx_wall = r >= R_camera - tol_wall;

for i = 1:size(P,1)

    if idx_wall(i)

        if r(i) > eps

            % Proiezione radiale sulla camera
            P_ref(i,:) = R_camera * P(i,:) / r(i);

        end

    end

end

end