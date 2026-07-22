close all
clear
clc

%% PREDESIGN
% dati prefissati a istante 0 
trust = 50000;   % N
eps = 200;   % Ae/At
pch_bar = 20;   % bar
pch = pch_bar*1e5;   % Pa
O_F = 2;
GOX = 500;   % 500 to 800 [kg/m2s]

% conviene partire dall'OF più basso possibile
% GOX tende a scendere nel tempo

% FUEL: HTPB
% regression rate: rf = a*Gox^n
rho_f = 920;   % kg/m3
a_rf = 0.027;   % (mm/s)/(kg/m2 s)^n
n_rf = 0.75;   % -
a_rf = a_rf*1e-3;   % m/s
% OXIDIZER: LOX
load("CEA_functions.mat");

%% assuming frozen
% pch, O/F --> chamber
Tch = T_fun_of_p(O_F, pch_bar);
k = k_fun_of_p(O_F, pch_bar);
R = R_fun_of_p(O_F, pch_bar);
K2 = k*((2/(k+1))^((k+1)/(k-1))); % squared Vandenkerckhove
cstar = sqrt(R*Tch/K2);
% eps --> Me
eps_fun = @(Me, k) 1/Me * sqrt(((1+0.5*(k-1)*(Me^2))/(0.5*(k+1)))^((k+1)/(k-1)));
Me = fzero(@(Me) eps-eps_fun(Me, k), 2); %Mach di uscita, parte provando da 2
% Me, pch --> pe
pe = pch/((1+0.5*(k-1)*(Me^2))^(k/(k-1)));
% Me, Tch --> Te
Te = Tch/(1+0.5*(k-1)*(Me^2));
% Te --> ae
ae = sqrt(k*R*Te);
% Me, ae --> ue
ue = Me*ae;
% thrust, ue, pe, eps, cstar, pch --> mdot
% thrust = mdot*ue + pe*eps*At
%       cstar = pch*At/mdot
%       --> At = (cstar/pch)*mdot
%   --> thrust = ue*mdot + pe*eps*(cstar/pch)*mdot
%   --> mdot = T/(ue+pe*eps*(cstar/pch))
mdot = Te/(ue+pe*eps*(cstar/pch));

% cstar, mdot, pch --> At
At = cstar*mdot/pch;
dt = 2*sqrt(At/pi);
% mdot, O/F, --> mox, fox
%   mdot = mox + fox
%       O_F = mox/fox
%       --> mox = O_F*fox
%   --> mdot = O_F*fox + fox
%   --> fox = mdot/(O_F+1)
fox = mdot/(O_F+1);
mox = O_F*fox;
% GOX, mox --> Ap
Ap = mox/GOX;
diam_cyl_ref = 2*sqrt(Ap/pi);
rad_cyl_ref = 0.5*diam_cyl_ref;
%
%% ---------WHICH GEOMETRY-------------
% caso = "cylinder";

% caso = "star";
%     n = 5;
%     radius = rad_cyl_ref*0.7;

caso = "RAT";
    radius = rad_cyl_ref*0.7;

switch caso
    case "star"
Area_inner_polygon = @(ri) n*(ri^2)*tan(pi/n);
% Perim_inner_polygon = 2*n*ri*tan(pi/n);
side_inner_polygon = @(ri) 2*ri*tan(pi/n);
Area_external_triangle = @(re, ri) 0.25*(re-ri)*side_inner_polygon(ri);
Area_external_triangles = @(re, ri) n*Area_external_triangle(re,ri);
Ap_fun = @(re, ri) Area_inner_polygon(ri) + Area_external_triangles(re,ri);

if radius < rad_cyl_ref
    %radius is the inner radius
    ri = radius;
    re = fzero(@(re) Ap-Ap_fun(re, ri), rad_cyl_ref);
elseif radius > rad_cyl_ref
    %radius is the outer radius
    re = radius;
    ri = fzero(@(ri) Ap-Ap_fun(re, ri), rad_cyl_ref);
else
    % this is a cylinder...
    ri = radius;
    re = radius;
end
l_triangle = sqrt((re-ri)^2+(0.5*side_inner_polygon(ri))^2);
Pb = 2*n*l_triangle;
L = fox/(rho_f*Pb*a_rf*(mox^n_rf)/Ap^n_rf);
diam_in = ri*2;
diam_e = re*2;

    case "RAT"
% ------- rod and tube --------
        %INPUT
        radius = rad_cyl_ref*0.7;
        %~~

        area_in = @(ri) pi*(ri^2);
        area_e = @(re) pi*(re^2);
        Ap_fun = @(re, ri) area_e(re) - area_in(ri);
        if radius <= rad_cyl_ref
            ri = radius;
            re = fzero(@(re) Ap-Ap_fun(re, ri), rad_cyl_ref);
        elseif radius > rad_cyl_ref
            re = radius;
            ri = fzero(@(ri) Ap-Ap_fun(re, ri), rad_cyl_ref);
        end
Pb = 2*pi*(re+ri);
L = fox/(rho_f*Pb*a_rf*(mox^n_rf)/Ap^n_rf);
diam_in = ri*2;
diam_e = re*2;

case "cylinder"
% ----- cylinder -----
    % fox = rho_f*Ab*fr
    % --> fox = rho_f*L*Pb*a_rf*(mox^n_rf)/Ap^n_rf
    % --> L = fox/(rho_f*Pb*a_rf*(mox^n_rf)/Ap^n_rf)
    diam = diam_cyl_ref;
    Pb = pi*diam;
    L = fox/(rho_f*Pb*a_rf*(mox^n_rf)/Ap^n_rf);
otherwise
        error("case not defined")
end