function [Z,properties] = Z_chamber_stst(pc, vars)
%% Extract data
% geometry
Ap = vars.geometry.port_area;
perimB = vars.geometry.burning_perimeter;
L = vars.geometry.grain_length;
At = vars.geometry.throat_area;

% fuel
a_rf = vars.fuel.a_rf;
n_rf = vars.fuel.n_rf;
rho_f = vars.fuel.rho_f;

% combustion
mdot_ox = vars.combustion.mdot_ox;
Tc_fun = vars.combustion.Tc_fun;
R_fun = vars.combustion.R_fun;
k_fun = vars.combustion.k_fun;


%% Calculate
% compute mdotin
Gox = mdot_ox/Ap;
rf = a_rf*(Gox^n_rf);
Ab = perimB*L;
mdot_f = rho_f*Ab*rf;
mdot_in = mdot_f + mdot_ox;

% compute mdotout
O_F = mdot_ox/mdot_f;
pc_bar = pc*1e-5;

Tc = Tc_fun(O_F, pc_bar);
R  = R_fun(O_F, pc_bar);
k  = k_fun(O_F, pc_bar);
K2=k*(2/(k+1))^((k+1)/(k-1));
c_star=sqrt(R*Tc/K2);
mdot_out=pc*At/c_star;

Z= mdot_in - mdot_out;

properties.O_F=O_F;
properties.Gox=Gox;
end