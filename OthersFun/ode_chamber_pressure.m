function dp_p = ode_chamber_pressure(mdot_in,mdot_out, rho_g, vars, fine_ode)

% ODE_CHAMBER_PRESSURE 
% calcola la variazione nel tempo della pressione in camera normalizzata
% 
% INPUT: 
% mdot_in: portata in massa tot in ingresso alla camera [kg/s]
% mdot_out: portata in massa tot in uscita dalla camera [kg/s]
% rho_g: densità gas [kg/m3]
% geometria
%         grain_lenght: lunghezz del grano [m]
%         port_area: area del porto [m2]
%         d_area_area: variazione area del porto nel tempo normalizzata
%         [1/s]
% fine_ode
%         n_rf: regression rate [-]
%         d_perim_perim: variazione del perimetro nel tempo normalizzata
%         [1/s]
%         deltapc: derivazione nel tempo press camera = p/(RT)*d(RT)/dp
%         deltaOF: derivazione nel tempo O/F = of/(RT)*d(RT)/dof

L = vars.geometry.grain_length;
Ap = vars.geometry.port_area;
d_area_area = vars.geometry.d_area_area;

if nargin<5
    n_rf = 0;
    d_perim_perim = 0;
    deltapc = 0;
    deltaOF = 0;
else 
    n_rf = fine_ode.n_rf;
    d_perim_perim = fine_ode.d_perim_perim;
    deltapc = fine_ode.deltapc;
    deltaOF = fine_ode.deltaOF;
end

term_mass_balance = (mdot_in-mdot_out) / ( L*Ap*rho_g);
term_of_perimeter = deltaOF * d_perim_perim;
term_area_coupling = (1-n_rf*deltaOF)*d_area_area;

% fprintf("\nDEBUG ode_chamber_pressure:\n");
% fprintf("mdot_in          = %.6e kg/s\n", mdot_in);
% fprintf("mdot_out         = %.6e kg/s\n", mdot_out);
% fprintf("mdot_diff        = %.6e kg/s\n", mdot_in - mdot_out);
% fprintf("L                = %.6e m\n", L);
% fprintf("Ap               = %.6e m^2\n", Ap);
% fprintf("rho_g            = %.6e kg/m^3\n", rho_g);
% fprintf("m_gas = L*Ap*rho = %.6e kg\n", L*Ap*rho_g);
% fprintf("term_mass        = %.6e 1/s\n", term_mass_balance);
% fprintf("term_area        = %.6e 1/s\n", term_area_coupling);
% fprintf("term_perim       = %.6e 1/s\n", term_of_perimeter);
% fprintf("deltapc          = %.6e\n", deltapc);
% fprintf("dp_p             = %.6e 1/s\n", ...
%     (1/(1-deltapc)) * (term_mass_balance-term_of_perimeter-term_area_coupling));

dp_p = (1/(1-deltapc)) * (term_mass_balance-term_of_perimeter-term_area_coupling);

end