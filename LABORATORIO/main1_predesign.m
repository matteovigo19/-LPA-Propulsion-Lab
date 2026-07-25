close all
clear
clc

%% PREDESIGN
% dati prefissati a istante 0 

thrust = 50000;   % N

% rapporto di espansione dell'ugello
eps = 200;   % Ae/At

% pressione in camera
pch_bar = 20;   % bar
pch = pch_bar*1e5;   % Pa

% rapporto ossidante/combustibile
O_F = 2;

% oxidizer mass flux (flusso massico di ossidante)
GOX = 500;   % 500 to 800 [kg/m2s]

% conviene partire dall'OF più basso possibile
% GOX tende a scendere nel tempo

% FUEL: HTPB
% regression rate: rf = a*Gox^n
rho_f = 920;   % kg/m3 --> densità combustibile
a_rf = 0.027;   % (mm/s)/(kg/m2 s)^n 
n_rf = 0.75;   % -
a_rf = a_rf*1e-3;   % m/s
% OXIDIZER: LOX
load("CEA_functions.mat"); % carica gli output ottenuti dal main_CEA

%% assuming frozen
% pch, O/F --> chamber

% dall'output del CEA:
Tch = T_fun_of_p(O_F, pch_bar);
k = k_fun_of_p(O_F, pch_bar);
R = R_fun_of_p(O_F, pch_bar);

K2 = k*((2/(k+1))^((k+1)/(k-1))); % squared Vandenkerckhove
cstar = sqrt(R*Tch/K2); % Velocità caratteristica (c*)

% eps --> Me (Mach in uscita dall'ugello)
eps_fun = @(Me, k) 1/Me * sqrt(((1+0.5*(k-1)*(Me^2))/(0.5*(k+1)))^((k+1)/(k-1)));
Me = fzero(@(Me) eps-eps_fun(Me, k), 2); % Mach di uscita, parte provando da 2

% Me, pch --> pe (Pressione in uscita dall'ugello)
pe = pch/((1+0.5*(k-1)*(Me^2))^(k/(k-1)));

% Me, Tch --> Te (Temperatura in uscita dall'ugello)
Te = Tch/(1+0.5*(k-1)*(Me^2));

% Te --> ae (Velocità sonica in uscita dall'ugello)
ae = sqrt(k*R*Te);

% Me, ae --> ue (Velocità in uscita dall'ugello)
ue = Me*ae;

% thrust, ue, pe, eps, cstar, pch --> mdot
% thrust = mdot*ue + pe*eps*At --> invertendo/sostituendo si trova mdot
%       cstar = pch*At/mdot
%       --> At = (cstar/pch)*mdot
%   --> thrust = ue*mdot + pe*eps*(cstar/pch)*mdot
%   --> mdot = T/(ue+pe*eps*(cstar/pch))
mdot = thrust/(ue+pe*eps*(cstar/pch)); % mdot_OX+mdot_F

% cstar, mdot, pch --> At (Area di gola)
At = cstar*mdot/pch; 
dt = 2*sqrt(At/pi); % diametro di gola

% mdot, O/F, --> mox, fox
%   mdot = mox + fox
%       O_F = mox/fox
%       --> mox = O_F*fox
%   --> mdot = O_F*fox + fox
%   --> fox = mdot/(O_F+1)
fox = mdot/(O_F+1); % fuel mass flow rate (mdot_F)
mox = O_F*fox; % oxidizer mass flow rate (mdot_OX)

% GOX, mox --> Ap (Area di porto)
Ap = mox/GOX;
diam_cyl_ref = 2*sqrt(Ap/pi);  % diametro cilindro equivalente con area Ap
rad_cyl_ref = 0.5*diam_cyl_ref;  % raggio cilindro equivalente con area Ap

% fissata l'area del porto si procede a svolgere i calcoli relativi alla
% geometria scelta
 
%% ---------WHICH GEOMETRY-------------

% SISTEMARE I COMMENTI IN BASE ALLA SCELTA DELLA GEOMETRIA

% caso = "cylinder";

caso = "star";
    n = 5; % numero di punte della stella
    radius = rad_cyl_ref*0.7; % il valore del fattore non è fisso (?)

%caso = "RAT"; % stella o RAT
%    radius = rad_cyl_ref*0.7; % il valore del fattore non è fisso (?)

%% PARTE CHE OPERA E SVOLGE I CALCOLI IN BASE ALLO SWITCH
switch caso
    case "star"
Area_inner_polygon = @(ri) n*(ri^2)*tan(pi/n); % Area poligono regolare interno, senza punte
% Perim_inner_polygon = 2*n*ri*tan(pi/n);
side_inner_polygon = @(ri) 2*ri*tan(pi/n); % Lato poligono regolare
Area_external_triangle = @(re, ri) 0.5*(re-ri)*side_inner_polygon(ri); % Area di una punta della stella
Area_external_triangles = @(re, ri) n*Area_external_triangle(re,ri); % Area di tutte le punte
Ap_fun = @(re, ri) Area_inner_polygon(ri) + Area_external_triangles(re,ri); % Area totale stella

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

% calcoli sulle dimensioni del grano
l_triangle = sqrt((re-ri)^2+(0.5*side_inner_polygon(ri))^2); % lunghezza lato punta
Pb = 2*n*l_triangle; % lunghezza del perimetro su cui avviene la combustione
L = fox/(rho_f*Pb*a_rf*(mox^n_rf)/Ap^n_rf); % lunghezza assiale del grano (?)

% diametro interno/esterno finale (su cui va costruita la mesh)
diam_in = ri*2;
diam_e = re*2;

    case "RAT"
% ------- rod and tube --------

area_in = @(ri) pi*(ri^2); % Area rod
area_e = @(re) pi*(re^2); % Area tube
Ap_fun = @(re, ri) area_e(re) - area_in(ri); % Area del porto (=corona circolare)
if radius <= rad_cyl_ref
    ri = radius;
    re = fzero(@(re) Ap-Ap_fun(re, ri), rad_cyl_ref);
elseif radius > rad_cyl_ref
    re = radius;
    ri = fzero(@(ri) Ap-Ap_fun(re, ri), rad_cyl_ref);
end

Pb = 2*pi*(re+ri); % lunghezza del perimetro su cui avviene la combustione
L = fox/(rho_f*Pb*a_rf*(mox^n_rf)/Ap^n_rf); % lunghezza assiale del grano (?)

% diametro rod/tube finale (su cui va costruita la mesh)
diam_in = ri*2;
diam_e = re*2;

case "cylinder"
 % fox = rho_f*Ab*fr
 % --> fox = rho_f*L*Pb*a_rf*(mox^n_rf)/Ap^n_rf
 % --> L = fox/(rho_f*Pb*a_rf*(mox^n_rf)/Ap^n_

 diam = diam_cyl_ref; % direttamente quello del porto
 Pb = pi*diam; % lunghezza del perimetro su cui avviene la combustione
 L = fox/(rho_f*Pb*a_rf*(mox^n_rf)/Ap^n_rf); % lunghezza assiale del grano (?)

otherwise
    error("case not defined")
end