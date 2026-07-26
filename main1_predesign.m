close all
clear
clc

%% PREDESIGN

thrust = 50000;   % N
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
mdot = thrust/(ue+pe*eps*(cstar/pch));

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
rad_cyl_ref = 0.5*diam_cyl_ref; %raggio del cilindro corrispondente tale che Ap sia coerente
%
%% --------- WHICH GEOMETRY -------------

% caso = "cylinder";

caso = "star";
n_tips = 5;

radius = rad_cyl_ref*0.7;

% caso = "RAT";
% radius = rad_cyl_ref*0.7;

switch caso

    case "star"

        % Angolo tra una punta interna e una punta esterna consecutiva
        alpha = pi/n_tips;

        %% FUNZIONI GEOMETRICHE

        % Area del poligono formato collegando le punte interne.
        % ri è il raggio circoscritto, non l'apotema.
        Area_inner_polygon = @(ri) ...
            0.5*n_tips*ri.^2*sin(2*alpha);

        % Lato del poligono che collega due punte interne consecutive
        side_inner_polygon = @(ri) ...
            2*ri.*sin(alpha);

        % Area di uno dei triangoli esterni.
        %
        % La base è il lato del poligono interno.
        % L'altezza è:
        % re - ri*cos(alpha)
        Area_external_triangle = @(re,ri) ...
            0.5.*side_inner_polygon(ri).* ...
            (re - ri.*cos(alpha));

        % Area complessiva dei triangoli esterni
        Area_external_triangles = @(re,ri) ...
            n_tips.*Area_external_triangle(re,ri);

        % Area complessiva della stella
        Ap_fun = @(re,ri) ...
            Area_inner_polygon(ri) + ...
            Area_external_triangles(re,ri);

        % Formula equivalente semplificata:
        Ap_fun_direct = @(re,ri) ...
            n_tips.*re.*ri.*sin(alpha);

        % Lunghezza del segmento tra una punta interna e una esterna
        star_side_fun = @(re,ri) ...
            sqrt(re.^2 + ri.^2 - ...
            2.*re.*ri.*cos(alpha));

        % Perimetro complessivo della stella
        Pb_fun = @(re,ri) ...
            2*n_tips.*star_side_fun(re,ri);


        %% DETERMINAZIONE DI ri E re

        if radius < rad_cyl_ref

            % radius rappresenta il raggio delle punte interne
            ri = radius;

            % Dall'equazione:
            % Ap = n_tips*ri*re*sin(alpha)
            re = Ap/(n_tips*ri*sin(alpha));

        elseif radius > rad_cyl_ref

            % radius rappresenta il raggio delle punte esterne
            re = radius;

            % Dall'equazione:
            % Ap = n_tips*ri*re*sin(alpha)
            ri = Ap/(n_tips*re*sin(alpha));

        else

            % Caso degenere: tutte le punte sono sulla stessa
            % circonferenza.
            ri = radius;
            re = radius;

        end


        %% CONTROLLI DI COERENZA

        if ri <= 0 || re <= 0
            error('I raggi ri e re devono essere positivi.');
        end

        if re < ri
            error(['Geometria non valida: il raggio delle punte esterne ' ...
                   're risulta minore del raggio interno ri.']);
        end


        %% AREA E PERIMETRO RISULTANTI

        Ap_geometry = Ap_fun(re,ri);
        Pb = Pb_fun(re,ri);

        % Controllo numerico dell'area
        area_relative_error = abs(Ap_geometry - Ap)/max(abs(Ap),eps);

        if area_relative_error > 1e-10
            warning(['L''area geometrica non coincide con Ap. ' ...
                     'Errore relativo: %.3e'], ...
                     area_relative_error);
        end


        %% LUNGHEZZA DEL GRANO

        L = fox/(rho_f*Pb*a_rf*(mox^n_rf)/Ap^n_rf);


        %% DIAMETRI CARATTERISTICI

        diam_i = 2*ri;
        diam_e = 2*re;


        %% STAMPA DI CONTROLLO

        fprintf('\n--- Geometria stella ---\n');
        fprintf('Numero punte:             %d\n', n_tips);
        fprintf('Raggio punte interne ri:  %.6f m\n', ri);
        fprintf('Raggio punte esterne re:  %.6f m\n', re);
        fprintf('Diametro interno:         %.6f m\n', diam_i);
        fprintf('Diametro esterno:         %.6f m\n', diam_e);
        fprintf('Area assegnata Ap:        %.8e m^2\n', Ap);
        fprintf('Area ricalcolata:         %.8e m^2\n', Ap_geometry);
        fprintf('Perimetro Pb:             %.8e m\n', Pb);
        fprintf('Lunghezza grano L:        %.8e m\n', L);

save('prevars.mat', ...
    'L', 'diam_i','diam_e','n_tips','mox','eps','Ap','Pb', ...
    '-mat');

end