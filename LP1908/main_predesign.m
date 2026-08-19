close all
clear
clc

%% PREDESIGN
% dati prefissati a istante 0 
% variabili ottimizzabili numericamente: 
%  - thrust iniziale
%  - O/F iniziale
%  - GOX iniziale
%  - geometria
%  - parametri fuel (non ne sono sicuro), magari si possono ricostruire le
%  percentuali a posteriori dopo aver ottimizzato i parametri relativi

% variabili non ottimizzabili numericamente:
%  - swirl injection
%  - cambio di oxidiser 


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

%% ---------WHICH GEOMETRY-------------

% SISTEMARE I COMMENTI IN BASE ALLA SCELTA DELLA GEOMETRIA

 %type = "cylinder";

 type = "star";
    n_tips = 5; % numero di punte della stella
    radius = rad_cyl_ref*0.7; % il valore del fattore non è fisso (?)

% type = "RAT"; % stella o RAT
%    radius = rad_cyl_ref*0.7; % il valore del fattore non è fisso (?)


% fissata l'area del porto si procede a svolgere i calcoli relativi alla
% geometria scelta
 

%% PARTE CHE OPERA E SVOLGE I CALCOLI IN BASE ALLO SWITCH
switch type
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
    'L', 'diam_i','diam_e','n_tips','mox','eps','Ap','Pb','pch_bar','type', ...
    '-mat');

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

save('prevars.mat', ...
    'L', 'diam','mox','eps','Ap','Pb','pch_bar','type', ...
    '-mat');

otherwise
    error("case not defined")
end