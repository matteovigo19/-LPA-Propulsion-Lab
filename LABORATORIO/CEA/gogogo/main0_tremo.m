clc;
close all;
clear;
%% 1 CREARE LA MESH 

% 1.1 INPUT DATI
geometry = "star"; % scegliere la geometria ("cylinder" o "star")

% SOLO PER GRANO CILINDRICO
vars.diametro = 100;

% SOLO PER GRANO A STELLA
vars.ntips = 8; % numero di punte della stella
vars.diametro_est = 120; % valore del diametro esterno
vars.diametro_int = 60; % valore del diametro interno

% NUMERO DI PUNTI DI DISCRETIZZAZIONE (CIL E STAR)
n = 200; % ATTENZIONE! La funzione calcola il numero di punti su ogni lato, se non è preciso aumenta il numero di punti fino a distribuire lo stesso numero di punti su ogni lato

% CERCA I VERTICI E I PUNTI DI DISCRETIZZAZIONE
[coordmesh, coord_vertici_cart] = make_mesh0_tremo(geometry, n, vars, "cartesiane");
[~,coord_vertici_pol] = make_mesh0_tremo(geometry, n, vars, "polari");

% PLOT CARTESIANE
figure()
plot(coord_vertici_cart(:,1),coord_vertici_cart(:,2), "k-",LineWidth=2); % plotta la figura dai vertici
hold on
plot(coordmesh(:,1),coordmesh(:,2),"r.",MarkerSize=15); % plotta i punti di discretizzazione
axis equal; grid on

%% STEP 2: VALUTARE PERIMETRO E AREA DELLA MESH

% COMPUTE PERIMETER AND AREA
[perimetro, area] = eval_mesh_tremo(coord_vertici_cart, "cartesian");
fprintf("perimetro = %f m\n", perimetro)
fprintf("area = %f m^2\n", area)

%% STEP 3: integrate the mesh
% 
rf = 2;         % Velocità di regressione [m/s]
t_fine = 6;        % Tempo totale di combustione [s] (parte da 0s)

Y0 = [coordmesh(:,1);coordmesh(:,2)]; % preparo l'imput a ode45

% LANCIO DI ode45 PER LA RISOLUZIONE DEL 
[t,sol] = ode45(@(t,y) ode_mesh_v1(t, y, rf), [0 t_fine], Y0);

% RISCOMPOSIZIONE OUTPUT ode45
x = sol(:, 1:0.5*size(sol,2));
y = sol(:, (0.5*size(sol,2) + 1):end);
% vettori x e y rappresentano l'andamento di ogni coordinata ai vari
% istanti di tempo del vettore t

% PLOT SERIO (istante iniziale e finale)
figure()
plot(x(1,:),y(1,:) ,"b.-",LineWidth=2, MarkerSize=15);
hold on
plot(x(end,:),y(end,:) ,"r.-",LineWidth=2, MarkerSize=15);
axis equal;grid on;

% PLOT EVOLUZIONE DELLA MESH (ATTIVA L'INTERRUTTORE)
interruttore = 1;
    if interruttore
    
     figure()
     for k = 1:size(x,1)
    
         cla % cancella l'iterazione precedente
    
         plot(x(k,:),y(k,:) ,"k.-",LineWidth=2, MarkerSize=15);
         axis equal;grid on;
         drawnow;
         pause(0.1)
     end
    end
%% PROSSIMA SECTION