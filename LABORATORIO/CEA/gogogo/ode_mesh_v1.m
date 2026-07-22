function dY = ode_mesh_v1(~, Y, rf)

% ~~~~~~~~~~~~~~~~~~~~~~~~~ %
% NAME: ode_mesh_v1(~, Y, rf)
% DESCRIPTION: create differential mesh in time dY/dt
%
% INPUT:
%
%   Y  --> coordinate MESH cartesiane come vettore perchè la ode
%   richiede il vettoree
%          Y = [x1 ... xn, y1 ... yn] (2nx1) [m]
%
%   rf --> regression rate [m/s]
%
% OUTPUT:
%
%   dY --> differential mesh [m]
%          dY = dY/dt
% ~~~~~~~~~~~~~~~~~~~~~~~~~ %


%

Y = Y(:);
n = 0.5*length(Y);

%separo le coordinate dal vettore riga 

x = Y(1:n);
y = Y(n+1:end);

%estendo il vettore coordinate per poter applicare le differenze finite

x = [x(2); x; x(end-1)];         % sono vettori colonna ; impila
y = [-y(2); y; -y(end-1)];

% creo 2 vettori con le coordinate dei vettori tangenti tra i punti che si
% susseguono

tx = x(3:end) - x(1:end-2);
ty = y(3:end) - y(1:end - 2);
mod = hypot(tx,ty);  %calcolo il modulo per normalizzare
tx = tx./mod;
ty = ty./mod;
nx = ty; % invece che usare la matrice di rotazione basta cambiare i segn correttamente 
ny = -tx;

norY = [nx,ny];  %affianco le coordinate dei versori normali, diventa n x 2

dY = norY.*rf;   % dY/dt trovo lo spostamento infinitesimo normale alla superficie
%rf può poi essere definito anche come funzione per questo non impilo prima
dY = dY(:); %impilo nx su ny
end