function [x_nodi] = nodi_CGL(n_nodi, I)
%nodi di Chebyshev-Gauss-Lobatto

%%input
%n_nodi = punti per l'interpolazione
%I = [a,b] intervallo di riferimento

%%output
%restituisce n + 1 nodi

a = I(1);
b = I(2);
n = n_nodi;
x_nodi = [];

for ii = 0 : n

    x_i = - cos(pi/n*ii);
    x_nodi(ii + 1) = (a + b)/2 + (b - a)/2*x_i;


end
end