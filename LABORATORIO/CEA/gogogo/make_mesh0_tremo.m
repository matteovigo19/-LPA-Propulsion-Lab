function  [coordmesh,coord_vertici]= make_mesh0_tremo(geometry,n,vars,type)

% il codice estrae le informazioni in coordinate POLARI per poi passarle in
% cartesiane

    switch lower(geometry)
    
        case "cylinder"
            % ALLOCA I DATI DA INPUT FUNZIONE
            diam = vars.diametro; % si può anche togliere, decidere

            theta = linspace(0,pi,n+1)'; % colonna

            % la funzione lavora in radianti su una semicirconferenza
            % quindi qui da 0 fino a pi inserisce pt.i di discretizzazione
            % +1 punti (non capisco il +1 perchè ma va bene) dove
            % theta(1)=0 e theta(n+1)=pi

            r = ones(size(theta)) * diam/2; % riga
            % qui semplicemente ho un vettore della stessa dimensione di
            % theta dove il raggio è ripetuto tutto uguale
            
            [x_mesh,y_mesh]=pol2cart(theta,r);
            coordmesh=[x_mesh,y_mesh];

        case "star"
            % ALLOCA I DATI DA INPUT FUNZIONE
            ntips = vars.ntips/2;
            % /2 per via della semicirconferenza

            theta = linspace(0,pi, 2 * ntips + 1)'; % colonna
            % ntips +1 funziona per geometria. Di questo vettore gli
            % elementi con indice dispari sono le punte e quelli con indice
            % pari sono le valli. Apprezzamento al fatto che, in questa
            % riga crea i theta esatti sia con un numero di punte pari che
            % dispari
            r = ones(size(theta))*vars.diametro_est/2; % colonna
            r(2 : 2 : end ) = vars.diametro_int/2; % colonna
            % crea il vettore con i raggi esterni e poi sostituisci agli
            % elementi di indice pari i raggi interni

            % procede a convertire le coordinate polari dei vertici in CARTESIANE
            [xv,yv] = pol2cart(theta,r); %trovo le coordinate dei vertici in cartesiane
            
            nlati = ntips*2;  % considerando la stella con lati uguali !!!
            n_meshlato = ceil((n-1)/nlati); % quanti punti di discretizzazione su ogni lato?

            x_fill = []; % inizializzazione vettori
            y_fill = [];

            for ii = 1: nlati

                x_segm = linspace(xv(ii),xv(ii+1),n_meshlato + 1)';
                y_segm = linspace(yv(ii),yv(ii+1),n_meshlato + 1)';
                % ad ogni iterazione creano vettori di punti dove
                % x/y_segm(1) rappresenta il vertice i-esimo e
                % x/y_segm(end) rappresenta il vertice (i+1)-esimo.
                % il secondo vertice coincide con il primo dell'iterata
                % successiva.
               
                x_fill = [x_fill; x_segm(1 : end-1)];
                y_fill = [y_fill; y_segm(1 : end-1)]; 
                % salva le coordinate ottenute ad ogni iterazione
                % escludendo l'ultimo di ogni iterata, che coincide con il
                % primo dell'iterata successiva.
            end

            x_mesh = [x_fill;xv(end)];
            y_mesh = [y_fill;yv(end)];
            % Alla fine del ciclo completa il vettore con il punto
            % posizioneto a pi radianti.

            coordmesh = [x_mesh,y_mesh]; % salva tutto in una matrice (n+1)x2
        
        case "flower"
        % ALLOCA I DATI DA INPUT FUNZIONE
            ntips = vars.ntips/2;
            % /2 per via della semicirconferenza

            theta = linspace(0,pi, 2 * ntips + 1)'; % colonna
            % ntips +1 funziona per geometria. Di questo vettore gli
            % elementi con indice dispari sono le punte e quelli con indice
            % pari sono le valli. Apprezzamento al fatto che, in questa
            % riga crea i theta esatti sia con un numero di punte pari che
            % dispari
            r = ones(size(theta))*vars.diametro_est/2; % colonna
            r(2 : 2 : end ) = vars.diametro_int/2; % colonna
            % crea il vettore con i raggi esterni e poi sostituisci agli
            % elementi di indice pari i raggi interni

    otherwise
    error("Unknown geometry type");
    end

    
    if type == "cartesiane"
        [x,y] = pol2cart(theta,r);
        % se voglio il risultato in cartesiane prende le coordinate dei
        % vertici (vale sia per il polare che per il cartesiano) e le
        % trasforma in coordinate CARTESIANE
        coord_vertici = [x,y];
    
    elseif type == "polari"
        coord_vertici = [r,theta];
    else 
        error("Coordinate sconosciute")
    end
    end