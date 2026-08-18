function [coord, coordmesh, coord_v, idx_v, idx_v_interni, idx_v_esterni] = make_mesh0(geometry,n,vars,type)

%la funzione ha lo scopo di effettuare e caratterizzare una mesh della geometria di grano desiderata
%n è il numeri di lati che voglio abbia la geometria
%npoint sono i punti con cui vado a discretizzare la geometria

% togliendo o mettendo le percentuali passo da 180 a 360 e viceversa

    switch lower(geometry)
    
        case "cylinder"

    diam = vars.geometry.diametro;
    theta = linspace(0, pi, n + 1)';
    r = ones(size(theta)) * diam/2;

    [x_mesh,y_mesh]=pol2cart(theta,r);
    coordmesh=[x_mesh,y_mesh];

    coord = [];
    coord_v = [];
    idx_v = [];
    idx_v_interni = [];
    idx_v_esterni = [];

case "star"

    ntips = vars.geometry.ntips/2;

    theta = linspace(0, pi, 2*ntips + 1)';

    r = ones(size(theta)) * vars.geometry.diametro_est/2;
    r(2:2:end) = vars.geometry.diametro_int/2;

    % Vertici della stella in cartesiane
    [xv, yv] = pol2cart(theta, r);
    coord_v = [xv, yv];

    nlati = ntips*2;
    n_meshlato = ceil(n/nlati);

    x_fill = [];
    y_fill = [];

    % Qui salvo gli indici dei vertici dentro coordmesh
    idx_v = zeros(nlati + 1, 1);

    for ii = 1:nlati

        % Il primo punto del lato ii è un vertice.
        % Siccome x_fill contiene già i punti dei lati precedenti,
        % il nuovo vertice sarà in posizione length(x_fill)+1.
        idx_v(ii) = length(x_fill) + 1;

        Ix = [xv(ii), xv(ii+1)];
        Iy = [yv(ii), yv(ii+1)];

        x_segm = nodi_CGL(n_meshlato + 1, Ix);
        y_segm = nodi_CGL(n_meshlato + 1, Iy);

        % Escludo l'ultimo punto perché coincide col primo del lato dopo
        x_fill = [x_fill; x_segm(1:end-1)'];
        y_fill = [y_fill; y_segm(1:end-1)'];

    end

    % Aggiungo l'ultimo vertice finale
    x_mesh = [x_fill; xv(end)];
    y_mesh = [y_fill; yv(end)];

    coordmesh = [x_mesh, y_mesh];

    % Indice dell'ultimo vertice
    idx_v(end) = length(x_mesh);

    %% Classificazione coerente dei vertici

    % Coordinate dei vertici dentro coordmesh
    xv_mesh = coordmesh(idx_v,1);
    yv_mesh = coordmesh(idx_v,2);

    r_v = hypot(xv_mesh, yv_mesh);

    idx_v_interni = [];
    idx_v_esterni = [];

    for kk = 1:length(idx_v)

        if kk == 1 || kk == length(idx_v)

            % Gli estremi della semi-geometria sono vertici esterni
            idx_v_esterni = [idx_v_esterni; idx_v(kk)];

        else

            if r_v(kk) < r_v(kk-1) && r_v(kk) < r_v(kk+1)

                % Vertice più vicino all'origine rispetto ai vicini:
                % punta interna
                idx_v_interni = [idx_v_interni; idx_v(kk)];

            else

                % Altrimenti lo considero vertice esterno
                idx_v_esterni = [idx_v_esterni; idx_v(kk)];

            end

        end

    end
    %passaggio da polari a cartesiane

    if type == "cartesiane"

        [x,y] = pol2cart(theta,r);
      
        % plot(x,y)
        % hold on
        % 
        % if geometry == "star"
        % 
        % plot(x_mesh,y_mesh,'.')
        % 
        % end

        % axis equal
        % xlabel('x')
        % ylabel('y')
        % grid on
        coord = [x,y];

    
    else

        x = theta;
        y = r;

        % plot(x,y)
        % hold on
        % axis equal
        % xlabel('theta')
        % ylabel('raggio')
        % grid on
        coord = [x,y];

    end

        case "fiore"

    % ============================================================
    %  MEZZA STELLA / FIORE APERTA
    %  NON collega l'ultimo punto al primo
    % ============================================================

    n_punte = vars.geometry.ntips;

    R_ext = vars.geometry.diametro_est / 2;
    R_int = vars.geometry.diametro_int / 2;

    curvatura = 0.25;

    % Vertici della mezza geometria:
    % esterno - interno - esterno - interno - ... - esterno
    n_vertici = 2*n_punte - 1;

    % Lati aperti tra vertici consecutivi
    n_lati = n_vertici - 1;

    if n < n_vertici
        error("Il numero di punti n deve essere almeno pari a 2*ntips - 1.");
    end

    % ============================================================
    %  1. VERTICI DELLA MEZZA GEOMETRIA
    % ============================================================

    theta = linspace(0, pi, n_vertici)';

    r = zeros(n_vertici, 1);
    r(1:2:end) = R_ext;
    r(2:2:end) = R_int;

    x_v = r .* cos(theta);
    y_v = r .* sin(theta);

    coord_v = [x_v y_v];

    idx_v_esterni_geom = (1:2:n_vertici)';
    idx_v_interni_geom = (2:2:n_vertici)';

    % ============================================================
    %  2. DISTRIBUZIONE DEI PUNTI
    %
    %  Voglio ottenere esattamente n punti.
    %  Distribuisco n-1 punti sui lati, poi aggiungo l'ultimo vertice.
    % ============================================================

    n_da_distribuire = n - 1;

    n_base = floor(n_da_distribuire / n_lati);
    resto = n_da_distribuire - n_base * n_lati;

    punti_per_lato = n_base * ones(n_lati, 1);
    punti_per_lato(1:resto) = punti_per_lato(1:resto) + 1;

    coordmesh = [];
    idx_v = zeros(n_vertici, 1);

    count = 0;

    for i = 1:n_lati

        P1 = coord_v(i, :);
        P2 = coord_v(i+1, :);

        idx_v(i) = count + 1;

        m = punti_per_lato(i);

        % Includo P1, escludo P2.
        % P2 sarà il primo punto del lato successivo.
        t = linspace(0, 1, m + 1)';
        t(end) = [];

        % Punto medio
        M = 0.5 * (P1 + P2);

        % Direzione lato
        d = P2 - P1;
        L = norm(d);

        % Normale
        normale = [d(2), -d(1)];
        normale = normale / norm(normale);

        % Normale verso l'esterno
        if dot(normale, M) < 0
            normale = -normale;
        end

        % Punto di controllo Bézier
        C = M + curvatura * L * normale;

        % Bézier quadratica
        B = (1 - t).^2 .* P1 + ...
            2 * (1 - t) .* t .* C + ...
            t.^2 .* P2;

        coordmesh = [coordmesh; B];

        count = size(coordmesh, 1);

    end

    % Aggiungo SOLO l'ultimo vertice finale.
    % Non aggiungo il primo punto alla fine.
    coordmesh = [coordmesh; coord_v(end,:)];

    idx_v(end) = size(coordmesh, 1);

    % ============================================================
    %  3. INDICI DELLE PUNTE
    % ============================================================

    idx_v_esterni = idx_v(idx_v_esterni_geom);
    idx_v_interni = idx_v(idx_v_interni_geom);

    idx_v = idx_v(:);

    % ============================================================
    %  4. OUTPUT
    % ============================================================

    x = coordmesh(:,1);
    y = coordmesh(:,2);

    coord = [x; y];

    end

end