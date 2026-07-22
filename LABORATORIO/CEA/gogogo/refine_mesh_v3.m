function P_ref = refine_mesh_v3(P, idx_v_interni)

P_ref = P;

n = size(P,1);

idx_v_interni = idx_v_interni(:);

tol = 1e-12;

for kk = 1:length(idx_v_interni)

    i_tip = idx_v_interni(kk);

    P_tip = P_ref(i_tip,:);
    r_tip = norm(P_tip);

    if r_tip < eps
        continue;
    end

    %% ============================================================
    %  1. CERCO BLOCCO NON VALIDO A SINISTRA DELLA PUNTA
    % ============================================================

    idx_left_invalid = [];

    i = i_tip - 1;

    while i >= 1

        r_i = norm(P_ref(i,:));

        if r_i < r_tip - tol
            idx_left_invalid = [i idx_left_invalid]; %#ok<AGROW>
            i = i - 1;
        else
            break;
        end

    end

    i_left_valid = i;

    %% ============================================================
    %  2. CERCO BLOCCO NON VALIDO A DESTRA DELLA PUNTA
    % ============================================================

    idx_right_invalid = [];

    i = i_tip + 1;

    while i <= n

        r_i = norm(P_ref(i,:));

        if r_i < r_tip - tol
            idx_right_invalid = [idx_right_invalid i]; %#ok<AGROW>
            i = i + 1;
        else
            break;
        end

    end

    i_right_valid = i;

    %% ============================================================
    %  3. CORREGGO SOLO SE HO DUE PUNTI VALIDI DI APPOGGIO
    % ============================================================

    if i_left_valid < 1 || i_right_valid > n
        continue;
    end

    P_left  = P_ref(i_left_valid,:);
    P_right = P_ref(i_right_valid,:);

    %% ============================================================
    %  4. RIPOSIZIONAMENTO LOCALE A SINISTRA
    % ============================================================

    m_left = length(idx_left_invalid);

    for j = 1:m_left

        alpha = j / (m_left + 1);

        P_new = (1 - alpha) * P_left + alpha * P_tip;

        r_new = norm(P_new);

        if r_new < r_tip
            P_new = P_new / max(r_new, eps) * r_tip;
        end

        P_ref(idx_left_invalid(j),:) = P_new;

    end

    %% ============================================================
    %  5. RIPOSIZIONAMENTO LOCALE A DESTRA
    % ============================================================

    m_right = length(idx_right_invalid);

    for j = 1:m_right

        alpha = j / (m_right + 1);

        P_new = (1 - alpha) * P_tip + alpha * P_right;

        r_new = norm(P_new);

        if r_new < r_tip
            P_new = P_new / max(r_new, eps) * r_tip;
        end

        P_ref(idx_right_invalid(j),:) = P_new;

    end

end

%% ============================================================
%  6. CORREZIONE SPECIALE PER ULTIMA PUNTA
% ============================================================
% Questa correzione viene applicata solo se l'ultima punta della mesh
% è una punta interna, cioè quando idx_v_interni contiene n.
%
% Serve per geometrie con numero dispari di punte, quando alcuni punti
% vicino alla chiusura della mezza mesh finiscono sotto l'asse x.
%
% Non elimino punti: sostituisco tutti i punti con y < 0 interpolandoli
% con la stessa ratio usata nella refine standard.

if any(idx_v_interni == n)

    P_ref = refine_last_tip_below_axis(P_ref);

end

end


%% ============================================================
%  FUNZIONE LOCALE: CORREZIONE PUNTI SOTTO ASSE X
% ============================================================

function P_ref = refine_last_tip_below_axis(P_ref)

n = size(P_ref,1);

if n < 3
    return;
end

y = P_ref(:,2);

scale_ref = max(vecnorm(P_ref, 2, 2));

if scale_ref < eps
    scale_ref = 1;
end

tol_axis = 1e-12 * scale_ref;

%% ============================================================
%  1. CONTROLLO SE ESISTONO PUNTI SOTTO L'ASSE X
% ============================================================

idx_below = find(y < -tol_axis);

if isempty(idx_below)
    return;
end

%% ============================================================
%  2. RICOSTRUISCO L'ULTIMA PUNTA SU y = 0
% ============================================================
% Cerco:
%   - l'ultimo punto della mesh con y > 0 prima della zona sotto asse;
%   - il primo punto della mesh con y < 0.
%
% Costruisco il segmento tra questi due punti e trovo
% l'intersezione con y = 0.
%
% Questo nuovo punto viene imposto come ultima punta:
%
%   P_ref(n,:) = [x_axis, 0]

idx_neg = idx_below(1);

idx_pos = find(P_ref(1:idx_neg-1,2) > tol_axis, 1, 'last');

if isempty(idx_pos)

    % Caso degenerato: almeno impongo la simmetria sull'ultima punta.
    P_ref(n,2) = 0;
    return;

end

P_pos = P_ref(idx_pos,:);
P_neg = P_ref(idx_neg,:);

dy_seg = P_neg(2) - P_pos(2);

if abs(dy_seg) < eps

    P_ref(n,2) = 0;
    return;

end

alpha_axis = -P_pos(2) / dy_seg;

alpha_axis = max(0, min(1, alpha_axis));

P_axis = P_pos + alpha_axis * (P_neg - P_pos);

P_ref(n,:) = [P_axis(1), 0];

%% ============================================================
%  3. SOSTITUZIONE DI TUTTI I PUNTI CON y < 0
% ============================================================
% Dopo aver ricostruito l'ultima punta, sostituisco ogni blocco
% consecutivo di punti sotto asse interpolando tra il punto valido
% precedente e il punto valido successivo.
%
% Se il blocco arriva fino alla fine, il punto valido successivo
% è l'ultima punta P_ref(n,:).

y = P_ref(:,2);

idx_bad = y < -tol_axis;

if ~any(idx_bad)
    return;
end

d = diff([false; idx_bad; false]);

idx_start = find(d == 1);
idx_end   = find(d == -1) - 1;

for b = 1:length(idx_start)

    i_start = idx_start(b);
    i_end   = idx_end(b);

    idx_block = i_start:i_end;
    m = length(idx_block);

    %% ------------------------------------------------------------
    % Punto di appoggio sinistro
    %% ------------------------------------------------------------

    if i_start > 1
        P_left = P_ref(i_start - 1,:);
    else
        P_left = P_ref(1,:);
    end

    %% ------------------------------------------------------------
    % Punto di appoggio destro
    %% ------------------------------------------------------------

    if i_end < n
        P_right = P_ref(i_end + 1,:);
    else
        P_right = P_ref(n,:);
    end

    %% ------------------------------------------------------------
    % Riposizionamento con la stessa ratio della refine standard
    %% ------------------------------------------------------------

    for j = 1:m

        alpha = j / (m + 1);

        P_new = (1 - alpha) * P_left + alpha * P_right;

        % Sicurezza: non permetto punti sotto l'asse
        if P_new(2) < 0
            P_new(2) = 0;
        end

        P_ref(idx_block(j),:) = P_new;

    end

end

%% ============================================================
%  4. PULIZIA NUMERICA FINALE
% ============================================================

idx_axis = abs(P_ref(:,2)) <= tol_axis;

P_ref(idx_axis,2) = 0;

P_ref(n,2) = 0;

end