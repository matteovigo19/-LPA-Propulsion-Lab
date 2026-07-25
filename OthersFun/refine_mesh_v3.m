function P_ref = refine_mesh_v3(P, idx_v_interni)

P_ref = P;

n = size(P,1);

idx_v_interni = idx_v_interni(:);

tol = 1e-12;

%% ============================================================
%  1. CORREZIONE STANDARD DELLE PUNTE INTERNE
% ============================================================

for kk = 1:length(idx_v_interni)

    i_tip = idx_v_interni(kk);

    P_tip = P_ref(i_tip,:);
    r_tip = norm(P_tip);

    if r_tip < eps
        continue;
    end

    %% ============================================================
    %  CERCA BLOCCO INVALIDO A SINISTRA DELLA PUNTA
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
    %  CERCA BLOCCO INVALIDO A DESTRA DELLA PUNTA
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
    %  SE NON HO DUE PUNTI VALIDI, SALTO
    % ============================================================

    if i_left_valid < 1 || i_right_valid > n
        continue;
    end

    P_left  = P_ref(i_left_valid,:);
    P_right = P_ref(i_right_valid,:);

    %% ============================================================
    %  SOSTITUZIONE BLOCCO SINISTRO
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
    %  SOSTITUZIONE BLOCCO DESTRO
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
%  2. CORREZIONE SEMPRE ATTIVA SU PRIMA E ULTIMA PUNTA
% ============================================================
% Non dipende più da idx_v_interni == n.
% Viene fatta sempre:
%   - vicino alla prima punta
%   - vicino all'ultima punta
%
% Se vicino a una delle due estremità ci sono punti con y < 0,
% questi vengono sostituiti tramite interpolazione tra punti validi.

P_ref = refine_first_tip_below_axis(P_ref);
P_ref = refine_last_tip_below_axis(P_ref);

end


%% ============================================================
%  FUNZIONE LOCALE: CORREZIONE VICINO ALLA PRIMA PUNTA
% ============================================================

function P_ref = refine_first_tip_below_axis(P_ref)

n = size(P_ref,1);

if n < 3
    return;
end

scale_ref = max(vecnorm(P_ref, 2, 2));

if scale_ref < eps
    scale_ref = 1;
end

tol_axis = 1e-12 * scale_ref;

% La prima punta deve stare sull'asse x.
P_ref(1,2) = 0;

y = P_ref(:,2);

idx_bad = y < -tol_axis;

if ~any(idx_bad)
    return;
end

% Considero solo punti negativi vicini alla prima punta.
boundary_window = max(5, ceil(0.15*n));

d = diff([false; idx_bad; false]);

idx_start = find(d == 1);
idx_end   = find(d == -1) - 1;

% Primo blocco negativo vicino all'inizio.
candidate = find(idx_start <= boundary_window, 1, 'first');

if isempty(candidate)
    return;
end

i_start = idx_start(candidate);
i_end   = idx_end(candidate);

%% ============================================================
%  CASO BLOCCO ATTACCATO ALLA PRIMA PUNTA
% ============================================================

if i_start <= 2

    i_right = i_end + 1;

    if i_right > n
        P_ref(1,2) = 0;
        return;
    end

    P_neg = P_ref(i_end,:);
    P_pos = P_ref(i_right,:);

    dy_seg = P_pos(2) - P_neg(2);

    if abs(dy_seg) > eps

        alpha_axis = -P_neg(2) / dy_seg;
        alpha_axis = max(0, min(1, alpha_axis));

        P_axis = P_neg + alpha_axis * (P_pos - P_neg);

        % Aggiorno la prima punta come intersezione con l'asse x.
        P_ref(1,:) = [P_axis(1), 0];

    else

        P_ref(1,2) = 0;

    end

    P_left  = P_ref(1,:);
    P_right = P_ref(i_right,:);

    idx_fix = 2:i_end;

    P_ref = replace_negative_block(P_ref, idx_fix, P_left, P_right);

else

    %% ============================================================
    %  CASO BLOCCO VICINO ALL'INIZIO MA NON ATTACCATO ALLA PUNTA
    % ============================================================

    i_left  = i_start - 1;
    i_right = i_end + 1;

    if i_left < 1 || i_right > n
        return;
    end

    P_left  = P_ref(i_left,:);
    P_right = P_ref(i_right,:);

    idx_fix = i_start:i_end;

    P_ref = replace_negative_block(P_ref, idx_fix, P_left, P_right);

end

% Sicurezza finale.
P_ref(1,2) = 0;

end


%% ============================================================
%  FUNZIONE LOCALE: CORREZIONE VICINO ALL'ULTIMA PUNTA
% ============================================================

function P_ref = refine_last_tip_below_axis(P_ref)

n = size(P_ref,1);

if n < 3
    return;
end

scale_ref = max(vecnorm(P_ref, 2, 2));

if scale_ref < eps
    scale_ref = 1;
end

tol_axis = 1e-12 * scale_ref;

% L'ultima punta deve stare sull'asse x.
P_ref(n,2) = 0;

y = P_ref(:,2);

idx_bad = y < -tol_axis;

if ~any(idx_bad)
    return;
end

% Considero solo punti negativi vicini all'ultima punta.
boundary_window = max(5, ceil(0.15*n));

d = diff([false; idx_bad; false]);

idx_start = find(d == 1);
idx_end   = find(d == -1) - 1;

% Ultimo blocco negativo vicino alla fine.
candidate = find(idx_end >= n - boundary_window + 1, 1, 'last');

if isempty(candidate)
    return;
end

i_start = idx_start(candidate);
i_end   = idx_end(candidate);

%% ============================================================
%  CASO BLOCCO ATTACCATO ALL'ULTIMA PUNTA
% ============================================================

if i_end >= n - 1

    i_left = i_start - 1;

    if i_left < 1
        P_ref(n,2) = 0;
        return;
    end

    P_pos = P_ref(i_left,:);
    P_neg = P_ref(i_start,:);

    dy_seg = P_neg(2) - P_pos(2);

    if abs(dy_seg) > eps

        alpha_axis = -P_pos(2) / dy_seg;
        alpha_axis = max(0, min(1, alpha_axis));

        P_axis = P_pos + alpha_axis * (P_neg - P_pos);

        % Aggiorno l'ultima punta come intersezione con l'asse x.
        P_ref(n,:) = [P_axis(1), 0];

    else

        P_ref(n,2) = 0;

    end

    P_left  = P_ref(i_left,:);
    P_right = P_ref(n,:);

    idx_fix = i_start:n-1;

    P_ref = replace_negative_block(P_ref, idx_fix, P_left, P_right);

else

    %% ============================================================
    %  CASO BLOCCO VICINO ALLA FINE MA NON ATTACCATO ALLA PUNTA
    % ============================================================

    i_left  = i_start - 1;
    i_right = i_end + 1;

    if i_left < 1 || i_right > n
        return;
    end

    P_left  = P_ref(i_left,:);
    P_right = P_ref(i_right,:);

    idx_fix = i_start:i_end;

    P_ref = replace_negative_block(P_ref, idx_fix, P_left, P_right);

end

% Sicurezza finale.
P_ref(n,2) = 0;

end


%% ============================================================
%  FUNZIONE LOCALE: SOSTITUZIONE BLOCCO NEGATIVO
% ============================================================

function P_ref = replace_negative_block(P_ref, idx_fix, P_left, P_right)

if isempty(idx_fix)
    return;
end

m = length(idx_fix);

for j = 1:m

    alpha = j / (m + 1);

    P_new = (1 - alpha) * P_left + alpha * P_right;

    % Sicurezza: nessun punto sotto l'asse x.
    if P_new(2) < 0
        P_new(2) = 0;
    end

    P_ref(idx_fix(j),:) = P_new;

end

end