function eval_jac = f_evaluate_jacobian(k, x, ...
                Sm, idx_basic_species, jac_v, cons_laws)


%% Step 1. Initialization
n_species = size(Sm, 1);
eval_jac = zeros(n_species);

%% Step 2. Evaluate jacobian of v(x; k)
idx = find(jac_v.Jv ~= 0); % troviamo gli indici (su matrice srotolata e quindi vettore)
    % delle entrate della jacobiana che saranno non nulle
eval_jac_v = zeros(size(jac_v.Jv)); % taglia = num reazioni x num variabili
aux_x = [x; 1]; %valori di x in cui calcolare la jacobiana
eval_jac_v(idx) = aux_x(jac_v.Jv(idx)); % nelle entrate della jacobiana che saranno non nulle
    % mettiamo i valori delle x
eval_jac_v = ( diag(jac_v.factor).*diag(k(jac_v.k)) ) *  eval_jac_v; % assemblo:
    % moltiplico i valori delle x già inseriti nella futura jacobiana per i
    % rispettivi fattori (1 o 2) e per i k corrispondenti a ciascuna
    % reazione

%% Step 3. Evaluate jacobian
%   3.a. Basic species
eval_jac(idx_basic_species, :) = cons_laws; % vado a sostituire alcune righe con le righe
    % della matrice N visto che parte della jacobiana è data da N stessa

%   3.b. All other species
aux = setdiff(1:n_species, idx_basic_species);
Sm_2 = Sm; Sm_2(idx_basic_species, :) = []; % creo un S_2 a partire da S,
    % in cui prendo nulle tutte le righe corrispondenti alle specie elementari
eval_jac(aux, :) = Sm_2*eval_jac_v; % moltiplico, come da derivata di S_2 * v(x,k)
    % e inserisco nelle righe corrispondenti alle specie non elementari

end