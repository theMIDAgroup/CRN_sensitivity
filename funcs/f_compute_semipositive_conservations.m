function cons_laws = f_compute_semipositive_conservations(stoich)

% For a given stoichiometric matrix S, this function returns a set of 
% generator of the convex polyhedral cone defining all possible 
% semi-positive conservation vectors.
% The function implements the algorithm described in [1].
%
% References
% -----------
% [1] Schister & Hofer. 1991 J.Chem. Soc. Faraday Trans.
%     Determining all extreme semi-positive conservation relations in
%     chemical reaction systems: a test criterion for conservativity.

% Input:
% stoich : array of double n x r 
%   Stoichiometric matric

% Output : array of int p x n 
%   Generator of the convex cone of semi-positive conservation vectors.

% (n=number of species, r=number of reactions, p=number of generators)

%% Step 1. Initialize the tableaux T^(0)
[n_sp, n_reac] = size(stoich);
tab = [stoich, eye(n_sp)];

%% Step 2. Iteration T^(j) -> T^(j+1)
for jj = 1:n_reac
    
    column_j = tab(:, jj);
   
    % Step 2.1. Find raws of T^(j) that already have zero on column j+1.
    aux_tab = tab(column_j==0, :);

    % Step 2.2. Define vectors theta.
    neg_coeff = find(column_j < 0);
    pos_coeff = find(column_j > 0);
    I_matr = tab(:, n_reac+1:end);
    
    for ii = neg_coeff'
        for kk = pos_coeff'
        
            % 2.2.a. Compute the intersection of S(ii) and S(kk)
            I_matr = tab(:, n_reac+1:end);
            idx_column = find(sum(I_matr([ii, kk], :), 1)==0);
            I_matr([ii, kk], :) = [];
            % 2.2.b. Check wheter such an intersection belong to some S(l)
            check = find(sum(I_matr(:, idx_column), 2) == 0);
            if isempty(check)
                new_vec = abs(tab(ii, jj))*tab(kk, :) ...
                        + abs(tab(kk, jj))*tab(ii, :);
                aux_tab = [aux_tab; new_vec]; 
            else
                fprintf('Discarding pair : (%d, %d)\n', ii, kk)
            end
        end
    end
    
    tab = aux_tab; clear aux_tab
end
cons_laws = tab(:, n_reac+1:end);

%% Step 3. Normalize each conservation law by the corresponding cgd
n_cons_laws = size(cons_laws, 1);
for ii = 1:n_cons_laws
    non_null_coeff = cons_laws(ii, cons_laws(ii, :)>0);
    temp_gcd = gcd(non_null_coeff(1), non_null_coeff(1));
    for jj = 2:numel(non_null_coeff)
        temp_gcd = gcd(temp_gcd, non_null_coeff(jj));
    end
    cons_laws(ii, :) = cons_laws(ii, :) / temp_gcd;
end

end

%% Something to try to improve speed
% (1) Posso fermarmi prima se a un certo punto S^(j) è già tutta nulla.
%     Nota S^(j) può cambiare numero di righe, ma il numero di colonne
%     dovrebbe rimanere invariato.
% (2) Simbiology uses sparse matrices.