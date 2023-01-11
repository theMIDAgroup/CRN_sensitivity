function eval_jac = f_evaluate_jacobian_neworder(k, x, ...
                Sm, idx_basic_species, jac_v, cons_laws)

%% This function permits to evaluate the jacobian of f with respect to x.
%% The result is in the form [S_2*Jacv; N].

%% Step 1. Initialization
n_species = size(Sm, 1);
eval_jac = zeros(n_species);

%% Step 2. Evaluate jacobian of v(x; k)
idx = find(jac_v.Jv ~= 0); 
eval_jac_v = zeros(size(jac_v.Jv)); 
aux_x = [x; 1]; 
eval_jac_v(idx) = aux_x(jac_v.Jv(idx));
eval_jac_v = ( diag(jac_v.factor).*diag(k(jac_v.k)) ) *  eval_jac_v; 

%% Step 3. Evaluate jacobian
Sm_2 = Sm; Sm_2(idx_basic_species, :) = []; 
eval_jac = [Sm_2*eval_jac_v; ...
            cons_laws];



end
