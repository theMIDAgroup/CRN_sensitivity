function f_x = f_evaluate_mim(rate_constants, x, idx_basic_species, ...
                                Nl, rho, S, v, ind_one)

%% Step 1. Initialization
aux_x = x;
aux_x(ind_one) = 1;
coeff_eval = rate_constants(v(:,1));
species1_eval = aux_x(v(:,2));
species2_eval = aux_x(v(:,3));

%% Step 2. Evaluate the function f 
f_x = S * (coeff_eval.*species1_eval.*species2_eval);
f_x(idx_basic_species) = Nl*x - rho;

end
