%% TODO per velocizzare.
% Io calcolo tutto il sistema S*v e poi sostituisco le leggi di cons.
% Potrei già calcolare il sistema ridotto S_2*v.

function f_x = f_evaluate_mim(rate_constants, x, idx_basic_species, ...
                                Nl, rho, S, v, ind_one)

aux_x = x;
aux_x(ind_one) = 1;
coeff_eval = rate_constants(v(:,1));
species1_eval = aux_x(v(:,2));
species2_eval = aux_x(v(:,3));

f_x = S * (coeff_eval.*species1_eval.*species2_eval);
f_x(idx_basic_species) = Nl*x - rho;

end