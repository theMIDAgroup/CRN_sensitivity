function jacobian = f_evaluate_jacobian_c(n_species, n_idx)

%% Compute the Jacobian of f with respect to vector c

% Case gof mutation

I_n_idx=eye(n_idx);
jacobian=[zeros(n_species-n_idx, n_idx); -I_n_idx];

end
