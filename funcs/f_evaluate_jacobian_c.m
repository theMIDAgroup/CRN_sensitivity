function jacobian = f_evaluate_jacobian_c(n_species, n_idx, lof)

%% Compute the Jacobian of f with respect to vector c

if lof==0

% Case gof mutation

I_n_idx=eye(n_idx);
jacobian=[zeros(n_species-n_idx, n_idx); -I_n_idx];

else

% Case lof mutation that involves one elemental species

n_idx=n_idx-lof;
I_n_idx=eye(n_idx);
jacobian=[zeros(n_species-n_idx, n_idx); -I_n_idx];

end
end
