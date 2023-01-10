function jacobian = f_compute_analytic_jacobian_v(vm, n_species, idx_one)

n_reactions = size(vm, 1);
Jv_formula = zeros(n_reactions, n_species);
for ir = 1:n_reactions
    Jv_formula(ir, vm(ir, 2)') = vm(ir, 3);
    if vm(ir, 3) ~= idx_one
        Jv_formula(ir, vm(ir, 3)') = vm(ir, 2);
    end  
end

k_formula = vm(:, 1);

factor_formula = ones(n_reactions, 1);
factor_formula(vm(:, 2)==vm(:, 3)) = 2;

jacobian.Jv = Jv_formula;
jacobian.k = k_formula;
jacobian.factor = factor_formula; 

end

%% TODO possibili per velocizzare:
% (1) Eliminare il for nella definizione di Jv
% (2) Sfruttare le matrici sparse