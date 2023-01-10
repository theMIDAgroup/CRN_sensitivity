%jacobiano di v rispetto a x
function jacobian = f_compute_analytic_jacobian_v(vm, n_species, idx_one)

n_reactions = size(vm, 1); % numero reazioni è lo stesso di vm
Jv_formula = zeros(n_reactions, n_species); % creo matrice nulla per lo jacobiano
for ir = 1:n_reactions % per ogni riga metti in posizione corretta la concentrazione specie
    Jv_formula(ir, vm(ir, 2)') = vm(ir, 3); 
    if vm(ir, 3) ~= idx_one  % se nella stessa riga sonon coinvolte più reazioni, metti in posizione corretta anche l'altra
        Jv_formula(ir, vm(ir, 3)') = vm(ir, 2);
    end  
end

k_formula = vm(:, 1);

factor_formula = ones(n_reactions, 1);
factor_formula(vm(:, 2)==vm(:, 3)) = 2;

jacobian.Jv = Jv_formula; %matrice dello jacobiano di v rispetto a x
jacobian.k = k_formula; %coefficienti delle rate constants
jacobian.factor = factor_formula; %esponenti delle x

end

%% TODO possibili per velocizzare:
% (1) Eliminare il for nella definizione di Jv
% (2) Sfruttare le matrici sparse