function jacobian = f_evaluate_jacobian_k(vm, x, S, idx_basic_species)

%% This function permits to evaluate the jacobian of f=[S_2*v; Nx-c] with respect to k

%% Step 1. Initialization
n_reactions=size(vm,1);
num_cons_laws=numel(idx_basic_species);
S(idx_basic_species,:)=[]; 
Jv_formula_1=x(vm(:,2)); 
x_values=[x; 1];
Jv_formula_2=x_values(vm(:,3)); 

%% Step 2. Evaluate jacobian of f
jacobian=[(Jv_formula_1.*Jv_formula_2)'.*S; zeros(num_cons_laws, n_reactions)];

end
