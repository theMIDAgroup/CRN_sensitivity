function x_ss = f_draw_from_ssurf(cons_laws, rho, idx_basic_species, ...
                            par_loguni_sp_nocl)

%% Step 1. Initialize
[n_cons_law, n_species] = size(cons_laws);
x_ss = zeros(n_species, 1);
    

%% Step 2. Species that do not belong to any conservation law
idx_no_cl = find(sum(cons_laws, 1)==0);
x_ss(idx_no_cl) = 10.^((par_loguni_sp_nocl(2) - par_loguni_sp_nocl(1))...
        *rand(numel(idx_no_cl), 1)+par_loguni_sp_nocl(1));

%% Step 3. Species in cls but not basic
idx_in_cl = 1:n_species; idx_in_cl(idx_no_cl) = [];
idx_sp_todo = setdiff(idx_in_cl, idx_basic_species);
ordered_sp = idx_sp_todo(randperm(length(idx_sp_todo)));
aux_rho = rho;
for is = ordered_sp
    cl_is = find(cons_laws(:, is)>0);
    tmp_value = rand(1) * min(aux_rho(cl_is) ./ cons_laws(cl_is, is));
    x_ss(is) = tmp_value;
    aux_rho(cl_is) = aux_rho(cl_is) - cons_laws(cl_is, is) * tmp_value;
    if any(aux_rho<0)
        disp(is)
    end
end

%% Step 4. Basic species
for ir = 1:n_cons_law
    tmp_basic = intersect(idx_basic_species, find(cons_laws(ir, :)));
    x_ss(tmp_basic) = max(rho(ir) - cons_laws(ir, :)*x_ss, 0);
end



end