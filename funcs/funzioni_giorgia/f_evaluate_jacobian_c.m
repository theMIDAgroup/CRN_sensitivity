%jacobiano di f rspetto a c
function jacobian = f_evaluate_jacobian_c(n_species, n_idx, lof)
% ordine=zeros(n_idx,1);

% for j=1:n_idx
%     ordine(j) = find(Nl(:, idx_basic_species(j)));
% end
if lof==0
I_n_idx=eye(n_idx);
jacobian=[zeros(n_species-n_idx, n_idx); -I_n_idx];%(ordine,:)];
else
n_idx=n_idx-lof;
I_n_idx=eye(n_idx);
jacobian=[zeros(n_species-n_idx, n_idx); -I_n_idx];%(ordine,:)];

end
