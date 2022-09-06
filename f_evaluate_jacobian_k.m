%questo jacobiano si calcola con i valori di x
function jacobian = f_evaluate_jacobian_k(vm, x, S, idx_basic_species)

n_reactions=size(vm,1); %numero di reazioni

S(idx_basic_species,:)=[]; %creo S2, matrice con righe delle specie non elementari, cancellando righe elementari

Jv_formula_1=x(vm(:,2)); % considero solo i valori di x negli indici dati dalla seconda colonna del vettore vm

%Non so se idx_one=n_species+1 sempre, ma nel caso vale questo:
x_values=[x; 1];%creo un nuovo vettore che ha entrata idx_one (che qui dovrebbe essere n_species+1) pari a 1
Jv_formula_2=x_values(vm(:,3)); % considero i valori dati dalla terza colonna. Quando vm(i,3)=idx_one, allora x(idx_one)=1

num_cons_laws=numel(idx_basic_species); %conto quanto conservation laws ci sono, che mi dà p

%jacobian è una matrice zeros(n_species, n_reactions);
jacobian=[(Jv_formula_1.*Jv_formula_2)'.*S; zeros(num_cons_laws, n_reactions)];

%la prima parte di matrice dello jacobiano è data da S, dove la colonna j è
%moltiplicata per la componente j dello ja di v, vista come vettore.
%possiamo considerarlo vettore in quanto il vero jac è una matrice
%diagonale.
%la diagonale coincide proprio con il prodotto puntuale tra Jv_fromula_1 e
%Jv_formula_2
%non so quanto sia chiaro quello che ho scritto

end
