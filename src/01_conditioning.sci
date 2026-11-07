rand("seed", 0);
n = 100;    
p = 10;
puiss_bruit_col = 1e-5;            
puiss_bruit_y = 1e-3;   
// Génération de A avec colonnes corrélées 
A = zeros(n, p);
A(:,1) = rand(n, 1, 'normal');

for j = 2:p
    bruit = rand(n, 1, 'normal');
    A(:,j) = A(:,1) + puiss_bruit_col * bruit;
end
// Vérification du rang de A (doit être p)
rang_A = rank(A);
if rang_A < p then
    mprintf("Attention : A n''est pas de rang plein. Augmentez sigma.\n");
end

//  Vecteur beta_exact 
beta_exact = (1:p)';  
//  Construction de y = A*beta_exact + bruit 
bruit_y = puiss_bruit_y * rand(n, 1, 'normal');
y = A * beta_exact + bruit_y;

// Calcul des conditionnements;
cond_A=cond(A);
cond_ATA = cond(A' * A);

//  Affichage des résultats 
mprintf("puissance bruit (bruit entre colonnes) = %g\n", puiss_bruit_col);
mprintf("puiss_bruit_y = %g\n", puiss_bruit_y);
mprintf("cond(A)     = %g\n", cond_A);
mprintf("cond(A''*A) = %g\n", cond_ATA);

