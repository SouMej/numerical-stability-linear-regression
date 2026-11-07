

clear;
clc;
rand("seed", 0);

// ---------- Paramètres ----------
n = 100;
p = 10;

puiss_bruit_col = 0.1;
puiss_bruit_y   = 0.05;

// ---------- Génération de A corrélée ----------
A = zeros(n, p);
A(:,1) = rand(n,1,'normal');

for j = 2:p
    A(:,j) = A(:,1) + puiss_bruit_col * rand(n,1,'normal');
end

// ---------- Vrai beta ----------
beta_exact = (1:p)';
disp(beta_exact)

// ---------- y ----------
y = A * beta_exact + puiss_bruit_y * rand(n,1,'normal');

// =================================================
// 1. Changement d’échelle (IMPORTANT)
// =================================================

A_scaled = A;

// On modifie fortement certaines colonnes
A_scaled(:,1) = 1e6 * A_scaled(:,1);   // énorme
A_scaled(:,2) = 1e-3 * A_scaled(:,2);  // très petit

// =================================================
// 2. Résolution AVANT normalisation
// =================================================

beta_scaled = inv(A_scaled'*A_scaled) * A_scaled' * y;

// erreurs et résidus
err_scaled = norm(beta_scaled - beta_exact);
res_scaled = norm(A_scaled * beta_scaled - y);

// conditionnement
cond_scaled = cond(A_scaled);

// =================================================
// 3. Normalisation des données
// =================================================

A_norm = A_scaled;

for j = 1:p
    mu = mean(A_scaled(:,j));
    sigma = stdev(A_scaled(:,j));
    A_norm(:,j) = (A_scaled(:,j) - mu) / sigma;
end

// =================================================
// 4. Résolution APRÈS normalisation
// =================================================

beta_norm = inv(A_norm'*A_norm) * A_norm' * y;

// erreurs et résidus
err_norm = norm(beta_norm - beta_exact);
res_norm = norm(A_norm * beta_norm - y);

// conditionnement
cond_norm = cond(A_norm);

// =================================================
// 5. AFFICHAGE - CORRECTION ICI
// =================================================

mprintf("\n===== SANS NORMALISATION =====\n");
mprintf("Conditionnement = %g\n", cond_scaled);
mprintf("Erreur beta     = %g\n", err_scaled);
mprintf("Résidu          = %g\n\n", res_scaled);

mprintf("===== AVEC NORMALISATION =====\n");
mprintf("Conditionnement = %g\n", cond_norm);
mprintf("Erreur beta     = %g\n", err_norm);
mprintf("Résidu          = %g\n", res_norm);


// =================================================
// 6. ETUDE : variation LARGE échelle (log)
// =================================================

// facteurs de 10^-6 à 10^6
facteurs = 10.^(-6:1:6);

erreurs = zeros(facteurs);
erreurs_norm = zeros(facteurs);

for k = 1:length(facteurs)
    
    facteur = facteurs(k);
    
    // =========================
    // CAS NON NORMALISÉ
    // =========================
    A_test = A;
    
    A_test(:,1) = facteur * A_test(:,1);
    A_test(:,2) = (1/facteur) * A_test(:,2);
    
    beta_test = inv(A_test'*A_test) * A_test' * y;
    
    erreurs(k) = norm(beta_test - beta_exact) / norm(beta_exact); // erreur relative
    
    
    // =========================
    // CAS NORMALISÉ
    // =========================
    A_test_norm = A_test;
    
    for j = 1:p
        mu = mean(A_test(:,j));
        sigma = stdev(A_test(:,j));
        A_test_norm(:,j) = (A_test(:,j) - mu) / sigma;
    end
    
    beta_test_norm = inv(A_test_norm'*A_test_norm) * A_test_norm' * y;
    
    erreurs_norm(k) = norm(beta_test_norm - beta_exact) / norm(beta_exact);
end

// =================================================
// 7. GRAPHIQUE COMPARATIF
// =================================================

clf();
plot(facteurs, erreurs, '-o');         // sans normalisation
plot(facteurs, erreurs_norm, '-r');    // avec normalisation
xgrid();

xlabel("Facteur déchelle (log)");
ylabel("Erreur relative sur beta");
title("Impact du changement déchelle (avec / sans normalisation)");

legend(["Sans normalisation", "Avec normalisation"]);

// échelle log-log
set(gca(), "log_flags", "ll");
