// =============================================
// PARTIE 3 : Le paradoxe du "bon résidu"
// =============================================

clear;
clc;
rand("seed", 0);

// ---------- Paramètres ----------
n = 100;    
p = 10;
puiss_bruit_y = 0.05;              // bruit sur y (fixe)
beta_exact = (1:p)';               // vrai beta

// ---------- Vecteur des puissances de bruit à tester ----------
// Plus sigma est petit, plus les colonnes sont corrélées
sigmas = [1.0, 0.5, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001, 1e-5, 1e-6];

// ---------- Initialisation des tableaux de résultats ----------
n_sigmas = length(sigmas);
cond_A = zeros(1, n_sigmas);
cond_ATA = zeros(1, n_sigmas);
erreur_beta = zeros(1, n_sigmas);
residu = zeros(1, n_sigmas);

// ---------- Boucle sur les différentes corrélations ----------
for k = 1:n_sigmas
    
    sigma = sigmas(k);
    
    // Génération de A avec corrélation contrôlée par sigma
    A = zeros(n, p);
    A(:,1) = rand(n, 1, 'normal');
    
    for j = 2:p
        bruit = rand(n, 1, 'normal');
        A(:,j) = A(:,1) + sigma * bruit;
    end
    
    // Construction de y
    bruit_y = puiss_bruit_y * rand(n, 1, 'normal');
    y = A * beta_exact + bruit_y;
    
    // Résolution par QR (méthode la plus stable)
    [Q, R] = qr(A);
    beta_estime = R \ (Q' * y);
    
    // Calcul des métriques
    cond_A(k) = cond(A);
    cond_ATA(k) = cond(A' * A);
    erreur_beta(k) = norm(beta_estime - beta_exact);
    residu(k) = norm(A * beta_estime - y);
    
end

// =============================================
// AFFICHAGE DES RÉSULTATS
// =============================================

mprintf("\n");
mprintf("=====================================================================\n");
mprintf("PARTIE 3 : Évolution des métriques avec l''augmentation de corrélation\n");
mprintf("=====================================================================\n\n");

mprintf("%-12s %-15s %-15s %-15s %-15s\n", ...
        "sigma", "cond(A)", "cond(A''A)", "Erreur beta", "Résidu");
mprintf("---------------------------------------------------------------------\n");

for k = 1:n_sigmas
    mprintf("%-12.6f %-15.4e %-15.4e %-15.4e %-15.4e\n", ...
            sigmas(k), cond_A(k), cond_ATA(k), erreur_beta(k), residu(k));
end

// =============================================
// GRAPHIQUES
// =============================================

// Figure 1 : Conditionnement
scf(1);
clf(1);
subplot(2,1,1);
plot(log10(sigmas), log10(cond_A), 'b-o', 'LineWidth', 2);
xlabel('log10(sigma)');
ylabel('log10(cond(A))');
title('Évolution du conditionnement de A avec la corrélation');
xgrid(1);

subplot(2,1,2);
plot(log10(sigmas), log10(cond_ATA), 'r-o', 'LineWidth', 2);
xlabel('log10(sigma)');
ylabel('log10(cond(A''A))');
title('Évolution du conditionnement de A''A avec la corrélation');
xgrid(1);

// Figure 2 : Erreur sur beta et Résidu
scf(2);
clf(2);
subplot(2,1,1);
plot(log10(sigmas), log10(erreur_beta), 'g-o', 'LineWidth', 2);
xlabel('log10(sigma)');
ylabel('log10(Erreur sur beta)');
title('Évolution de l''erreur sur beta avec la corrélation');
xgrid(1);

subplot(2,1,2);
plot(log10(sigmas), log10(residu), 'm-o', 'LineWidth', 2);
xlabel('log10(sigma)');
ylabel('log10(Résidu)');
title('Évolution du résidu avec la corrélation');
xgrid(1);
// Figure 3 : Comparaison directe Erreur vs Résidu
scf(3);
clf(3);
plot(log10(sigmas), log10(erreur_beta), 'g-o', 'LineWidth', 2);
plot(log10(sigmas), log10(residu), 'm--s', 'LineWidth', 2);
xlabel('log10(sigma)');
ylabel('log10(Valeur)');
title('Comparaison : Erreur sur beta (trait plein) vs Résidu (pointillés)');
legend(['Erreur beta', 'Résidu'], 'location', 'northwest');
xgrid(1);



