// ---------- Méthode 1 : Équations normales (inversion directe) ----------
ATA = A' * A;
ATy = A' * y;
beta1 = inv(ATA) * ATy;

// ---------- Méthode 2 : Résolution directe du système ----------
beta2 = ATA \ ATy;

// ---------- Méthode 3 : Décomposition QR ----------
[Q, R] = qr(A);
beta3 = R \ (Q' * y);

// ---------- Calcul des erreurs sur beta ----------
err1 = norm(beta1 - beta_exact);
err2 = norm(beta2 - beta_exact);
err3 = norm(beta3 - beta_exact);

// ---------- Calcul des résidus ----------
res1 = norm(A * beta1 - y);
res2 = norm(A * beta2 - y);
res3 = norm(A * beta3 - y);

// ---------- Affichage ----------
mprintf("\n========== Résultats Partie 2 ==========\n");

mprintf("\n--- Erreurs |beta_i - beta_exact ---\n");
mprintf("Méthode 1 (inv(A''A))  : %e\n", err1);
mprintf("Méthode 2 (système \) : %e\n", err2);
mprintf("Méthode 3 (QR)        : %e\n", err3);

mprintf("\n--- Résidus |A*beta_i - y| ---\n");
mprintf("Méthode 1 (inv(A''A))  : %e\n", res1);
mprintf("Méthode 2 (système \) : %e\n", res2);
mprintf("Méthode 3 (QR)        : %e\n", res3);

mprintf("\n--- Comparaison des beta (5 premières composantes) ---\n");
mprintf("%-5s %-12s %-12s %-12s %-12s\n", "i", "beta_exact", "beta1", "beta2", "beta3");
for i = 1:10
    mprintf("%-5d %-12.6f %-12.6f %-12.6f %-12.6f\n", ...
        i, beta_exact(i), beta1(i), beta2(i), beta3(i));
end
