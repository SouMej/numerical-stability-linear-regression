clc;
clear;
rand("seed", 0);
n = 100;
p = 10;
// Matrice mal conditionnée en changeant valeur singulière la minimale 
A = rand(n,p);
[U,S,V] = svd(A, "e");
s = diag(S);
s(p) = 1e-12; 
S_mod = diag(s);

A_bad = U * S_mod * V';
// Beta exact
beta_exact = (1:p)';
// Bruit (logarithmique)
noise_levels = logspace(-20, log10(5), 40);
// Stockage erreurs
err_NE = zeros(noise_levels);
err_QR = zeros(noise_levels);
// Boucle
for i = 1:length(noise_levels)
    sigma = noise_levels(i);
    bruit_y = sigma * rand(n,1,'normal');
    y = A_bad * beta_exact + bruit_y;
    // NE
    beta_NE = (A_bad' * A_bad) \ (A_bad' * y);
    // QR
    [Q,R] = qr(A_bad,0);
    beta_QR = R \ Q' * y;
    // erreurs
    err_NE(i) = norm(beta_NE - beta_exact);
    err_QR(i) = norm(beta_QR - beta_exact);
end
// Affichage de courbes
clf();

plot(log10(noise_levels), log10(err_NE), '-o');
plot(log10(noise_levels), log10(err_QR), '-x');

legend("Normal Equations", "QR");

xlabel("log10(bruit y)");
ylabel("log10(erreur beta)");
title("Instabilité NE vs QR (bruit variable)");

