clc;
clear;
rand("seed", 0);
n = 100;
p = 10;
// Matrice mal conditionnée( par la méthode des valeurs singulières )
A = rand(n,p);
[U,S,V] = svd(A, "e");
s = diag(S);
s(p) = 1e-10; 
S_mod = diag(s);
A = U * S_mod * V';
// Beta exact
beta_exact = (1:p)';
// y bruité ( mais très faible bruit )
bruit_y = 1e-20 * rand(n,1,'normal');
y = A * beta_exact + bruit_y;
// Affichage A
disp("Extrait de A :");
disp(A(1:5,1:10));
// QR (Householder)
[Q,R] = qr(A,0);
beta_QR = R \ (Q' * y);
// Semi-normales (RᵀR)
beta_semi = (R'*R) \ (A' * y);
// Gram-Schmidt classique
Q_gs = zeros(n,p);
R_gs = zeros(p,p);

for j = 1:p
    v = A(:,j);

    for i = 1:j-1
        R_gs(i,j) = Q_gs(:,i)' * A(:,j);
        v = v - R_gs(i,j) * Q_gs(:,i);
    end

    R_gs(j,j) = norm(v);
    Q_gs(:,j) = v / R_gs(j,j);
end
beta_GS = R_gs \ (Q_gs' * y);

// Tableau pour les betas 
disp("----------------------------------------------");
disp(" i    Exact      QR        Semi        GS");
disp("----------------------------------------------");
for i = 1:p
    mprintf("%2d  %8.4f  %8.4f  %10.2e  %8.4f\n", ...
        i, beta_exact(i), beta_QR(i), beta_semi(i), beta_GS(i));
end
// Erreurs
err_QR   = norm(beta_QR   - beta_exact);
err_semi = norm(beta_semi - beta_exact);
err_GS   = norm(beta_GS   - beta_exact);
disp("----------------------------------------------");
disp("Erreurs : [QR   Semi   GS]");
disp([err_QR err_semi err_GS]);
// Graphique des erreurs en batons 
clf();
errors = [err_QR err_semi err_GS];
// échelle log pour lisibilité
bar(log10(errors));
a = gca();
a.x_ticks = tlist(["ticks","locations","labels"], ...
                  [1 2 3], ...
                  ["QR","Semi","GS"]);

ylabel("log10(erreur)");
title("Comparaison des méthodes");
xgrid();
