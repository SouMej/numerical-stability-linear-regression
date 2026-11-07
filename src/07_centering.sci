clc;
clear;
rand("seed", 0);

n = 100;
p = 10;

//  Données mal centrées
A = rand(n,p,'normal');

for j = 1:p
    A(:,j) = A(:,j) + 1e10 * j;
end

beta_exact = (1:p)';
bruit = 1e-20 * rand(n,1,'normal');
y = A * beta_exact + bruit;

//  Affichage extrait

disp("Extrait de A (mal centré) :");
disp(A(1:5,1:10));
disp("conditonnement de A ",cond(A));
//  QR + Equations normales

[Q,R] = qr(A,0);
beta_QR = R \ (Q' * y);

beta_NE = (A' * A) \ (A' * y);

//  Données centrées

A_c = zeros(n,p);

for j = 1:p
    A_c(:,j) = A(:,j) - mean(A(:,j));
end

y_c = A_c * beta_exact + bruit;

[Q_c,R_c] = qr(A_c,0);
beta_QR_c = R_c \ (Q_c' * y_c);


//  Affichage des betas
disp("----------------------------------------------");
disp(" i    Exact      QR        NE        Centré");
disp("----------------------------------------------");

for i = 1:p
    mprintf("%2d  %8.4f  %8.4f  %10.2e  %8.4f\n", ...
        i, beta_exact(i), beta_QR(i), beta_NE(i), beta_QR_c(i));
end

//  Erreurs

err_QR   = norm(beta_QR   - beta_exact);
err_NE   = norm(beta_NE   - beta_exact);
err_c    = norm(beta_QR_c - beta_exact);


//  Conditionnement

cond_A   = cond(A);
cond_c   = cond(A_c);


//  Affichage

disp("----------------------------------------------");
disp("Conditionnement : [Mal centré, Centré]");
disp([cond_A cond_c]);

disp("Erreurs : [QR, NE, Centré]");
disp([err_QR err_NE err_c]);


//  Graphiques
clf();

// Conditionnement 
subplot(1,2,1);
bar(log10([cond_A cond_c]));
a = gca();
a.x_ticks = tlist(["ticks","locations","labels"], ...
                  [1 2], ...
                  ["Mal centré","Centré"]);
title("Conditionnement (log10)");
ylabel("log10(cond(A))");

//  Erreurs 
subplot(1,2,2);
bar(log10([err_QR err_NE err_c]));
a = gca();
a.x_ticks = tlist(["ticks","locations","labels"], ...
                  [1 2 3], ...
                  ["QR","NE","Centré"]);
title("Erreur sur beta");
ylabel("log10(erreur)");

xgrid();
