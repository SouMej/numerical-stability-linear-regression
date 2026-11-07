# Stabilité numérique et illusions algorithmiques en régression linéaire

**Projet d'Analyse Numérique Matricielle — INSAT**

Étude expérimentale sous Scilab de la stabilité numérique des méthodes de résolution
des problèmes de moindres carrés, avec un intérêt particulier pour le conditionnement,
la multicolinéarité et le paradoxe du « bon résidu ».

---

## Table des matières

1. [Présentation](#présentation)
2. [Objectifs](#objectifs)
3. [Organisation du projet](#organisation-du-projet)
   - [1. Génération des données et étude du conditionnement](#1-génération-des-données-et-étude-du-conditionnement)
   - [2. Comparaison des méthodes de résolution des moindres carrés](#2-comparaison-des-méthodes-de-résolution-des-moindres-carrés)
   - [3. Paradoxe du bon résidu et multicolinéarité](#3-paradoxe-du-bon-résidu-et-multicolinéarité)
   - [4. Changement d'échelle et normalisation](#4-changement-déchelle-et-normalisation)
   - [5. Expériences avancées d'instabilité numérique](#5-expériences-avancées-dinstabilité-numérique)
4. [Résultats principaux](#résultats-principaux)
5. [Conclusion générale](#conclusion-générale)
6. [Références](#références)

---

## Présentation

La régression linéaire constitue l'un des outils fondamentaux utilisés en data science.
Elle permet de modéliser la relation entre plusieurs variables explicatives et une
variable observée en estimant un vecteur de paramètres minimisant l'erreur entre les
valeurs observées et les valeurs prédites.

Le problème des moindres carrés s'écrit :

```
min_β  ||y - Aβ||₂
```

et sa solution vérifie les équations normales :

```
(AᵀA) β = Aᵀy
```

En arithmétique exacte, plusieurs formulations permettent de résoudre ce problème de
manière mathématiquement équivalente. Cependant, en arithmétique flottante, ces
formulations peuvent présenter des comportements numériques très différents.

Ce projet étudie expérimentalement ces différences et montre qu'un algorithme peut
produire un faible résidu tout en fournissant une estimation très éloignée des
véritables paramètres.

La question centrale du projet est donc :

> Peut-on obtenir un excellent ajustement des données tout en ayant un modèle
> numériquement faux ?

Les expériences réalisées montrent que la réponse est **oui**.

---

## Objectifs

Ce projet a pour objectifs de :

- manipuler des matrices de données sous Scilab ;
- générer des données synthétiques contrôlées ;
- étudier l'impact du conditionnement sur les problèmes de moindres carrés ;
- analyser l'effet de la multicolinéarité ;
- comparer différentes formulations mathématiques d'un même problème ;
- mesurer les erreurs sur les paramètres estimés ;
- comparer les erreurs sur les paramètres aux résidus ;
- mettre en évidence le paradoxe du « bon résidu » ;
- étudier l'influence du changement d'échelle des variables ;
- analyser l'intérêt de la normalisation et du centrage ;
- étudier l'influence du bruit sur les solutions numériques ;
- comparer différentes méthodes de factorisation QR ;
- construire des matrices mal conditionnées à l'aide de la SVD ;
- développer une analyse critique de la fiabilité des algorithmes numériques.

---

## Organisation du projet

Le projet est organisé autour de cinq axes principaux.

### 1. Génération des données et étude du conditionnement

Une matrice `A ∈ R^(100×10)` est générée avec des colonnes fortement corrélées afin de
créer une situation de multicolinéarité.

Le vecteur exact des paramètres est défini par :

```
β_exact = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]ᵀ
```

Les observations sont ensuite générées selon :

```
y = A β_exact + ε
```

où `ε` représente un bruit de faible amplitude.

La corrélation entre les colonnes est contrôlée par un paramètre `σ`. Plus `σ` diminue,
plus les colonnes deviennent similaires et plus la matrice devient mal conditionnée.

Le conditionnement de `A` et celui de `AᵀA` sont ensuite étudiés. Une observation
fondamentale est :

```
κ(AᵀA) ≈ κ(A)²
```

Ainsi, la formation de `AᵀA` peut amplifier fortement les problèmes de conditionnement.

### 2. Comparaison des méthodes de résolution des moindres carrés

Trois formulations du problème sont étudiées.

**Méthode 1 — Équations normales avec inversion explicite**

```
β₁ = (AᵀA)⁻¹ Aᵀy
```

Cette méthode consiste à calculer explicitement l'inverse de `AᵀA`.

**Méthode 2 — Résolution directe du système**

```
(AᵀA) β₂ = Aᵀy
```

Le système est directement résolu à l'aide de l'opérateur de résolution de systèmes
linéaires de Scilab (opérateur `\`, factorisation LU interne).

**Méthode 3 — Décomposition QR**

La matrice est factorisée sous la forme :

```
A = QR
```

puis le système réduit est résolu :

```
R β₃ = Qᵀy
```

Pour chaque méthode, deux grandeurs principales sont comparées :

- Erreur sur les paramètres : `||βᵢ - β_exact||`
- Résidu : `||A βᵢ - y||`

Cette comparaison permet d'observer que des méthodes mathématiquement équivalentes
peuvent présenter des comportements numériques très différents. Dans les expériences
réalisées, les trois méthodes produisent une erreur d'environ 18 sur `β` (alors que
`β_exact` vaut entre 1 et 10), avec des résidus tous très faibles (de l'ordre de
`10⁻³`) — ce qui illustre déjà l'illusion algorithmique : le modèle semble bien
ajuster les données alors que `β` est en réalité complètement faux.

### 3. Paradoxe du bon résidu et multicolinéarité

Cette expérience constitue l'un des résultats centraux du projet.

Le niveau de corrélation entre les colonnes de `A` est progressivement augmenté en
diminuant `σ`. Pour chaque valeur de `σ`, les grandeurs suivantes sont calculées :

- `κ(A)`
- `κ(AᵀA)`
- l'erreur sur `β`
- le résidu

Les résultats montrent un phénomène particulièrement important :

> Lorsque le conditionnement devient très élevé, l'erreur sur les coefficients peut
> augmenter de manière spectaculaire alors que le résidu reste faible.

Dans l'une des expériences, pour `σ = 10⁻⁴`, l'erreur sur les coefficients atteint
environ 209 alors que le résidu reste proche de 0,55.

Cela montre qu'un faible résidu ne garantit pas une estimation fiable des paramètres :

```
Faible résidu  =/=>  Paramètres fiables
```

Le résidu mesure la qualité de l'ajustement aux données observées, mais ne permet pas
à lui seul de mesurer la précision des paramètres estimés.

D'après la théorie des perturbations du problème des moindres carrés, l'erreur
relative sur la solution vérifie une relation du type :

```
||β̂ - β|| / ||β||  ≤  κ(A) · ( ||r|| / ||Aβ|| + ||δA|| / ||A|| )
```

où `r = Aβ̂ - y` est le résidu et `δA` une perturbation de la matrice `A`. Ainsi,
lorsque `A` est mal conditionnée, `κ(A)` amplifie fortement les erreurs présentes dans
les données, même si le résidu reste faible.

### 4. Changement d'échelle et normalisation

Cette partie étudie l'influence de l'échelle des variables explicatives sur la
stabilité numérique. Le facteur d'échelle appliqué à certaines colonnes de `A` varie
sur plusieurs ordres de grandeur, notamment entre `10⁻⁶` et `10⁶`.

Sans normalisation, une différence importante d'échelle entre les variables peut
fortement dégrader le conditionnement et provoquer une augmentation importante de
l'erreur sur les paramètres. Avec normalisation, le comportement numérique devient
plus stable.

Dans un cas concret (première colonne multipliée par `10⁶`, deuxième par `10⁻³`) :

| | Conditionnement | Erreur sur β | Résidu |
|---|---|---|---|
| Sans normalisation | ≈ 1,03 × 10¹⁰ | 1990,33 | 0,483 |
| Avec normalisation | 108,89 | 1,21 | 14,85 |

Cette expérience met en évidence l'importance du prétraitement des données : un
faible résidu ne suffit pas à valider une solution, la normalisation est
indispensable pour garantir la fiabilité des paramètres estimés.

### 5. Expériences avancées d'instabilité numérique

**5.1 Construction d'une matrice mal conditionnée par SVD**

Une matrice est construite à partir de sa décomposition en valeurs singulières :

```
A = U Σ Vᵀ
```

La plus petite valeur singulière est volontairement réduite à `σ_p = 10⁻¹⁰` afin de
rendre la matrice presque singulière et de contrôler directement son conditionnement.

**5.2 Effet du bruit sur les solutions**

Le niveau de bruit ajouté au vecteur `y` est progressivement augmenté et deux méthodes
sont comparées : équations normales et décomposition QR. Une perturbation `δy`
entraîne une erreur sur la solution approximativement donnée par :

```
||β̃ - β|| / ||β||  ≈  κ(A) · ||δy|| / ||y||
```

Des solutions de stabilisation sont également étudiées : SVD, SVD tronquée,
régularisation de Tikhonov / Ridge (`β_λ = (AᵀA + λI)⁻¹ Aᵀy`), et moindres carrés
pondérés.

**5.3 Gram-Schmidt, QR de Householder et équations semi-normales**

Trois approches sont comparées : QR par transformations de Householder, Gram-Schmidt
classique, et équations semi-normales.

- La méthode de Gram-Schmidt classique devient instable sur les matrices fortement
  mal conditionnées en raison de la perte d'orthogonalité (`QᵀQ ≠ I`).
- Les équations semi-normales (`RᵀR β = Aᵀy`) présentent également une amplification
  du conditionnement : `κ(RᵀR) ≈ κ(A)²`.
- La factorisation QR par transformations de Householder apparaît comme la méthode la
  plus fiable parmi les approches étudiées, car elle évite la formation de `AᵀA` et
  préserve numériquement l'orthogonalité de `Q`.

**5.4 Effet du mauvais centrage des données**

Des constantes de très grande amplitude sont ajoutées aux colonnes de la matrice :

```
A(:, j) = A₀(:, j) + 10¹⁰ · j
```

Cette opération crée une matrice fortement mal conditionnée. Après centrage :

```
A_c(:, j) = A(:, j) - mean(A(:, j))
```

le conditionnement diminue fortement et la solution obtenue devient beaucoup plus
proche de la solution exacte.

---

## Résultats principaux

- **Le conditionnement joue un rôle fondamental** : plus une matrice est mal
  conditionnée, plus la solution devient sensible aux perturbations des données et
  aux erreurs d'arrondi.
- **Les équations normales peuvent amplifier les erreurs** : la formation de `AᵀA`
  entraîne approximativement `κ(AᵀA) = κ(A)²`, ce qui rend cette approche
  particulièrement instable lorsque `A` est mal conditionnée. La borne d'erreur
  associée est de la forme `||β̂ - β|| / ||β|| ≲ c · κ(A)² · u`, où `u` est la
  précision machine — contre une dépendance en `κ(A)` seulement pour la méthode QR.
- **Un faible résidu ne garantit pas une bonne solution** : un modèle peut présenter
  un résidu très faible tout en ayant des coefficients très éloignés des valeurs
  exactes.
- **QR (Householder) est la méthode la plus fiable** parmi celles étudiées, car elle
  évite la formation de `AᵀA` et préserve mieux l'orthogonalité.
- **Le prétraitement des données est essentiel** : le centrage et la normalisation
  peuvent améliorer significativement le conditionnement et la fiabilité des calculs.
- **La convergence ne garantit pas la fiabilité** : un algorithme peut converger vers
  une solution présentant un faible résidu sans que les paramètres obtenus soient
  proches de la solution exacte.

---

## Conclusion générale

Ce projet met en évidence un point essentiel : un algorithme peut être
mathématiquement correct sans être numériquement fiable. En pratique, le
conditionnement du problème joue un rôle déterminant sur la qualité des résultats.

La décomposition QR apparaît comme la méthode la plus stable pour la résolution des
moindres carrés, tandis que les équations normales peuvent conduire à des erreurs
importantes, malgré un bon ajustement apparent.

Ces observations montrent que, en data science, il ne suffit pas d'obtenir un bon
résidu pour valider un modèle. Il est indispensable d'analyser le conditionnement, de
choisir une formulation adaptée et de préparer les données (centrage, mise à
l'échelle) afin de garantir des résultats fiables.

---

## Références

Les analyses théoriques présentées dans ce projet, en particulier celles relatives à
l'origine numérique des instabilités (conditionnement, amplification des erreurs,
stabilité des méthodes), s'appuient sur des résultats classiques issus de la
littérature en analyse numérique, notamment un chapitre de référence consacré au
problème des moindres carrés (conditionnement, stabilité numérique, comparaison des
méthodes de résolution, stabilité backward des algorithmes).
