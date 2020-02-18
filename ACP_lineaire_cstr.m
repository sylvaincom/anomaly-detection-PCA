% Script Matlab pour la détection de défaillance : ACP
% Version adaptée aux données du CSTR (non confidentielles)
% Version Matlab : R2018b
% Auteurs : Sylvain COMBETTES et Houssam L'GHOUL
% Dernière modification : 02/06/2019

clear all % pour supprimer toutes les variables contenues dans le Workspace
close all % pour fermer les figures précédentes
clc % clear command window

%% Importation des données (sans défaut)

% cd : pour voir la structure de l'arborescence des fichiers
% changement de répertoire vers le dossier où se situe le fichier de données que l'on souhaite importer
name = 'datav3.mat' ; % nom du fichier de données
data0 = load(name) ;
data = data0.cstr' ; % chaque variable est rangée dans une colonne
t0 = 100 ; % indice de début du régime permanent
data = data(t0+1:end, :) ; % matrice des données sans défaut en régime permanent

% On nomme les variables :
colNames = {'indice', 'Tc', 'T0', 'CAA', 'CAS', 'Fs', 'Fc', 'CA', 'T', 'FA'} ;
dataname = array2table(data, 'VariableNames', colNames) ;

dataname = head(dataname,5) % pour visualiser les 5 premières lignes du tableau data
% summary(T) % pour visualiser un résumé des données du tableau data

clear name data0 colNames

%% Traitement et nettoyage des données : obtention de la matrice des données Xb (sans défaut)

% Suppression des colonnes et lignes non pertinentes
data(:,1) = [] ;

% Suppression des lignes comportant des données manquantes
data = rmmissing(data,1) ;

Xb = data ; % matrice des données sans défaut
[N, m] = size(Xb) ; % N est le nombre d'observations et m est le nombre de variables

j = 6 ; % tester quelques valeurs de j (avec j compris entre 1 et m)
% Tracé de la variable j en fonction de l'indice de l'observation
figPos = get(0, 'defaultfigureposition') ;
width = 900 ;
height = 700 ;
figure('Position', [figPos(1), figPos(2), width, height]) ;
plot([1:N], Xb(:,j), 'o')
txt2 = 'en fonction de l''index de temps' ;
txt1 = sprintf('Tracé de la variable %d %s', j, txt2);
title({txt1, '(données sans défaut, non centrées et non réduites)'})
xlabel('Indice de l''observation')
ylabel(['Variable ', num2str(j)])

% % Boucle qui permet de tracer successivement chaque variable :
% for j=1:m
%     figure('Position', [figPos(1), figPos(2), width, height]) ;
%     plot([1:N], Xb(:,j), 'o')
%     txt2 = 'en fonction de l''index de temps' ;
%     txt1 = sprintf('Tracé de la variable %d %s', j, txt2);
%     title({txt1, '(données sans défaut, non centrées et non réduites)'})
%     xlabel('Indice de l''observation')
%     ylabel(['Variable ', num2str(j)])
% end

clear data j txt2 txt1

%% Vérification de l'hypothèse : les données suivent une loi normale multivariable

j = 6 ; % tester quelques valeurs de j (avec j compris entre 1 et m)
% Vérification de l'hypothèse : la variable j suit une loi normale
figure('Position', [figPos(1), figPos(2), width, height]) ;
histfit(Xb(:,j))
title({ sprintf('Histogramme de la variable %d',j) ; '(données sans défaut, non centrées et non réduites)' })
legend([{sprintf(' histogramme de la variable %d', j); ' densité normale estimée à partir des données'}])

% % Boucle qui permet de tracer successivement les histogrammes de chaque variable
% for j=1:m
%    figure('Position', [figPos(1), figPos(2), width, height]) ;
%    histfit(Xb(:,j))
%    title({ sprintf('Histogramme de la variable %d',j) ; '(données sans défaut, non centrées et non réduites)' })
%    legend([{sprintf(' histogramme de la variable %d', j); ' densité normale estimée à partir des données'}])
% end

r = 0.05 ; % risque
pval = [] ;
hres = [] ;
for j=1:m
    x = Xb(:,j) ;
    [h, p] = chi2gof(x, 'Alpha', r) ;
    pval = [pval ; p ] ;
    hres = [hres ; h ] ;
end
figure('Position', [figPos(1), figPos(2), width, height]) ;
plot([1:m], pval, 'ob', [1:m], hres, 'xr')
title({ sprintf('Test du khi-deux avec le risque %d %%', roundn(r*100,-4)) ; '(données sans défaut)' })
legend([{'p-value' ; 'conclusion'}])
% The returned value h = 0 indicates that chi2gof does not reject the null hypothesis at the default 5% significance level. 

pd = makedist('Normal');
rng default;  % for reproducibility
x = random(pd,100,1);
h = chi2gof(x)

%% Obtention de X : matrice des données centrées et réduites (sans défaut)

X = zeros(N,m) ; % initialisation de X qui est la matrice des données centrée réduite
mean0 = mean(Xb) ; % vecteur qui contient les valeurs pour centrer la matrice des données
std0 = std(Xb)' ; % vecteur qui contient les valeurs pour réduire la matrice des données

for j = 1:m
    X(:,j) = (Xb(:,j)-mean0(j))./(std0(j)) ;
end
% mean(X(:,2)) % -6.2937e-12 qui est le zéro informatique
% std(X(:,2)) % = 1.0000
% X est bien centrée (moyenne nulle) et réduite (écart-type valant 1)

clear j

%% Obtention de l'espace principal Pp et de l'espace résiduel Pr

S = (1/(N-1))*X'*X ; % matrice de variance / covariance

[V0,D0] = eig(S) ; % V0 est la matrice des vecteurs propres (non triés), D0 est la matrice des valeurs propres (non triées)
D0 = real(D0) ; % c'est plus lisible car par exemple la composante 5.1798 + 0.0000i devient 5.1798
vectD0 = diag(D0) ; % vecteur des valeurs propres (non triees)
[vectD, indice_sort] = sort(vectD0, 'descend') ;
% vectD : vecteur des valeurs propres triées par ordre décroissant
% indice_sort contient les indices des valeurs triées, de sorte à reproduire le tri
V = V0(:,indice_sort) ; % vecteurs propres triés par ordre décroissant selon leurs valeurs propres correspondantes

clear S V0 D0 vectD0 indice_sort

l = 0 ; % initialisation de l
% l est le nombre de composantes retenues dans le modèle ACP (dimension du sous-espace des composantes principales
PCV = 0 ; % initialisation du pourcentage cumulé de la variance totale (PCV)
info = 0.9 ; % on choisit info arbitrairement (pour que l'espace résiduel contienne au moins un vecteur)
while PCV < info
    l = l+1 ;
    PCV = (sum(vectD(1:l)))./(sum(vectD)) ;
end
l
PCV

Pp = V(:, 1:l) ; % espace principal
Pr = V(:, l+1:m) ; % espace résiduel

% Tracé de la PCV en fonction de l
l2 = 0 ;
PCV2 = 0 ;
abs_l2 = [] ;
ord_PCV2 = [] ;
for l2 = 1:m
    PCV2 = (sum(vectD(1:l2)))./(sum(vectD)) ;
    abs_l2 = [abs_l2, l2] ;
    ord_PCV2 = [ord_PCV2, PCV2] ;
end
figure('Position', [figPos(1), figPos(2), width, height]) ;
plot(abs_l2, ord_PCV2, 'o')
yline(info, 'r', 'LineWidth', 1) ;
title({'Evolution du PCV avec le nombre de composantes principales l' ; '(données sans défaut)'})
xlabel('Nombre de composantes principales l')
ylabel('Pourcentage cumulé de la variance totale (PCV)')
legend([{' PCV', sprintf(' cible : %d%% de l''information', round(info*100,2))}])

clear l2 PCV2 abs_l2 ord_PCV2

%% Calcul de Xp qui est l'estimation de X par l'ACP : Xp est la projection de X sur l'espace principal Pp

Xp = zeros(N,m) ; % initialisation
for i = 1:N
    Y = 0 ;
    for j = 1:l
        Y = Y + Pp(:,j)*X(i,:)*Pp(:,j) ;
    end
    Xp(i,:) = Y' ;
end

clear i j Y

j = 6 ; % tester quelques valeurs de j (avec j compris entre 1 et m)
% Tracé de l'erreur de l'estimation de la variable j par sa projection sur l'espace propre Pp
figure('Position', [figPos(1), figPos(2), width, height]) ;
plot([1:N], X(:,j)-Xp(:,j), 'o')
yline(0,'r','LineWidth',2) ;
txt1 = sprintf('Erreur de l''estimation de la variable X(:,%d) par sa projection Xp(:,%d) sur l''espace propre', j, j);
title({txt1, '(données sans défaut)'})
xlabel('Indice de l''observation')
ylabel([sprintf('Erreur X(:,%d)-Xp(:,%d)', j, j)])
legend(sprintf(' erreur X(:,%d)-Xp(:,%d)', j, j), ' droite : sans erreur')

clear j txt1

%% Calcul de la SPE (erreur quadratique de prédiction)

SPE = zeros(N,1) ; % initialisation de SPE : vecteur colonne dont chaque ligne
% contient la SPE de l'observation numero i
for i = 1:N
    SPE(i) = (X(i,:)-Xp(i,:))*(X(i,:)-Xp(i,:))' ;
end

clear i

% SPE = diag((X-Xp)*(X-Xp)') ; % mais c'est plus coûteux en espace mémoire

% Tracé de la SPE en fonction de l'indice de l'observation
figure('Position', [figPos(1), figPos(2), width, height]) ;
plot([1:N], SPE, 'o')
title({'Tracé de la SPE en fonction de l''indice de l''observation' ; '(données sans défaut)'})
xlabel('Indice de l''observation')
ylabel('SPE')

% xlswrite('spe.xlsx',SPE); % exportation du vecteur SPE en fichier spe.xlsx (Excel)

%% Calcul du seuil de détection delta2

theta = zeros(3,1) ; % initialisation
lambdar = vectD(l+1:end) ; % vecteur des valeurs propres residuelles de S triée
for i = 1:3
    theta(i) = sum(lambdar.^i) ;
end

h0 = 1 - [(2*theta(1)*theta(3))/(3*theta(2)^2)] ;
calpha = theta(1) * [(mean(SPE)^2/theta(1))^h0 - 1 - (theta(2)*h0*(h0-1)/theta(1)^2)] / [sqrt(2*theta(2)*h0^2)] ;
delta2 = theta(1) * [ [calpha*sqrt(2*theta(2)*h0^2)]/theta(1) + 1 + [theta(2)*h0*(h0-1)]/(theta(1)^2) ]^(1/h0) 

clear theta lambdar i h0 calpha

% Comparaison de la SPE avec le seuil de détection
figure('Position', [figPos(1), figPos(2), width, height]) ;
ylim = Inf ; % valeur maximale des ordonnées que l'on affiche
plot([1:N], SPE, 'o')
yline(delta2, 'r', 'LineWidth', 2) ;
axis([-Inf Inf -Inf ylim])
title({'Comparaison de la SPE avec le seuil de détection' ; '(données sans défaut)'})
xlabel('Indice de l''observation')
ylabel('SPE')
legend(' SPE', ' seuil de détection')

%% Importation et traitement des données défaillantes

name2 = 'dataDefautv3.mat' ; % nom du fichier comportant des données défaillantes
data02 = load(name2) ;
data2 = data02.cstr1' ; % chaque variable est rangée dans une ligne
% rappel : t0 = 200 : indice de début du régime permanent
data2 = data2(t0+1:end, :) ; % matrice des données avec défaut en régime permanent

% T2 = array2table(data2) ;
% head(T2,10) % pour visualiser les 10 premières lignes du tableau data
% summary(T2) % pour visualiser un résumé des données du tableau data

data2(:,1) = [] ;
data2 = rmmissing(data2,1) ;

Xb2 = data2 ; % matrice des données avec défaut
[N2, m2] = size(Xb2) ; % N2 est le nombre d'observations et m2 est le nombre de variables

% Les matrices Xb2 et Xb ont-elles les mêmes dimensions ?
N2 == N ;
m2 == m ;

j = 6 ; % tester quelques valeurs de j (avec j compris entre 1 et m)
% Tracé de la variable j en fonction de l'indice de l'observation
figure('Position', [figPos(1), figPos(2), width, height]) ;
plot([1:N], Xb2(:,j), 'o')
txt2 = 'en fonction de l''index de temps' ;
txt1 = sprintf('Tracé de la variable %d %s', j, txt2);
title({txt1, '(données avec défaut)'})
xlabel('Indice de l''observation i')
ylabel(['Variable ', num2str(j)])

% for j=1:m2
%     figure('Position', [figPos(1), figPos(2), width, height]) ;
%     plot([1:N], Xb2(:,j), 'o')
%     txt2 = 'en fonction de l''index de temps' ;
%     txt1 = sprintf('Tracé de la variable %d %s', j, txt2);
%     title({txt1, '(données avec défaut, non centrées et non réduites)'})
%     xlabel('Indice de l''observation i')
%     ylabel(['Variable ', num2str(j)])
% end

% On centre et on réduit avec les moyennes et écart-types des données sans
% défaut mean0 et std0 :
X2 = zeros(N2,m2) ;
for j = 1:m2
    X2(:,j) = (Xb2(:,j)-mean0(j))./(std0(j)) ;
end

clear data02 data2 j

%% Détection des observations défaillantes

Xp2 = zeros(N2,m2) ;
for i = 1:N2
    Y = 0 ;
    for j = 1:l
        Y = Y + Pp(:,j)*X2(i,:)*Pp(:,j) ;
    end
    Xp2(i,:) = Y ;
end

clear Y

j = 6 ; % tester quelques valeurs de j (avec j compris entre 1 et m)
% Tracé de l'erreur de l'estimation de la variable j par sa projection sur l'espace propre Pp
figure('Position', [figPos(1), figPos(2), width, height]) ;
plot([1:N], X2(:,j)-Xp2(:,j), 'o')
yline(0,'r','LineWidth',2) ;
txt1 = sprintf('Erreur de l''estimation de la variable X(:,%d) par sa projection Xp(:,%d) sur l''espace propre', j, j) ;
title({txt1, '(données avec défaut)'})
xlabel('Indice de l''observation')
ylabel([sprintf('Erreur X(:,%d)-Xp(:,%d)', j, j)])
legend(sprintf(' erreur X(:,%d)-Xp(:,%d)', j, j), ' droite : sans erreur')

SPE2 = zeros(N2,1) ; % initialisation de SPE : vecteur colonne dont chaque ligne contient la SPE de l'observation numéro i
for i = 1:N2
    SPE2(i) = (X2(i,:)-Xp2(i,:))*(X2(i,:)-Xp2(i,:))' ;
end

clear i

% Tracé : comparaison de la SPE avec le seuil de détection
figure('Position', [figPos(1), figPos(2), width, height]) ;
ylim = Inf ; % valeur maximale des ordonnées que l'on affiche, prendre Inf puis 3
plot([1:N2], SPE2, 'o')
yline(delta2, 'r', 'LineWidth', 2) ;
axis([-Inf Inf -Inf ylim])
title({'Comparaison de la SPE avec le seuil de détection' ;  '(données avec défaut)'})
xlabel('Indice de l''observation')
ylabel('SPE')
legend('SPE', 'seuil de détection')

clear ylim

defaillance = [] ;
for i = 1:N2
    if SPE2(i) > delta2
        defaillance = [defaillance; i, SPE2(i), Xb2(i,6)] ;
    end
end

clear i

pourcentage_de_mesures_defaillantes = round(length(defaillance)/N2*100,0) % 93

% defaillants = array2table(defaillance) ;
% head(defaillants) ;
% summary(defaillants) ;

% figure('Position', [figPos(1), figPos(2), width, height]) ;
% plot(defaillance(:,1), defaillance(:,3), 'o')
% title('Réprésentation des données défaillantes')
% xlabel('Indice de l''observation i')
% ylabel('Valeur mesurée FA')

%% D'où provient la défaillance ?

% Rappel : il faut regarder les données centrées et réduites.

abs = [] ;
moy = [] ;
moy2 = [] ;

for j = 1:m
    abs = [abs, j] ;
    moy = [moy, mean(X(:,j))] ;
    moy2 = [moy2, mean(X2(:,j))] ;
end

clear j

% Tracé de la moyenne en fonction de l'index de variable
figure('Position', [figPos(1), figPos(2), width, height]) ;
plot(abs, moy, 'og', abs, moy2, 'xr')
title({sprintf('Tracé de la moyenne en fonction de l''index de variable') ; '(données centrées et réduites)'})
xlabel('Numéro de variable')
ylabel('Moyenne')
legend('données sans défaut', 'données avec défaut')

j = 6 ; % tester j=6, j=7 et j=8
% Superposition des graphes (données centrées et réduites)
figure('Position', [figPos(1), figPos(2), width, height]) ;
plot([1:N], X(:,j), 'og', [1:N], X2(:,j), 'xr')
title({sprintf('Tracé de la variable %d %s', j, 'en fonction de l''index de temps') ; '(données centrées et réduites)'})
xlabel('Indice de l''observation')
ylabel(['Variable ', num2str(j)])
legend('données sans défaut', 'données avec défaut')

% for j=1:m2
%     figure('Position', [figPos(1), figPos(2), width, height]) ;
%     plot([1:N], X(:,j), 'og', [1:N], X2(:,j), 'xr')
%     title({sprintf('Tracé de la variable %d %s', j, 'en fonction de l''index de temps') ; '(données centrées et réduites)'})
%     xlabel('Indice de l''observation')
%     ylabel(['Variable ', num2str(j)])
%     legend('données sans défaut', 'données avec défaut')
% end

clear j