% Script Matlab pour la détection de défaillance
% ACP à noyau gaussien
% Version adaptée aux données du CSTR (non confidentielles)
% Version Matlab : R2018b
% Auteurs : Sylvain COMBETTES et Houssam L'GHOUL
% Dernière modification : 02/06/2019

clear all % pour supprimer toutes les variables contenues dans le Workspace
close all % pour fermer les figures précédentes
clc % clear command window

%__________________________________________________________________________
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

%__________________________________________________________________________
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

%__________________________________________________________________________
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

%__________________________________________________________________________
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

%__________________________________________________________________________
%% Définition du noyau gaussien pour les données sans défaut (noté K)

sigma = 3 ; % il faut jouer sur cette valeur
% Exemples :
% - lorsque sigma augmente, l diminue fortement (tout le reste constant)
% - lorsque sigma augmente, le conditionnement augmente fortement

K = exp( -(dist(X,X').^2)./(2*sigma^2) ) ; % calcul matriciel

% Bon choix de sigma ?
% mean(mean(K))
% median(median(K))

% conditionnement de K
cond_K = norm(inv(K))*norm(K)

%__________________________________________________________________________
%% Obtention de l'espace principal Pp et de l'espace résiduel Pr

[V0,D0] = eig(K/(N-1)) ; % V0 est la matrice des vecteurs propres (non triés), D0 est la matrice des valeurs propres (non triées)
D0 = real(D0) ; % c'est plus lisible car par exemple la composante 5.1798 + 0.0000i devient 5.1798
vectD0 = diag(D0) ; % vecteur des valeurs propres (non triees)
% cond_K_2 = max(abs(vectD0))/min(abs(vectD0)) % conditionnement de K (autre méthode de calcul)
[vectD, indice_sort] = sort(vectD0, 'descend') ;
% vectD : vecteur des valeurs propres triées par ordre décroissant
% indice_sort contient les indices des valeurs triées, de sorte à reproduire le tri
V = V0(:,indice_sort) ; % vecteurs propres triés par ordre décroissant selon leurs valeurs propres correspondantes

clear S V0 D0 vectD0 indice_sort

l = 0 ; % initialisation de l
% l est le nombre de composantes retenues dans le modèle ACP (dimension du sous-espace des composantes principales
PCV = 0 ; % initialisation du pourcentage cumulé de la variance totale (PCV)
info = 0.95 ; % on choisit info arbitrairement (pour que l'espace résiduel contienne au moins un vecteur)
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
for l2 = 1:N
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

%__________________________________________________________________________
%% Calcul de la SPE (erreur quadratique de prédiction)

% vectD(1:l) : vecteur des l plus grandes valeurs propres triées par ordre decroissant
vp = 1./vectD(1:l) ;
LAMBDA = diag(vp)/(N-1) ;

SPE = zeros(N,1) ; % initialisation de SPE : vecteur colonne dont chaque ligne
% contient la SPE de l'observation numero i.

for i = 1:N
    SPE(i) = K(i,i) - K(i,:)*Pp*LAMBDA*Pp'*K(i,:)' ;
end

clear i

% Tracé de la SPE en fonction de l'indice de l'observation
figure('Position', [figPos(1), figPos(2), width, height]) ;
plot([1:N], SPE, 'o')
title({'Tracé de la SPE en fonction de l''indice de l''observation' ; '(données sans défaut)'})
xlabel('Indice de l''observation')
ylabel('SPE')

% xlswrite('spe.xlsx',SPE); % exportation du vecteur SPE en fichier spe.xlsx (Excel)

%__________________________________________________________________________
%% Choix du seuil de détection

centile = 95 ;
seuil = prctile(SPE, centile)
% 95% des valeurs de la SPE sont inférieures au seuil.

% Tracé de l'histogramme de la SPE et du seuil choisi
figure('Position', [figPos(1), figPos(2), width, height]) ;
histogram(SPE)
title({'Choix du seuil de détection' ; ' (données sans défaut)'})
xline(seuil, 'r', 'LineWidth', 2) ;
legend([{' histogramme de la SPE (données sans défaut)' ; sprintf(' seuil : %de centile (données sans défaut)', centile)}])

% Comparaison de la SPE avec le seuil de détection
figure('Position', [figPos(1), figPos(2), width, height]) ;
plot([1:N], SPE, 'o')
yline(seuil, 'r', 'LineWidth', 2) ;
title({'Comparaison de la SPE avec le seuil de détection' ; '(données sans défaut)'})
xlabel('Indice de l''observation')
ylabel('SPE')
legend([{' SPE', sprintf(' seuil : %de centile', centile)}])

%__________________________________________________________________________
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
title({txt1, '(données avec défaut, non centrées et non réduites)'})
xlabel('Indice de l''observation i')
ylabel(['Variable ', num2str(j)])

% On centre et on réduit avec les moyennes et écart-types des données sans
% défaut mean0 et std0 :
X2 = zeros(N2,m2) ;
for j = 1:m2
    X2(:,j) = (Xb2(:,j)-mean0(j))./(std0(j)) ;
end

clear data02 data2 j

%__________________________________________________________________________
%% Définition du noyau gaussien pour les données avec défaut (noté K12)

sigma = 3 ; % il faut jouer sur cette valeur
% Exemples :
% - lorsque sigma augmente, l diminue fortement (tout le reste constant)
% - lorsque sigma augmente, le conditionnement augmente fortement

% calcul matriciel entre la base de données avec défaut par rapport à la base de données sans défaut
K12 = exp( -(dist(X2,X').^2)./(2*sigma^2) ) ;
% calcul matriciel entre la base de données en défaut et elle-même

% K2 = exp( -(dist(X2,X2').^2)./(2*sigma^2) ) ;
% pas besoin de calculer car on n'a besoin que de K2(i,i) qui vaut 1
% ainsi, on économise beaucoup en complexité

% Bon choix de sigma ?
% mean(mean(K))
% median(median(K))

% conditionnement de K
cond_K12 = norm(inv(K12))*norm(K12)

%__________________________________________________________________________
%% Détection des observations défaillantes

SPE2 = zeros(N2,1) ; % initialisation de SPE : vecteur colonne dont chaque ligne contient la SPE de l'observation numéro i
for i = 1:N2
    SPE2(i) = 1 - K12(i,:)*Pp*LAMBDA*Pp'*K12(i,:)' ;
    % SPE2(i) = K2(i,i) - K12(i,:)*Pp*LAMBDA*Pp'*K12(i,:)' ; % mais K2(i,i)vaut toujours 1
end

clear i

% Tracé : comparaison de la SPE avec le seuil de détection
figure('Position', [figPos(1), figPos(2), width, height]) ;
ylim = 1.1 ; % valeur maximale des ordonnées que l'on affiche, prendre Inf puis 3
plot([1:N2], SPE2, 'o')
yline(seuil, 'r', 'LineWidth', 2) ;
axis([-Inf Inf -Inf ylim])
title({'Comparaison de la SPE avec le seuil de détection' ;  '(données avec défaut)'})
xlabel('Indice de l''observation')
ylabel('SPE')
legend('SPE', 'seuil de détection')

clear ylim

% Tracé de l'histogramme de la SPE et du seuil choisi
figure('Position', [figPos(1), figPos(2), width, height]) ;
histogram(SPE2)
title({'Comparaison de la SPE avec le seuil de détection' ; ' (données avec défaut)'})
xline(seuil, 'r', 'LineWidth', 2) ;
legend([{' histogramme de la SPE (données avec défaut)' ; sprintf(' seuil : %de centile (données sans défaut)', centile)}])

% defaillance = [] ;
% for i = 1:N2
%     if SPE2(i) > seuil
%         defaillance = [defaillance; i, SPE2(i), Xb2(i,6)] ;
%     end
% end
% 
% clear i

% pourcentage_de_mesures_defaillantes = round(length(defaillance)/N2*100,0) % 93

% defaillants = array2table(defaillance) ;
% head(defaillants) ;
% summary(defaillants) ;

% figure('Position', [figPos(1), figPos(2), width, height]) ;
% plot(defaillance(:,1), defaillance(:,3), 'o')
% title('Réprésentation des données défaillantes')
% xlabel('Indice de l''observation i')
% ylabel('Valeur mesurée FA')

%__________________________________________________________________________
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

j = 6 ; % tester j=2, j=6, j=7, j=8 (surtout 6, 7 et 8) (les autres sont identiques)

% Superposition des graphes (données centrées et réduites)
figure('Position', [figPos(1), figPos(2), width, height]) ;
plot([1:N], X(:,j), 'og', [1:N], X2(:,j), 'xr')
title({sprintf('Tracé de la variable %d %s', j, 'en fonction de l''index de temps') ; '(données centrées et réduites)'})
xlabel('Indice de l''observation')
ylabel(['Variable ', num2str(j)])
legend('données sans défaut', 'données avec défaut')

clear j