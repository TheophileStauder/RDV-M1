function my_tests()
% calcul des descripteurs de Fourier de la base de données
img_db_path = './db/';
img_db_list = glob([img_db_path, '*.gif']);
img_db = cell(1);
label_db = cell(1);
fd_db = cell(1);
for im = 1:numel(img_db_list);
    img_db{im} = logical(imread(img_db_list{im}));
    label_db{im} = get_label(img_db_list{im});
    disp(label_db{im}); 
    [fd_db{im},~,~,~] = compute_fd(img_db{im});
end

% importation des images de requête dans une liste
img_path = './dbq/';
img_list = glob([img_path, '*.gif']);
t=tic()

% pour chaque image de la liste...
for im = 1:numel(img_list)
   
    % calcul du descripteur de Fourier de l'image
    img = logical(imread(img_list{im}));
    [fd,r,m,poly] = compute_fd(img);
       
    % calcul et tri des scores de distance aux descripteurs de la base
    for i = 1:length(fd_db)
        scores(i) = norm(fd-fd_db{i});
    end
    [scores, I] = sort(scores);
       
    % affichage des résultats    
    close all;
    figure(1);
    top = 5; % taille du top-rank affiché
    subplot(2,top,1);
    imshow(img); hold on;
    plot(m(1),m(2),'+b'); % affichage du barycentre
    plot(poly(:,1),poly(:,2),'v-g','MarkerSize',1,'LineWidth',1); % affichage du contour calculé
    subplot(2,top,2:top);
    plot(r); % affichage du profil de forme
    for i = 1:top
        subplot(2,top,top+i);
        imshow(img_db{I(i)}); % affichage des top plus proches images
    end
    drawnow();
    waitforbuttonpress();
end
end

function [fd,r,m,poly] = compute_fd(img)
% Ici c'est le nombre de points du contour
N = 110;
M = 110;
h = size(img,1);
w = size(img,2);

% x sont les lignes 
% y sont toutes les colonnes
[x,y]=find(img);
% mx est la somme de toutes les colonnes avec des pixel blanc
mx=round(mean(y));
% my est la somme de toutes les lignes avec des pixel blanc
my=round(mean(x));
% m est le barycentre pixel blanc
m = [mx my];


% Initilisation
t = linspace(0,2*pi,N);
R = min(h,w)/2;
poly = zeros(N, 2);
r = zeros(1,N);

%On va parcourir chaque angle 
for i = 1:N
    tmp=1;
    tmpMx=mx;
    tmpMy=my;
    % Tant quand on est pas arrivé au bord de l'img
    % on cherchera le point qui se trouve sur la droite de l'angle au bord de l'img
    while ((tmpMx>1 & tmpMy>1) & (tmpMx<w & tmpMy<h))
        tmpMx=round(mx+tmp*cos(t(1,i)));
        tmpMy=round(my+tmp*sin(t(1,i)));
        %bord de l'img
        bordX=tmpMx;
        bordY=tmpMy;
        tmp=tmp+1;    
    end
    % Calcule de la distance au barycentre en partant du bord de l'img
    for tmp = 1:tmp
    tmpMx=round(bordX-tmp*cos(t(1,i)));
    tmpMy=round(bordY-tmp*sin(t(1,i)));
    % On s'arrete si on trouve un pixel blanc et on met le contour dans le
    % polygone 
    if (img(tmpMy,tmpMx) == 1)
        x=tmpMx-mx;
        y=tmpMy-my;
        r(1,i)=sqrt(x^2+y^2);
        poly(i,1)=tmpMx;
        poly(i,2)=tmpMy;
        break;
    end
    x=bordX-mx;
    y=bordY-my;
    r(1,i)=sqrt(x^2+y^2);
    poly(i,1)=bordX;
    poly(i,2)=bordY;
    end
end
%fd descripteu  de Fourier 
fd = zeros(1,N);
%R la transformé de Fourier 
R = zeros(1, M) ;
for i = 1 : M
    % On a R la transformation de Fourier de r
    R(1,i) = fft(r(1,i));
    % Le vecteur fd formé par les M premiers coefficients de R(f)/R(1)
    fd(1, i) = abs(R(1, i))/ abs(R(1, 1)) ;
end
end
