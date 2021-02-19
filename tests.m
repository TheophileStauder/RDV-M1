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
N = 80;
M = 80;
h = size(img,1);
w = size(img,2);

[x,y]=find(img);
mx=round(mean(y));
my=round(mean(x));
m = [mx my];


% Initilisation
t = linspace(0,2*pi,N);
R = min(h,w)/2;
poly = zeros(N, 2);
r = zeros(1,N);

%On va parcourir chaque angle 
for i = 1:N
    
    tmpMx = mx;
    tmpMy = my;
    % On cherche le premier pixel blanc de la droite de cet angle
    while ((tmpMx>1 & tmpMy>1) &(tmpMx<w & tmpMy<h))
       % On s'arrete si on trouve un pixel blanc et on met le contours de l'img dans tmpMx et tmpMy
       if (img(tmpMy,tmpMx) == 1)
           dx=tmpMx-mx;
           dy=tmpMy-my;
           r(1,i)=sqrt(dx^2+dy^2);
           poly(i,1)=tmpMx;
           poly(i,2)=tmpMy;
           
       end
       
       floatX = tmpMx+i*cos(t(1,i));
       floatY = tmpMy+i*sin(t(1,i));
       tmpMx=round(floatX);
       tmpMy=round(floatY);
    end

end
fd = zeros(1,N);

for i = 1 : M
    % On a R la transformation de Fourier de r
    R(1,i) = fft(r(1,i));
    % Le vecteur fd formé par les M premiers coefficients de R(f)/R(1)
    fd(1,i) = R(1,i) / R(1,1);
end
end
