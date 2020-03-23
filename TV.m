%% Projet - inpainting TV
% �cole centrale de Lyon
% Auteur : Bastien Laville

% Inpainting 3 couleurs - importation et d�tection masque

% Ajout des dossiers de travail
addpath('images_alterees') 

title_file = 'broken';
type_file = '.png';

I = imread([title_file type_file]);

[Nx, Ny, ~] = size(I);
% N = numel(I(:,:,1));

% Bruitage de l'image
sigma = 0; % Niveau de bruit � injecter dans toute l'image (E = D^c)

% Pour tester la reconstruction + d�bruitage
% prendre title_file='oiseau'
% lambda_smooth = 8, a = 5e1 et Nit = 200

% Canaux
ImR = I(:,:,1);
ImRd = double(ImR) + sigma*rand(Nx,Ny);
ImV = I(:,:,2);
ImVd = double(ImV) + sigma*rand(Nx,Ny);
ImB = I(:,:,3);
ImBd = double(ImB) + sigma*rand(Nx,Ny);

% Image alt�r�
ImAltere = zeros(Nx,Ny,3,'uint8');
ImAltere(:,:,1) = uint8(ImRd);
ImAltere(:,:,2) = uint8(ImVd);
ImAltere(:,:,3) = uint8(ImBd);

% D�tection des parties alt�r�es par seuillage
limitRouge = 251;
gamma = double((ImV >= limitRouge)); % Domaine gamma
gammat = double((ImV < limitRouge)); % Domaine omega priv� de gamma


figure('Name', "Image alt�r�e et masque des portions � reconstruire")
subplot(1,2,1)
imagesc(ImAltere)
axis('on', 'image');
subplot(1,2,2)
imshow(gamma)

%% Reconstruction 3 canaux

E = double(gammat); % en g�n�ral les parties bruit�es sont en fait toute
                    % l'image.
lambda_smooth = 8e2; % multiplicateur de lagrange (contrainte sur bruit).
a = 5e1; % Terme pour non-annulation du gradient.

% a plus grand acc�l�re la convergence, mais rapproche la reconstruction
% TV de la reconstruction harmonique

% IMPORTANT : comme l'image est encod�e en 8 bits, on travaille avec une 
% image u qui prend des valeurs entre 0 et 255, l� o� l'article utilise u 
% prenant ses valeurs entre 0 � 1.
% Dans nos programmes, a sera de l'ordre 1e0 � 1e2

% Initialisation
u_reconstruitR = tv_o(ImRd,a,lambda_smooth,E);
u_reconstruitV = tv_o(ImVd,a,lambda_smooth,E);
u_reconstruitB = tv_o(ImBd,a,lambda_smooth,E);


% Affichage en simultan�
figure('Name', "Reconstruction de l'image")
subplot(1,2,1)
imagesc(ImAltere)
axis('on', 'image');
title('Image alt�r�e', 'FontSize', 18)


% Application du sch�ma TV
Nit = 4e1;
coutVecteurR = zeros(1,Nit);
coutVecteurV = zeros(1,Nit);
coutVecteurB = zeros(1,Nit);
for i=1:Nit
    [ReconstruiR,coutR] = tv_o(u_reconstruitR,a,lambda_smooth,E);
    u_reconstruitR = ReconstruiR;
    [ReconstruiV,coutV] = tv_o(u_reconstruitV,a,lambda_smooth,E);
    u_reconstruitV = ReconstruiV;
    [ReconstruiB,coutB] = tv_o(u_reconstruitB,a,lambda_smooth,E);
    u_reconstruitB = ReconstruiB;
    
    % Cr�ation de l'image obtenue par application d'une it�ration de l'algo
    ImTV = zeros(Nx,Ny,3,'uint8');
    ImTV(:,:,1) = uint8(u_reconstruitR);
    ImTV(:,:,2) = uint8(u_reconstruitV);
    ImTV(:,:,3) = uint8(u_reconstruitB);
    
    % Affichage
    subplot(1,2,2)
    imagesc(ImTV)
    axis('on', 'image');
    title(['It�ration : ', num2str(i)], 'FontSize', 18)
    drawnow;
    
    % Mise � jour du co�t
    coutVecteurR(i) = coutR;
    coutVecteurV(i) = coutV;
    coutVecteurB(i) = coutB;
end


% Sauvegardons nos images
imwrite(ImTV, ['images_reparees/' title_file '_reparee.png'], 'png')
imwrite(ImAltere, ['images_reparees/' title_file '.png'], 'png')


% Affichage de l'�nergie de chaque canal de couleur
figure('Name', "�nergie de l'image")
semilogy(coutVecteurR,'-or','LineWidth',1.5)
hold on
semilogy(coutVecteurV,'-og','LineWidth',1.5)
hold on
semilogy(coutVecteurB,'-ob','LineWidth',1.5)
legend('Rouge','Vert','Bleu')
grid on
xlabel("Nombre d'it�rations",'FontSize',16)
ylabel('$\mathcal{J}_\lambda^a(\mathbf{u})$',...
        'FontSize',20,'Interpreter','latex')
title(sprintf('�volution �nergie par canal (a = %ld, \\lambda = %.0g)',...
                a, lambda_smooth),'FontSize',16)
saveas(gcf,['fig/' title_file '-energie.png'])




% %% Canny detection
% 
% subplot(2, 2, 1);
% imshow(I)
% axis('on', 'image');
% title('Original Image')
% 
% % Convert to gray scale.
% grayImage = rgb2gray(I);
% subplot(2, 2, 2);
% imshow(grayImage)
% axis('on', 'image');
% title('Grey Scale Image')
% 
% % Get edges
% Canny_img = edge(grayImage, 'Canny');
% subplot(2, 2, 3);
% imshow(Canny_img, [])
% axis('on', 'image');
% title('Edge Detected Image')
% 
% % Enlarge figure to full screen.
% % set(gcf, 'Units', 'Normalized', 'Outerposition', [0, 0.05, 1, 0.95]);
% 
% 
% % %% Bruits
% % 
% % Nsinc = 4e4;
% % X = linspace(-10,10,Nsinc);
% % Y = sinc(X);
% % Ybruite = Y + 0.1*rand(1,Nsinc);
% % 
% % 
% % subplot(1, 2, 1);
% % plot(X,Ybruite)
% % axis([-10 10 -0.4 1.2])
% % grid on
% % 
% % Yreconstruit = tvd_mm(Ybruite,0.7,1000);
% % 
% % subplot(1, 2, 2);
% % plot(X,Yreconstruit)
% % axis([-10 10 -0.4 1.2])
% % grid on

