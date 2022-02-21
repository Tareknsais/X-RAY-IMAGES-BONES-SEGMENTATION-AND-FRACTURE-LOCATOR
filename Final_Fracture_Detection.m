
clear all;
close all;
%1st CT 8th sem
%img = imread('http://i.stack.imgur.com/mHo7s.jpg');
%img = imread('frac3.jpg');
img = imread('leg_XRay.jpg');
figure(1)
imshow(img);
title('Input X Ray Image');

%% Important parameters

ImgBlurSigma = 2; % Quantité de débruitage de l'image d'entrée
MinHoughPeakDistance = 5; % Distance entre les pics dans la détection de l'angle de la transformée de Hough
HoughConvolutionLength = 40; % Longueur de la ligne à utiliser pour détecter les régions osseuses
HoughConvolutionDilate = 2; % Montant de la dilatation du noyau pour la détection des os
BreakLineTolerance = 0.25; % Tolérance pour la détection de l'extrémité de l'os
breakPointDilate = 6; % Quantité de dilatation des points finaux osseux détectés

%%%%%%%%%%%%%%%%%%%%%%%

img_gray = (rgb2gray(img)); % Load image
figure(2)
imshow(img_gray);
title('Gray Scale X Ray Image');
img_filtered = imfilter(img_gray, fspecial('gaussian', 10, ImgBlurSigma), 'symmetric'); % Débruitage
figure(3)
imshow(img_filtered);
title('image radiographique à échelle de gris débruitée');

% Détection des contours pour trouver les contours des os dans l'image
% Filtrer tout sauf les deux lignes les plus longues
% Cette caractéristique peut devoir être modifiée si la fracture n'est pas au milieu de l'os.
boneEdges = edge(img_filtered, 'canny');
figure(4)
imshow(boneEdges);
title('Edges of the bones');

boneEdges1 = bwmorph(boneEdges, 'close');
figure(5)
imshow(boneEdges1);
title('Opération morphologique sur les bords des os:edge detection ');

edgeRegs = regionprops(boneEdges1, 'Area', 'PixelIdxList');
AreaList = sort(vertcat(edgeRegs.Area), 'descend');
edgeRegs(~ismember(vertcat(edgeRegs.Area), AreaList(1:2))) = [];
edgeImg = zeros(size(img_filtered, 1), size(img_filtered,2));
edgeImg(vertcat(edgeRegs.PixelIdxList)) = 1;

%Effectuer une transformation de Hough sur l'image du bord pour trouver les angles auxquels les morceaux d'os sont trouvés.
%Utilisez la valeur maximale de la transformée de Hough par rapport à l'angle pour trouver les angles auxquels les lignes sont orientées.  S'il y a plus d'une contribution d'angle majeur, il y aura deux pics détectés, mais un seul pic s'il n'y a qu'une seule contribution d'angle majeur (c'est-à-dire que les pics ici = nombre d'os localisés = nombre de % de cassures + 1).

[H,T,R] = hough(edgeImg,'RhoResolution',1,'Theta',-90:2:89.5);
maxHough = max(H, [], 1);
HoughThresh = (max(maxHough) - min(maxHough))/2 + min(maxHough);
[~, HoughPeaks] = findpeaks(maxHough,'MINPEAKHEIGHT',HoughThresh, 'MinPeakDistance', MinHoughPeakDistance);

% Plot Hough detection results
figure(6)
plot(T, maxHough);
hold on
plot([min(T) max(T)], [HoughThresh, HoughThresh], 'g');
plot(T(HoughPeaks), maxHough(HoughPeaks), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
hold off
xlabel('Theta Value'); ylabel('Max Hough Transform');
legend({'Max Hough Transform', 'Hough Peak Threshold', 'Detected Peak'});
title('Hough Detection Plot : Max Hough transform vs Theta');



%%2nd CT 8th sem

% Locate site of break
if numel(HoughPeaks) > 1;
    BreakStack = zeros(size(img_filtered, 1), size(img_filtered, 2), numel(HoughPeaks));
    % Convolute edge image with line of detected angle from hough transform
    for m = 1:numel(HoughPeaks);

        boneKernel = strel('line', HoughConvolutionLength, T(HoughPeaks(m)));
        kern = double(bwmorph(boneKernel.getnhood(), 'dilate', HoughConvolutionDilate));
        BreakStack(:,:,m) = imfilter(edgeImg, kern).*edgeImg;
        figure()
        imshow(BreakStack(:,:,m));
        
    end

    % Prenez la différence entre les images de convolution.  L'endroit où cette différence croise le zéro (dans les limites de la tolérance) devrait être celui où se trouve la cassure.  
% Il faut filtrer nos régions ailleurs où l'os se termine simplement.
    
    
    brImg = abs(diff(BreakStack, 1, 3)) < BreakLineTolerance*max(BreakStack(:)) & edgeImg > 0;
    [BpY, BpX] = find(abs(diff(BreakStack, 1, 3)) < BreakLineTolerance*max(BreakStack(:)) & edgeImg > 0);
    brImg = bwmorph(brImg, 'dilate', breakPointDilate);
    figure(9);
    imshow(brImg);
    brReg = regionprops(brImg, 'Area', 'MajorAxisLength', 'MinorAxisLength', ...
        'Orientation', 'Centroid');
    brReg(vertcat(brReg.Area) ~= max(vertcat(brReg.Area))) = [];

    % Calculer l'ellipse de délimitation
    brReg.EllipseCoords = zeros(100, 2);
    t = linspace(0, 2*pi, 100);
    brReg.EllipseCoords(:,1) = brReg.Centroid(1) + brReg.MajorAxisLength/2*cos(t - brReg.Orientation);
    brReg.EllipseCoords(:,2) = brReg.Centroid(2) + brReg.MinorAxisLength/2*sin(t - brReg.Orientation);

else
    brReg = [];      %% Il n'y a pas de points de fracture

end

% Dessiner une ellipse autour de l'emplacement de la rupture
figure(10)
imshow(img)
hold on
colormap('gray')
if ~isempty(brReg)
    plot(brReg.EllipseCoords(:,1), brReg.EllipseCoords(:,2), 'r');
end
hold off
    
