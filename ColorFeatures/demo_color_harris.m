% Example of color Harris detector and of color boosted Harris detector.
%
% Note that a straight forward extention from luminance to color yields only
% few changes (compare fig.2 & 4) This is caused by the fact that also in the RGB-image most
% contrast changes are in the luminance direction. Only in the case of
% iso-luminance large differences are expected.
%
% By applying color boosting the color information content of detected points is
% increased. For more information on color boosting see:
%
% LITERATURE:
% J. van de Weijer, Th. Gevers, J-M Geusebroek
% "Boosting Color Saliency in Image Feature Detection"
% IEEE Trans. Pattern Analysis and Machine Intelligence,
% vol. 27 (4), April 2005.
% clear all;
close all;
[file, path]=uigetfile({...
'*.jpg;*.bmp;*.png;*.tif','All supported format(*.jpg;*.bmp;*.png;*.tif)';...
'*.bmp','Image format bmp (*.bmp)';...          % file selemean_seg_pixelstion with multiple
'*.jpg','Image format jpg (*.jpg)';...
'*.png','Image format png (*.png)';...
'*.tif','Image format tif (*.tif)';...
},'MultiSelemean_seg_pixelst','on');
if ~path                                  % si rien n'est selemean_seg_pixelstionner
        return;
end
input_im = double(imread([path,file])); 
% input_im=double(imread('castle.tif'));          % corel image
% cform = makecform('srgb2lab');
% input_im = applycform(input_im,cform);
sigma_g=1.5;
sigma_a=5;
nPoints=30;
si=size(input_im);
%% compute RGB-Harris Detector
[EnIm]= ColorHarris(input_im,sigma_g,sigma_a,0.04,1);

% extract corners in total nPoints
[x_max,y_max,corner_im,num_max]=getmaxpoints(EnIm,nPoints);

% visualize corners 
output_im=visualize_corners(input_im,corner_im);

figure(1);imshow(uint8(input_im));
figure(2);imshow(uint8(output_im));title('color Harris');

%% compute color boosted Harris Detector
% compute boosting matrix
% here the matrix is based on a single image 
% (for retrieval the Boosting matrix could also be pre-computed on a data set)
Mboost = BoostMatrix(input_im);

% apply matrix to image
boost_im= BoostImage(input_im,Mboost);

% Optional check if boosting matrix is identity matrix after color boosting
% Mboost2 = BoostMatrix(boost_im)

% compute Harris Energy
[EnIm]= ColorHarris(boost_im,sigma_g,sigma_a,0.04,1);

% extract corners in total nPoints
[x_max,y_max,corner_im2,num_max]=getmaxpoints(EnIm,nPoints);

% % add restraint (post processing)
% Ymsk=(y_max<round(si(1)-0.05*si(1))).*(y_max>round(0.05*si(1)));
% Xmsk=(x_max<round(si(2)-0.05*si(2))).*(x_max>round(0.05*si(2)));
% msk=logical(Xmsk.*Ymsk);
% deleteX=x_max(~msk);
% deleteY=y_max(~msk);
% x_max=x_max(msk);
% y_max=y_max(msk);
% convex hull
DT = delaunayTriangulation(x_max,y_max);
k = convexHull(DT);
%  k = boundary(x_max,y_max);
% visualize corners 
% corner_im2(deleteY,deleteX)=0;
output_im2=visualize_corners(input_im,corner_im2);
figure(3);imshow(uint8(output_im2));title('color boosted Harris');
hold on,plot(DT.Points(k,1),DT.Points(k,2),'r','LineWidth',2);
%% compute color luminance Harris Detector
luminance_im=make_image(sum(input_im,3),sum(input_im,3),sum(input_im,3));
[EnIm]= ColorHarris(luminance_im,sigma_g,sigma_a,0.04,1);

% extract corners in total nPoints
[x_max,y_max,corner_im3,num_max]=getmaxpoints(EnIm,nPoints);

% visualize corners 
output_im3=visualize_corners(input_im,corner_im3);

figure(4);imshow(uint8(output_im3));title('Luminance Image');