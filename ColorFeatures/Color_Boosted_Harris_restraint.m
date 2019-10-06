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
function [hull_msk,border_sp]=Color_Boosted_Harris_restraint(input_im,sup_image)
% [file, path]=uigetfile({...
% '*.jpg;*.bmp;*.png;*.tif','All supported format(*.jpg;*.bmp;*.png;*.tif)';...
% '*.bmp','Image format bmp (*.bmp)';...          % file selemean_seg_pixelstion with multiple
% '*.jpg','Image format jpg (*.jpg)';...
% '*.png','Image format png (*.png)';...
% '*.tif','Image format tif (*.tif)';...
% },'MultiSelemean_seg_pixelst','on');
% if ~path                                  % si rien n'est selemean_seg_pixelstionner
%         return;
% end
% input_im = double(imread([path,file])); 
sigma_g=1.5;
sigma_a=5;
nPoints=30;
si=size(input_im);
%% compute RGB-Harris Detector
% [EnIm]= ColorHarris(input_im,sigma_g,sigma_a,0.04,1);
% 
% % extract corners in total nPoints
% [x_max,y_max,corner_im,num_max]=getmaxpoints(EnIm,nPoints);
% 
% % add restraint (post processing)
% Ymsk=(y_max<round(si(1)-0.05*si(1))).*(y_max>round(0.05*si(1)));
% Xmsk=(x_max<round(si(2)-0.05*si(2))).*(x_max>round(0.05*si(2)));
% msk=logical(Xmsk.*Ymsk);
% deleteX=x_max(~msk);
% deleteY=y_max(~msk);
% x_max=x_max(msk);
% y_max=y_max(msk);
% % convex hull
% DT = delaunayTriangulation(x_max,y_max);
% k = convexHull(DT);
%  k = boundary(x_max,y_max);

% % visualize corners 
% output_im=visualize_corners(input_im,corner_im);
% 
% figure(1);imshow(uint8(input_im));
% figure(2);imshow(uint8(output_im));title('color Harris');

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
% add restraint (post processing)
Ymsk=(y_max<round(si(1)-0.05*si(1))).*(y_max>round(0.05*si(1)));
Xmsk=(x_max<round(si(2)-0.05*si(2))).*(x_max>round(0.05*si(2)));
msk=logical(Xmsk.*Ymsk);
% deleteX=x_max(~msk);
% deleteY=y_max(~msk);
x_max=x_max(msk);
y_max=y_max(msk);
% convex hull
DT = delaunayTriangulation(x_max,y_max);
if size(DT,1)==0
    luminance_im=make_image(sum(input_im,3),sum(input_im,3),sum(input_im,3));
[EnIm]= ColorHarris(luminance_im,sigma_g,sigma_a,0.04,1);

% extract corners in total nPoints
[x_max,y_max,corner_im2,num_max]=getmaxpoints(EnIm,nPoints);

% add restraint (post processing)
Ymsk=(y_max<round(si(1)-0.05*si(1))).*(y_max>round(0.05*si(1)));
Xmsk=(x_max<round(si(2)-0.05*si(2))).*(x_max>round(0.05*si(2)));
msk=logical(Xmsk.*Ymsk);
% deleteX=x_max(~msk);
% deleteY=y_max(~msk);
x_max=x_max(msk);
y_max=y_max(msk);
% convex hull
DT = delaunayTriangulation(x_max,y_max);
end
k = convexHull(DT);
%  k = boundary(x_max,y_max);

% corner_im2(deleteY,deleteX)=0;
% output_im2=visualize_corners(input_im,corner_im2);
% figure(1);
% imshow(uint8(output_im2));title('color boosted Harris');
% hold on,plot(DT.Points(k,1),DT.Points(k,2),'r');
% frame=getframe(1);
% close;
% display=frame2im(frame);

hull_msk=roipoly(rgb2gray(input_im),x_max(k),y_max(k));
if isempty(sup_image)
    border_sp=[];
else
AA2=unique(sup_image);
BB2=unique(sup_image(hull_msk));
border_sp=setdiff(AA2,BB2);
while isempty(border_sp)
    x_max(k)=[];
    y_max(k)=[];
    k = convhull(x_max,y_max);
%     k = boundary(x_max,y_max);
    hull_msk=roipoly(rgb2gray(input_im),x_max(k),y_max(k));
AA2=unique(sup_image);
BB2=unique(sup_image(hull_msk));
border_sp=setdiff(AA2,BB2);
end
end
% 
% visualize corners 
% output_im2=visualize_corners(input_im,corner_im2);
% figure(3);imshow(uint8(output_im2));title('color boosted Harris');
% hold on,plot(DT.Points(k,1),DT.Points(k,2),'r');

%% compute color luminance Harris Detector
% luminance_im=make_image(sum(input_im,3),sum(input_im,3),sum(input_im,3));
% [EnIm]= ColorHarris(luminance_im,sigma_g,sigma_a,0.04,1);
% 
% % extract corners in total nPoints
% [x_max,y_max,corner_im2,num_max]=getmaxpoints(EnIm,nPoints);
% 
% % add restraint (post processing)
% Ymsk=(y_max<round(si(1)-0.05*si(1))).*(y_max>round(0.05*si(1)));
% Xmsk=(x_max<round(si(2)-0.05*si(2))).*(x_max>round(0.05*si(2)));
% msk=logical(Xmsk.*Ymsk);
% deleteX=x_max(~msk);
% deleteY=y_max(~msk);
% x_max=x_max(msk);
% y_max=y_max(msk);
% % convex hull
% DT = delaunayTriangulation(x_max,y_max);
% k = convexHull(DT);
% %  k = boundary(x_max,y_max);
% hull_msk=roipoly(rgb2gray(input_im),x_max(k),y_max(k));
% % visualize corners 
% output_im3=visualize_corners(input_im,corner_im2);
% 
% figure(4);imshow(uint8(output_im3));title('Luminance Image');
end