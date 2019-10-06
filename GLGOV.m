% This code is for [1], and can only be used for non-comercial purpose. If
% you use our code, please cite [1].
% Code Author: Yijun Yan
% Email: yijun.yan@strath.ac.uk
% Date: 6/10//2019
% Note: some functions are employed from [2]

% [1] Yijun.Yan, Jinchang Ren, et al., Unsupervised image saliency
% detection with Gestalt-laws guided optimization and visual attention
% based refinement. Pattern Recognition, 2018.
% [2] Wangjiang Zhu, Shuang Liang, Yichen Wei, and Jian Sun. Saliency
% Optimization from Robust Background Detection. In CVPR, 2014.


%%
clear, clc, 
close all
srcImg = imread('test.jpg');
[h, w, ~] = size(srcImg);
%% Segment input rgb image into superpixels
pixNumInSP = 200;                           %pixels in each superpixel
spnumber = round( h * w / pixNumInSP );     %super-pixel number for current image
%     spnumber = 500;
[idxImg, adjcMatrix, pixelList] = SLIC_Split(srcImg, spnumber);
%% Get super-pixel properties
cform = makecform('srgb2lab');
lab = applycform(srcImg,cform);
[label_map,HP]=rgb2lab_quantization_JSfast(lab);
spNum = size(adjcMatrix, 1);
meanRgbCol = GetMeanColor(srcImg, pixelList);
meanLabCol = colorspace('Lab<-', double(meanRgbCol)/255);
bdIds = GetBndPatchIds(idxImg);
colDistM = GetDistanceMatrix(meanLabCol);
[clipVal, geoSigma, neiSigma] = EstimateDynamicParas(adjcMatrix, colDistM);

%% Saliency Optimization
[Rsal]=GLGO(idxImg,srcImg,spNum,label_map,HP); %% bottom up
[bgProb, bdCon, bgWeight] = EstimateBgProb(colDistM, adjcMatrix, bdIds, clipVal, geoSigma);
Final_result = Cost_function(adjcMatrix, bdIds, colDistM, neiSigma, bgWeight, Rsal); %% top down
Img = CreateImageFromSPs(Final_result, pixelList, h, w, true);