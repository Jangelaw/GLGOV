
function [lab_labelmap,C,mask_fore,border_sp]=MPEG7_LAB_CBH(lab,Area,mean_SP,Img,sup_image)
% Code Author: Yijun Yan
% Email: yijun.yan@strath.ac.uk
Lpixel=lab(:,1);
Apixel=lab(:,2);
Bpixel=lab(:,3);
[mask_fore,border_sp]=Color_Boosted_Harris_restraint(double(Img),sup_image); 
foremap=lab(mask_fore,:);
backmap=lab(~mask_fore,:);
foremean=mean(foremap,1);
backmean=mean(backmap,1);
% forestd=std(foremap,0,1);
% backstd=std(backmap,0,1);
D=sqrt(pdist2(foremean,backmean));
% thresh=D/2;
if D<=10
thresh=5;
elseif D>10&&D<=15
    thresh=10;
elseif D>15
    thresh=15;
end

threshL=0.5*(mean(Lpixel(mask_fore))+mean(Lpixel(~mask_fore)));
threshA=0.5*(mean(Apixel(mask_fore))+mean(Apixel(~mask_fore)));
threshB=0.5*(mean(Bpixel(mask_fore))+mean(Bpixel(~mask_fore)));
L=mean_SP(:,1);
A=mean_SP(:,2);
B=mean_SP(:,3);
l=L>threshL;
a=A>threshA;
b=B>threshB;
l=double(l);
b=double(b);
a=double(a);
b(b==1)=2;
a(a==1)=4;
lab_labelmap=l+a+b+1;                          %%rgb space is splitted into 8 bins, and 1-8 are the value of each bins

Sum = zeros(8, 5);
idx=zeros(length(mean_SP),1);
for e=1:length(mean_SP)
    Sum(lab_labelmap(e),1)=Sum(lab_labelmap(e),1)+L(e);  %% compute the sum of L value in each bins
    Sum(lab_labelmap(e),2)=Sum(lab_labelmap(e),2)+A(e);  %% compute the sum of A value in each bins
    Sum(lab_labelmap(e),3)=Sum(lab_labelmap(e),3)+B(e);  %% compute the sum of B value in each bins
    Sum(lab_labelmap(e),4)=Sum(lab_labelmap(e),4)+Area(e);     %% calculate the number of pixel points in each bins
    Sum(lab_labelmap(e),5)=Sum(lab_labelmap(e),5)+1; 
end

ori_idx=find(any(Sum,2)>0);
for re_idx=1:length(ori_idx)
    lab_labelmap(lab_labelmap==ori_idx(re_idx))=re_idx;
end

Sum(any(Sum,2)==0,:)=[];
C=zeros(size(Sum,1),4);

for k=1:size(C,1)
    C(k,1:3) = Sum(k,1:3)./Sum(k,5);
    C(k,4)=Sum(k,4)./(sum(Area));
end
dis=zeros(size(C,1));
for k=1:size(C,1)
    dis(k,:)=sqrt((C(k,1)-C(:,1)).^2+(C(k,2)-C(:,2)).^2+(C(k,3)-C(:,3)).^2); %%calculate mutual distance of two adjacent C, and then
end
% dis=triu(dis);
while find((dis<thresh)&(dis>0))
    [m, n]=find((dis<thresh)&(dis>0));
    C(m(1),1:3)=C(m(1),1:3).*(C(m(1),4)./(C(m(1),4)+C(n(1),4)))+C(n(1),1:3).*(C(n(1),4)./(C(m(1),4)+C(n(1),4)));
    C(m(1),4)=C(m(1),4)+C(n(1),4);
    
    lab_labelmap(lab_labelmap==n(1))=m(1);  %% re arrange the index
    lab_labelmap(lab_labelmap>n(1))=lab_labelmap(lab_labelmap>n(1))-1;
    
    C(n(1),:)=[];
    dis=zeros(size(C,1));
    for k=1:size(C,1)
        dis(k,:)=sqrt((C(k,1)-C(:,1)).^2+(C(k,2)-C(:,2)).^2+(C(k,3)-C(:,3)).^2);
    end
    %     dis=tril(dis);
    % q=q+1;
end
%    dis(dis==0)=NaN;
% while find((C(:,4)<0.06))
%     [m, ~]=find((C(:,4)<0.06));                                              %%if the color cluster take 1% in the whole RGB space, combine this cluster with the nearest color cluster.
%     [~, y]=find(dis(m(1),:)==min(dis(m(1),:)));
%     C(y(1),1:3)=C(m(1),1:3).*(C(m(1),4)./(C(m(1),4)+C(y(1),4)))+C(y(1),1:3).*(C(y(1),4)./(C(m(1),4)+C(y(1),4)));
%     C(y(1),4)=C(m(1),4)+C(y(1),4);
% 
%     lab_labelmap(lab_labelmap==m(1))=y(1);  %% re arrange the index
%     lab_labelmap(lab_labelmap>m(1))=lab_labelmap(lab_labelmap>m(1))-1;
% 
% 
%     C(m(1),:)=[];
%     dis=zeros(size(C,1));
%     for k=1:size(C,1)
%     dis(k,:)=sqrt((C(k,1)-C(:,1)).^2+(C(k,2)-C(:,2)).^2+(C(k,3)-C(:,3)).^2); 
%     end
%     %     dis=tril(dis);
%     dis(dis==0)=NaN;
% %     q=q+1;
% end

C(:,1:3)=round(C(:,1:3));
C(:,4)=str2num(sprintf('%7.4f', C(:,4)))';
end
