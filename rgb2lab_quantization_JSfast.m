function [label_map,HP]=rgb2lab_quantization_JSfast(lab)
% Code Author: Yijun Yan
% Email: yijun.yan@strath.ac.uk
t1=clock;
[aa, bb, cc]=size(lab);
lab_matrix=double(reshape(lab,[aa*bb,cc]));
l = double(lab(:,:,1)); 
a = double(lab(:,:,2)); 
b = double(lab(:,:,3)); 
numOfbins=12;
binranges_l = (min(l(:)):(max(l(:))-min(l(:)))/numOfbins:max(l(:)));
binranges_a = (min(a(:)):(max(a(:))-min(a(:)))/numOfbins:max(a(:)));
binranges_b = (min(b(:)):(max(b(:))-min(b(:)))/numOfbins:max(b(:)));

if isempty(binranges_l) | isnan(binranges_l)
    index_l=0;
else
binranges_l(1)=[];
binranges_l(end)=[];
index_l = quantiz(l(:),binranges_l);
end
if isempty(binranges_a) | isnan(binranges_a)
        index_a=0;
else
binranges_a(1)=[];
binranges_a(end)=[];
index_a = quantiz(a(:),binranges_a);
end
if isempty(binranges_b) | isnan(binranges_b)
    index_b=0;
else
binranges_b(1)=[];
binranges_b(end)=[];
index_b = quantiz(b(:),binranges_b);
end
map=(numOfbins.^2).*index_l+numOfbins.*index_a+index_b+1;  %% label index
map=reshape(map,[aa,bb]);
[mc,ind]=hist(map,unique(map));
summc=sum(mc,2);
t2=clock;
time1=etime(t2,t1);
Sum=zeros(length(ind),3);
[~,newmap]=histc(map,unique(map));
for e=1:numel(newmap)
    Sum(newmap(e),1:3)=Sum(newmap(e),1:3)+lab_matrix(e,:);
end
t3=clock;
time2=etime(t3,t2);
Sum=Sum./repmat(summc,[1,3]);
Sum(:,4)=summc./sum(summc);
[num,ind2]=sort(Sum(:,4),'descend');  %% sort the C into descend order.
i=1;
per=0;
while per<0.95
    per=per+num(i);
    i=i+1;
end
ind_hp=ind2(1:i);
HP=Sum(ind_hp,:); % the number of high probability bins, mean value and probability of each bins
ind_lp=ind2(i+1:end);
LP=Sum(ind_lp,:);
dis=pdist2(HP(:,1:3),LP(:,1:3));
%% normal version
n=zeros(size(LP,1),1);
for iMerge=1:size(LP,1)
    [~,n(iMerge)]=min(dis(:,1));
    HP(n(iMerge),1:3)=(HP(n(iMerge),1:3).*HP(n(iMerge),4)+LP(1,1:3).*LP(1,4))/(HP(n(iMerge),4)+LP(1,4)); %% may not need
    HP(n(iMerge),4)=HP(n(iMerge),4)+LP(1,4);
    LP(1,:)=[];
    dis=pdist2(HP(:,1:3),LP(:,1:3));
%     map(map==ind_lp(iMerge))=ind_hp(n);  %% replace the label of LP with that of HP (the max value in map is the length of ind)
end
%% fast version
% [~,n]=min(dis,[],1);
% HP(n,1:3)=(HP(n,1:3).*repmat(HP(n,4),[1,3])+LP(:,1:3).*repmat(LP(:,4),[1,3]))./repmat((HP(n,4)+LP(:,4)),[1,3]); %% may not need
% HP(n,4)=HP(n,4)+LP(:,4);
%%
for idx=1:length(ind_lp)
    newmap(newmap==ind_lp(idx))=ind_hp(n(idx));
end
label_map=zeros(aa,bb);
for idx=1:length(ind_hp)
label_map(newmap==ind_hp(idx))=idx;
end
HP(isnan(HP))=0;
end
