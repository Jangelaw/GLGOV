function [Sal_map2]=GLGO(sup_image,Input_im,num_sp,label_map,HP)
[aa, bb, cc]=size(Input_im);
stats=regionprops(sup_image,'Area','Centroid'); 
Centroid=cat(1,stats(:).Centroid);
Area=cat(1,stats(:).Area);
lab_matrix=HP(label_map(:),1:3);
temp_c=zeros(num_sp,num_sp);
mean_SP=zeros(num_sp,3);
spLabel=reshape(sup_image,[aa*bb,1]);
H_lab=cell(num_sp,1);
for idx=1:num_sp
    H_basic=zeros(size(HP,1),1);
    msk = spLabel ==idx;
    mean_SP(idx,:)=mean(lab_matrix(msk,:),1); 
    Sp_hist=label_map(msk);
    [num,~]=unique(Sp_hist);
    [x,y]=histc(Sp_hist,num);
    H_basic(num)=x./sum(x);
    H_lab{idx} = H_basic;
end
%% SP spatial distance
D_aa=min(abs(Centroid(:,2)-1),abs(Centroid(:,2)-aa));
D_bb=min(abs(Centroid(:,1)-1),abs(Centroid(:,1)-bb));
Ds=min(D_aa,D_bb);

%% sim_d sim_c
for iSP=1:num_sp-1
    for jSP=(iSP+1):num_sp
    temp_c(iSP,jSP) = sum((1-abs(H_lab{iSP}-H_lab{jSP})).*min(H_lab{iSP},H_lab{jSP}));
    end
end 
sim_d=pdist2(Centroid,Centroid);
sim_d= (sim_d-max(sim_d(:)))/(min(sim_d(:))-max(sim_d(:)));
sim_c=temp_c+eye(num_sp,num_sp)+temp_c';
sim=sim_d.*sim_c;
%% calculate the NGC,NSS,VS_i
A=Area./sum(Area);
Area_mat=repmat(A,1,num_sp)';
W=sim_d.*Area_mat;     %% weight
dist=pdist2(mean_SP,mean_SP);
GC_j=W.*dist;  %% global contrast
GC_i=sum(GC_j,2);
SS_i= sim*(Ds)./sum(sim,2); 
NGC_i= (GC_i-min(GC_i))/(max(GC_i)-min(GC_i));  %% normalized global color contrast measure
NSS_i= (SS_i-min(SS_i))/(max(SS_i)-min(SS_i));
RGC=sim*NGC_i./sum(sim,2);
RSS=sim*NSS_i./sum(sim,2);
sal=RGC.*RSS;
sal=(sal-min(sal))/(max(sal)-min(sal));

%% color quantization
[index1,O1,~,border_sp] = MPEG7_LAB_CBH(lab_matrix,Area,mean_SP,Input_im,sup_image);
center=zeros(2*size(O1,1),2);
%% first color map
Color_map1=zeros(aa,bb);
for reidx1=1:length(index1)
Color_map1(sup_image==reidx1)=index1(reidx1);
end
%%
num_quant=size(O1,1);
L=cell(num_quant,1);
Color_label=zeros(aa,bb);
for k=1:num_quant
    temp=Color_map1==k;
    L{k}=bwlabel(temp);
    Color_label=Color_label+(L{k}+length(unique(Color_label))-1).*(L{k}~=0);
end

%% calculate the border population for each color region
label_cl=zeros(num_sp,1);
for idx=1:num_sp
    msk = sup_image ==idx;
    label_cl(idx)=unique(Color_label(msk));
end
%% background supression
border_region=label_cl(border_sp);
border_label=zeros(num_sp,1);
border_label(border_sp)=label_cl(border_sp);
[sum_cl,idx_cl]=hist(border_region,unique(border_region));
BorderProbability=zeros(num_sp,1);
for iBP=1:length(idx_cl)
num_border_cl=label_cl==idx_cl(iBP);
BorderProbability(num_border_cl)=sum_cl(iBP)/sum( num_border_cl);

end
%%
BP=(BorderProbability-min(BorderProbability))/(max(BorderProbability)-min(BorderProbability));
ObjectProbability=1-BP;
%% VSS
index_OP=index1.*(ObjectProbability==1);
index_BP=index1.*(ObjectProbability~=1);
OP_group=cell(num_quant,1);
BP_group=cell(num_quant,1);
group_label=zeros(num_sp,1);
for k=1:num_quant
    msk_op=index_OP==k;
    msk_bp=index_BP==k;
    OP_group{k}=msk_op;
    BP_group{k}=msk_bp;
    group_label(msk_op)=k;
    group_label(msk_bp)=k+num_quant;
end
Regroup=cat(1,OP_group(:),BP_group(:));
Var=zeros(length(Regroup),1);
for k=1:length(Regroup)
    R=sum(Regroup{k});
    temp_C=Centroid(Regroup{k},:);
    temp_A=Area(Regroup{k});
    Aloc=temp_A/sum(temp_A);
    center(k,:)=(temp_C'*Aloc)';
    Var(k)=sum(pdist2(temp_C,center(k,:)).*temp_A)/R;  %% temp_A / A?
end
num_labels=length(unique(Color_label));
Var_map=Var(group_label);
VS_i=Var_map;
VSS_i= (VS_i-max(VS_i))/(min(VS_i)-max(VS_i));
ROP=sim*ObjectProbability./sum(sim,2);
NOP=(ROP-min(ROP))/(max(ROP)-min(ROP));
VSS=sim*VSS_i./sum(sim,2);
%% GP based color space smooth
sal_smooth=CS_smooth(num_sp,sim,sal,Area);
Nsal=(sal_smooth-min(sal_smooth))/(max(sal_smooth)-min(sal_smooth));
Nvar=(VSS-min(VSS))/(max(VSS)-min(VSS));
Sal=Nsal.*Nvar.*NOP; %% the saliency measure for each superpixel 
Sal2=sim*Sal./sum(sim,2);
Sal_map2=(Sal2-min(Sal2))/(max(Sal2)-min(Sal2));