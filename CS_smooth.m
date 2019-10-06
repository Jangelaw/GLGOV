function sal_smooth=CS_smooth(num_sp,sim,sal,Area)
% Code Author: Yijun Yan
% Email: yijun.yan@strath.ac.uk
K=double(round(num_sp/8));
[T_matirx,order_matrix] = sort(sim,2,'descend');
T_tot=T_matirx(:,2:K);
order_tot=order_matrix(:,2:K);
T=sum(T_tot,2);
temp_sal=sal';
% temp_vop=VSS';
temp_area=Area';
Fsc_sal=(repmat(T,1,K-1)-T_tot).*temp_sal(order_tot).*temp_area(order_tot);
sal_smooth=sum(Fsc_sal,2)./(T*(K-1));