function [y] = quick_combinejaclays_make_profile_spectral_trends(jac0,jaclayers,final_numlays);

indices = [];

iNum = length(jaclayers);
ichunks = ceil(iNum/final_numlays);

for ii = final_numlays : -1 : 1
  iaList = (1:ichunks) + (ii-1)*ichunks;
  oo = find(iaList >= 1 & iaList <= iNum);
  iaList = iaList(oo);
  if ii == 1
    iaList = 1 : max(iaList);
  end
  iaList = jaclayers(iaList);
  if length(iaList) > 0  
    %fprintf(1,'ii = %3i iaList(min/max/len) = %3i -> %3i / %2i \n',ii,iaList(1),iaList(end),length(iaList))
    indices{ii} = iaList;
    y(ii,:,:) = squeeze(sum(jac0(iaList,:,:),1));
  end
end
