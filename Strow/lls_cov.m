num_fixed = length(find(driver.jacobian.qstYesOrNo == 1));
if driver.jacobian.numQlays ~= 1
   disp('Error, only T and Q profiles together for now:')
end
cmat_size = num_fixed + 2*driver.jacobian.numlays;

