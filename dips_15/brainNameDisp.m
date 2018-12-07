function brainNameDisp(theta,name)
%------------------------------------------------------------------------%
% display node-names in each edge provided in theta in descending order
% Input:
%     + theta [1 nFea]: 0 and non-0 coeffs of Brain's edge-dual graph
%     + name [1 nArea]: a cell stores names of nodes in the brain network
% Output: None
% (c) 2014 Quang-Anh Dang - UCSB
%------------------------------------------------------------------------%



feaSelIdx   = find(theta~=0);
% feaSelCoef  = theta(feaSelIdx);
% [feaSelIdx feaSelCoef]

[ordFeaCoef idx]  = sort(abs(theta(feaSelIdx)),'descend');
ordFeaIdx = feaSelIdx(idx);
% [ordFeaIdx ordFeaCoef] 

edgeMat = triu(gen_map(length(name)));
for i=1:length(ordFeaIdx)
    [r c] = find(edgeMat == ordFeaIdx(i));
    fprintf('%4i: %s -- %s \n',ordFeaIdx(i), name{r}, name{c});
end


end