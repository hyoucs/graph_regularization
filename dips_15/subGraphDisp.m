function subGraphDisp(W,opt)
%------------------------------------------------------------------------%
% To display subnetworks extracted from affinity matrix W
% Input:
%     + W: affinity matrix defining entire network
%     + opt: structure 
%           .topFea: indices of top selected Fea/Nodes outputed by an algorithm
%           .gndFea (optional, 0/1-vector): ground truth Feas/Nodes, color them if provided
%           .feaName (optional): used to display fea/node names if provided
% Output:
%     + None
% Example:
% (c) 2014 Quang-Anh Dang - UCSB
%------------------------------------------------------------------------%

if (length(opt.topFea) > 600)
    error('Unable to visualize due to large subnetwork');
end

W=full(W);
W=W(opt.topFea,opt.topFea);
W=triu(W,1);

if isfield(opt,'feaName'),    
    feaName = opt.feaName(opt.topFea); 
else
%     feaName= strread(num2str(opt.topFea),'%s');
    feaName=cellstr(num2str(opt.topFea))';
end

gObj = biograph(W,feaName,'ShowArrows','off') 
gObj = view(gObj);

%opt.topFea'

if isfield(opt,'gndFea'),    
    %-- color GTnodes
    gtNodeIdx= find(opt.gndFea==1);
    [tf,loc]= ismember(gtNodeIdx,opt.topFea); %-- ismember(B,A): are bi's members of A
    loc(find(loc==0))=[]; %-- excluding Nodes not as GTnodes
    set(gObj.nodes(loc),'Color',[1 0 0],'size',[40 30]);
    dolayout(gObj);
end

end