function ImgNet(data,opt)
%------------------------------------------------------------------------%
% Show an image along with network visualization
% Input:
%     + data: struct data type, i.e., dataAll{3} in datSDM.mat file
%           .X [nFea nSmp]: dataset X (no need?)
%     + opt: structure
%           .gndNet: visualize groundtruth subnetworks?
%           .bkImg: background image, = 0 for clsMean,; > 0 for any image
%           .sparseImg: =1 to further plot sparse Image
% (c) 2014 Quang-Anh Dang - UCSB
%------------------------------------------------------------------------%

bkImg = 0;
if isfield(opt,'bkImg'),
    bkImg = opt.bkImg;
end

gndNet = 0;
if isfield(opt,'gndNet'),
    gndNet = opt.gndNet;
end

sparseImg = 0;
if isfield(opt,'sparseImg'),
    sparseImg = opt.sparseImg;
end


%-- find indices of edges from network W
halfW = triu(data.W);
[row col] = find(halfW>0);
x = data.syntX(:,row);
y = data.syntX(:,col);

if (bkImg==0)
    img = data.clsMean(:,2);
else
    %bkImg = 4; %previously
    img = data.Xfull(:,bkImg); %bkImg = 4 previously
end

figure;
h1 = subplot(1,2,1); imshow(reshape(img,data.r,data.c),[]); hold on
imgX = ones(data.r*data.c,1)*255; %-- initial img is white
imgX(data.idxSele) = img(data.idxSele); %-- replace with selected nodes
if sparseImg
    h2 = subplot(1,2,2);imshow(reshape(imgX,data.r,data.c),[]);
end

for i=1:size(x,2)
    subplot(h1); line([x(1,i) y(1,i)],[x(2,i)*(-1) y(2,i)*(-1)],'LineStyle','-');
    if sparseImg
        subplot(h2); line([x(1,i) y(1,i)],[x(2,i)*(-1) y(2,i)*(-1)],'LineStyle','-');
    end
end


if gndNet
    % remove non-gndFea from network W
    data.Wgnd = data.W;
    notGndFeaIdx = find(data.gndFea~=1);
    data.Wgnd(notGndFeaIdx,:) = 0;
    data.Wgnd(:,notGndFeaIdx) = 0;
%     data.Wgnd(notGndFeaIdx,notGndFeaIdx) = 0; %-- wrong! as pickup elements, not entire row/col 

    % create self-connection for isolate gndFea in network W
    diagMat = diag(data.gndFea);
    data.Wgnd = data.Wgnd + diagMat;
    
    
    halfWgnd = triu(data.Wgnd);
    [row col] = find(halfWgnd>0);
    xgnd = data.syntX(:,row);
    ygnd = data.syntX(:,col);
    
    for i=1:size(xgnd,2)
        subplot(h1); line([xgnd(1,i) ygnd(1,i)],[xgnd(2,i)*(-1) ygnd(2,i)*(-1)],...
            'Marker','.','LineStyle','-','LineWidth',1,'Color','r');
        if sparseImg
            subplot(h2); line([xgnd(1,i) ygnd(1,i)],[xgnd(2,i)*(-1) ygnd(2,i)*(-1)],...
                'Marker','.','LineStyle','-','LineWidth',1,'Color','r');
        end
    end
end

