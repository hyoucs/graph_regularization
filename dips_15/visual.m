function visual(data, m)

    Whalf = triu(data.W);
    [row, col] = find(Whalf>0);

    x = data.syntX(:,row);
    y = data.syntX(:,col);

    % img = data.clsMean(:,2);
    img = data.Xfull(:,2);

    h = figure;
    h1 = subplot(1,3,1); PlotX(data.syntX,'','','',''); axis square
    h2 = subplot(1,3,2); imshow(reshape(img,data.r,data.c),[]); hold on
    imgX = ones(data.r*data.c,1)*256;
    imgX(data.idxSele) = img(data.idxSele);
    % imgX(data.idxSele) = data.clsMean(data.idxSele,2);
    h3 = subplot(1,3,3);imshow(reshape(imgX,data.r,data.c),[]); hold on

    for i=1:size(x,2)
        subplot(h1);
        line([x(1,i) y(1,i)],[x(2,i) y(2,i)],'LineStyle','-');
        subplot(h2);
        line([x(1,i) y(1,i)],[x(2,i)*(-1) y(2,i)*(-1)],'Marker','.','LineStyle','-');
    end


    idxSeleDISM = find(m.theta~=0);
    idxUnseleDISM = setdiff(1:size(data.X,1), idxSeleDISM);
    % spy(Whalf(idxSeleDISM,idxSeleDISM));

    WhalfSele = Whalf;
    WhalfSele(idxUnseleDISM, :) = 0;
    WhalfSele(:, idxUnseleDISM) = 0;
    [rowSele, colSele] = find(WhalfSele>0);

    xSeleDISM = data.syntX(:, rowSele);
    ySeleDISM = data.syntX(:, colSele);
    pSeleDISM = data.syntX(:,idxSeleDISM);

    for i=1:size(xSeleDISM,2)
        subplot(h1);
        line([xSeleDISM(1,i) ySeleDISM(1,i)],[xSeleDISM(2,i) ySeleDISM(2,i)],'LineStyle','-','Color','b');
        plot(pSeleDISM(1,:),pSeleDISM(2,:),'bo','LineWidth',2);
        subplot(h2);
        line([xSeleDISM(1,i) ySeleDISM(1,i)],[xSeleDISM(2,i)*(-1) ySeleDISM(2,i)*(-1)],'Marker','.','LineStyle','-','LineWidth',2,'Color','r');
%         plot(pSeleDISM(1,:),pSeleDISM(2,:)*(-1),'bo','LineWidth',2);
        subplot(h3);
        line([xSeleDISM(1,i) ySeleDISM(1,i)],[xSeleDISM(2,i)*(-1) ySeleDISM(2,i)*(-1)],'Marker','.','LineStyle','-','LineWidth',2,'Color','r');
%         plot(pSeleDISM(1,:),pSeleDISM(2,:)*(-1),'bo','LineWidth',2);
    end
    
    

end