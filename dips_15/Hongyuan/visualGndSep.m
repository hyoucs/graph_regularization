function visualGndSep(data, idxSeleDISM, figName)

    Whalf = triu(data.W);
    [row, col] = find(Whalf>0);

    x = data.syntX(:,row);
    y = data.syntX(:,col);

    if isfield(data, 'Xfullnonoise')
        img = data.Xfullnonoise(:,2);
    else
        img = data.Xfull(:,2);
    end

    close all;

    % draw the standard image and the underlying network
    h = figure;
    imshow(reshape(img,data.r,data.c),[]); hold on

    for i=1:size(x,2)
        line([x(1,i) y(1,i)],[x(2,i)*(-1) y(2,i)*(-1)],'LineStyle','-');
    end

    % draw the selected features of the method on subplot 1 and subplot 3
    idxUnseleDISM = setdiff(1:size(data.X,1), idxSeleDISM);
    WhalfSele = Whalf;
    WhalfSele(idxUnseleDISM, :) = 0;
    WhalfSele(:, idxUnseleDISM) = 0;
    [rowSele, colSele] = find(WhalfSele>0);

    xSeleDISMs = data.syntX(:, rowSele);
    ySeleDISMs = data.syntX(:, colSele);

    for i=1:size(xSeleDISMs,2)
        line([xSeleDISMs(1,i) ySeleDISMs(1,i)],[xSeleDISMs(2,i)*(-1) ySeleDISMs(2,i)*(-1)], ...
            'Marker','.','LineStyle','-','LineWidth',1,'Color','r');
    end

    % draw the groundtruth features on subplot 2 and subplot 3
    idxSeleDISM = find(data.gndFea ~= 0);
    idxUnseleDISM = setdiff(1:size(data.X,1), idxSeleDISM);
    % spy(Whalf(idxSeleDISM,idxSeleDISM));

    WhalfSele = Whalf;
    WhalfSele(idxUnseleDISM, :) = 0;
    WhalfSele(:, idxUnseleDISM) = 0;
    [rowSele, colSele] = find(WhalfSele>0);

    xSeleDISMg = data.syntX(:, rowSele);
    ySeleDISMg = data.syntX(:, colSele);

    for i=1:size(xSeleDISMg,2)
        line([xSeleDISMg(1,i) ySeleDISMg(1,i)],[xSeleDISMg(2,i)*(-1) ySeleDISMg(2,i)*(-1)], ...
            'Marker','.','LineStyle','-','LineWidth',1,'Color','y');
    end

    disp(figName);
    saveas(gcf, [figName, '.png'], 'png');
    saveas(gcf, [figName, '.fig'], 'fig');
    print('-depsc',[figName, '.eps']);




    close all;

    % draw the standard image and the underlying network
    h = figure;
    imshow(reshape(img,data.r,data.c),[]); hold on

    for i=1:size(x,2)
        line([x(1,i) y(1,i)],[x(2,i)*(-1) y(2,i)*(-1)],'LineStyle','-');
    end

    for i=1:size(xSeleDISMg,2)
        line([xSeleDISMg(1,i) ySeleDISMg(1,i)],[xSeleDISMg(2,i)*(-1) ySeleDISMg(2,i)*(-1)], ...
            'Marker','.','LineStyle','-','LineWidth',1,'Color','y');
    end

    for i=1:size(xSeleDISMs,2)
        line([xSeleDISMs(1,i) ySeleDISMs(1,i)],[xSeleDISMs(2,i)*(-1) ySeleDISMs(2,i)*(-1)], ...
            'Marker','.','LineStyle','-','LineWidth',1,'Color','r');
    end

    disp(figName);
    saveas(gcf, [figName, '2.png'], 'png');
    saveas(gcf, [figName, '2.fig'], 'fig');
    print('-depsc',[figName, '2.eps']);


end