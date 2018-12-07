function visualAllPoses(data, idxSeleDISM, figName)

    Whalf = triu(data.W);
    [row, col] = find(Whalf>0);

    x = data.syntX(:,row);
    y = data.syntX(:,col);

    if isfield(data, 'Xfullnonoise')
        img1 = data.Xfullnonoise(:,2);
        img2 = data.Xfullnonoise(:,12);
        img3 = data.Xfullnonoise(:,20);
        img4 = data.Xfullnonoise(:,28);
    else
        img1 = data.Xfull(:,2);
        img2 = data.Xfull(:,12);
        img3 = data.Xfull(:,20);
        img4 = data.Xfull(:,28);
    end


    close all;
    % h = figure('visible', 'off');
    h = figure;
    h1 = subplot(2,2,1);
    imshow(reshape(img1,data.r,data.c),[]);     % left turn
    h2 = subplot(2,2,2);
    imshow(reshape(img2,data.r,data.c),[]);    % right turn
    h3 = subplot(2,2,3);
    imshow(reshape(img3,data.r,data.c),[]);    % straight turn
    h4 = subplot(2,2,4);
    imshow(reshape(img4,data.r,data.c),[]);    % up turn

    for i=1:size(x,2)
        subplot(h1);
        line([x(1,i) y(1,i)],[x(2,i)*(-1) y(2,i)*(-1)],'LineStyle','-');
        subplot(h2);
        line([x(1,i) y(1,i)],[x(2,i)*(-1) y(2,i)*(-1)],'LineStyle','-');
        subplot(h3);
        line([x(1,i) y(1,i)],[x(2,i)*(-1) y(2,i)*(-1)],'LineStyle','-');
        subplot(h4);
        line([x(1,i) y(1,i)],[x(2,i)*(-1) y(2,i)*(-1)],'LineStyle','-');
    end

    idxUnseleDISM = setdiff(1:size(data.X,1), idxSeleDISM);
    % spy(Whalf(idxSeleDISM,idxSeleDISM));

    WhalfSele = Whalf;
    WhalfSele(idxUnseleDISM, :) = 0;
    WhalfSele(:, idxUnseleDISM) = 0;
    [rowSele, colSele] = find(WhalfSele>0);

    xSeleDISM = data.syntX(:, rowSele);
    ySeleDISM = data.syntX(:, colSele);
    pSeleDISM = data.syntX(:, idxSeleDISM);

    for i=1:size(xSeleDISM,2)
        subplot(h1);
        line([xSeleDISM(1,i) ySeleDISM(1,i)],[xSeleDISM(2,i)*(-1) ySeleDISM(2,i)*(-1)], ...
            'Marker','.','LineStyle','-','LineWidth',1,'Color','r');
        subplot(h2);
        line([xSeleDISM(1,i) ySeleDISM(1,i)],[xSeleDISM(2,i)*(-1) ySeleDISM(2,i)*(-1)], ...
            'Marker','.','LineStyle','-','LineWidth',1,'Color','r');
        subplot(h3);
        line([xSeleDISM(1,i) ySeleDISM(1,i)],[xSeleDISM(2,i)*(-1) ySeleDISM(2,i)*(-1)], ...
            'Marker','.','LineStyle','-','LineWidth',1,'Color','r');
        subplot(h4);
        line([xSeleDISM(1,i) ySeleDISM(1,i)],[xSeleDISM(2,i)*(-1) ySeleDISM(2,i)*(-1)], ...
            'Marker','.','LineStyle','-','LineWidth',1,'Color','r');
        % plot(pSeleDISM(1,:),pSeleDISM(2,:)*(-1),'bo','LineWidth',2);
    end

    % disp(figName);
    saveas(gcf, [figName, '.pdf'], 'pdf');
    saveas(gcf, [figName, '.fig'], 'fig');
    saveas(gcf, [figName, '.png'], 'png');
end
