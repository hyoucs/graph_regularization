function visualStraight(data, idxSeleDISM, figName)

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
%     h = figure('visible', 'off');
    imshow(reshape(img,data.r,data.c),[]); hold on

    for i=1:size(x,2)
        line([x(1,i) y(1,i)],[x(2,i)*(-1) y(2,i)*(-1)],'LineStyle','-');
        % line([x(1,i) y(1,i)],[x(2,i)*(-1) y(2,i)*(-1)],'Marker','.','LineStyle','-');
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
        line([xSeleDISM(1,i) ySeleDISM(1,i)],[xSeleDISM(2,i)*(-1) ySeleDISM(2,i)*(-1)], ...
            'Marker','.','LineStyle','-','LineWidth',1,'Color','r');
        % line([xSeleDISM(1,i) ySeleDISM(1,i)],[xSeleDISM(2,i)*(-1) ySeleDISM(2,i)*(-1)], ...
        %     'Marker','.','LineStyle','-','LineWidth',2,'Color','r');
        % plot(pSeleDISM(1,:),pSeleDISM(2,:)*(-1),'bo','LineWidth',2);
    end

%     % disp(figName);
    saveas(gcf, [figName, '.pdf'], 'pdf');
    saveas(gcf, [figName, '.fig'], 'fig');
    saveas(gcf, [figName, '.png'], 'png');

end