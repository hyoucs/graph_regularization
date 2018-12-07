function visual(data, idxSeleDISM, figName)

    switch data.posetype
        case 'straight'
            visualStraight(data, idxSeleDISM, figName);
        case '3subplot'
        	visual3Subplot(data, idxSeleDISM, figName);
        case 'all'
        	visualAllPoses(data, idxSeleDISM, figName);
        case 'gnd'
            visualGnd(data, idxSeleDISM, figName);
        case 'gndsep'
            visualGndSep(data, idxSeleDISM, figName);
        otherwise
            disp('no visualization needed');
    end

end

