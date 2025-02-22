function smodel = DISMvis1( data,opt )
%DISMVIS1 Summary of this function goes here
%   Detailed explanation goes here

    if (isfield(opt, 'output') && opt.output) || ~isfield(opt, 'output')
        fprintf('%6s %6s %3s %11s %6s %6s %6s %6s %6s %6s %6s\n',...
            'lmbd1','lmbd2','k','avgFea','culFea','avgAC','stdAC','PrecF','RecF','PrecU','RecU');
    end

    smodel = {};
    lambda1Array = opt.lambda1s;
    lambda2Array = opt.lambda2s;
    knnArray = opt.k;
    opt.verbose = 0;

    for idx_lmd2 = 1:length(lambda2Array)
            opt.lambda2s = lambda2Array(idx_lmd2);

        for idx_lmd1 = 1:length(lambda1Array)
            opt.lambda1s = lambda1Array(idx_lmd1);

            for idx_knn = 1:length(knnArray)
                opt.k = knnArray(idx_knn);

                idxSeleDISM = [];
                idxSeleDISMfold = {};
                idxSeleDISMtheta = {};
                idxSeleDISMsize = zeros(1,length(unique(data.testSet)));
                accDISM = zeros(1,length(unique(data.testSet)));

                if isfield(data, 'gndFea')
                    idxSeleDISMprec = [];
                    idxSeleDISMrec = [];
                end
                idxSeleDISMcom = [];

                for idx_fold = 1:length(unique(data.testSet))
                    m = DISMfoldvisual(data, opt, idx_fold);
                    idxSeleDISMlocal = find(m.theta~=0)';
                    idxSeleDISM = [idxSeleDISM, idxSeleDISMlocal];
                    idxSeleDISMsize(idx_fold) = length(idxSeleDISMlocal);
                    idxSeleDISMfold{idx_fold} = idxSeleDISMlocal;
                    idxSeleDISMtheta{idx_fold} = m.theta;
                    if idx_fold == 1
                        idxSeleDISMcom = idxSeleDISMlocal;
                    else
                        idxSeleDISMcom = unique(intersect(idxSeleDISMcom, idxSeleDISMlocal));
                    end

                    if isfield(data, 'gndFea')
                        gndFeaSelenNum = sum(data.gndFea(idxSeleDISMlocal));
                        gndFeaNum = sum(data.gndFea);
                        idxSeleDISMprec(idx_fold) = gndFeaSelenNum / length(idxSeleDISMlocal);
                        idxSeleDISMrec(idx_fold) = gndFeaSelenNum / gndFeaNum;
                    end

                    accDISM(idx_fold) = mean(m.acc);

                    if isfield(opt, 'visfold') && opt.visfold
                        figName = sprintf('%s/fold_%d_%6.3f_%6.3f_%d_fea_%6.3f_%6.3f', ...
                            data.figName, idx_fold, opt.lambda1s, opt.lambda2s, opt.k, ...
                            length(idxSeleDISMlocal), accDISM(idx_fold));
                        visual(data, idxSeleDISMlocal, figName);
                    end
                end
                idxSeleDISM = unique(idxSeleDISM);
                
                % union precision and recall
                if isfield(data, 'gndFea')
                        gndFeaSelenNum = sum(data.gndFea(idxSeleDISM));
                        gndFeaNum = sum(data.gndFea);
                        model.idxSeleDISMprec_union = gndFeaSelenNum / length(idxSeleDISM);
                        model.idxSeleDISMrec_union = gndFeaSelenNum / gndFeaNum;
                end

                model.lambda1 = opt.lambda1s;
                model.lambda2 = opt.lambda2s;
                model.k = opt.k;
                model.acc = accDISM;
                model.accuCV = mean(accDISM);
                model.stdCV = std(accDISM);
                model.idxSeleDISM = idxSeleDISM;
                model.idxSeleDISMfold = idxSeleDISMfold;
                model.idxSeleDISMtheta = idxSeleDISMtheta;
                model.idxSeleDISMsize = idxSeleDISMsize;
                model.idxSeleDISMcom = idxSeleDISMcom;
                model.feaNum = mean(idxSeleDISMsize);

                if isfield(data, 'gndFea')
                    model.idxSeleDISMprec = idxSeleDISMprec;
                    model.idxSeleDISMrec = idxSeleDISMrec;
                end
                smodel = [smodel model];

                if (isfield(opt, 'output') && opt.output) || ~isfield(opt, 'output')
                    if isfield(data, 'gndFea')
                        fprintf('%6.3f %6.3f %d       %6.3f %d %6.3f %6.3f      %6.3f %6.3f %6.3f %6.3f\n', ...
                                opt.lambda1s, opt.lambda2s, opt.k, ...
                                mean(idxSeleDISMsize), length(idxSeleDISM), mean(accDISM), std(accDISM),...
                                mean(model.idxSeleDISMprec), mean(model.idxSeleDISMrec), ...
                                model.idxSeleDISMprec_union, model.idxSeleDISMrec_union);
                    else
                        fprintf('%6.3f %6.3f %d       %6.3f %d %6.3f %6.3f\n', ...
                                opt.lambda1s, opt.lambda2s, opt.k, ...
                                mean(idxSeleDISMsize), length(idxSeleDISM), mean(accDISM), std(accDISM));
                    end
                end

                if isfield(opt, 'visall') && opt.visall
                    figName = sprintf('%s/%6.3f_%6.3f_%d_avgFea_%6.3f_culFea_%d_%6.3f', ...
                            data.figName, opt.lambda1s, opt.lambda2s, opt.k, ...
                            mean(idxSeleDISMsize), length(idxSeleDISM), mean(accDISM));
                    visual(data, idxSeleDISM, figName);
                end


            end
        end
    end

end

