function smodel = netLassoGridSearch( data,opt )
%netlassoRunCV Summary of this function goes here
%   Detailed explanation goes here

    if (isfield(opt, 'output') && opt.output) || ~isfield(opt, 'output')
        fprintf('%6s %6s %3s %11s %6s %6s %6s %6s %6s %6s %6s\n',...
            'lmbd1','lmbd2','k',...
            'avgFea','culFea','avgAC','stdAC',...
            'PrecF','RecF','PrecU','RecU');
    end
    
    lmbd1_range = opt.lambda1s;
    lmbd2_range = opt.lambda2s;
    knn_range 	= opt.k;
    opt.verbose = 0;
    smodel = {};

    % - - - - - - - - - - - Grid search for all possible parameter combinations - - - - - - - - - - - -
    for idx_lmd2 = 1:length(lmbd2_range)
            opt.lambda2s = lmbd2_range(idx_lmd2);

        for idx_lmd1 = 1:length(lmbd1_range)
            opt.lambda1s = lmbd1_range(idx_lmd1);

            for idx_knn = 1:length(knn_range)
                opt.k = knn_range(idx_knn);

                % 
                idx_sel_DISM_ = [];
                idx_sel_DISM_fold = {};
                idx_sel_DISM_theta = {};
                idx_sel_DISM_size = zeros(1,length(unique(data.testSet)));
                accDISM = zeros(1,length(unique(data.testSet)));

                if isfield(data, 'gndFea')
                    idx_sel_DISM_prec = [];
                    idx_sel_DISM_rec = [];
                end
                idx_sel_DISM_com = [];

                % - - - - - - perform main algorithm in each fold - - - - - - 
                for idx_fold = 1:length(unique(data.testSet))

                	% - - - - - - train a model m & test on the rest - - - - - - 
                    m = DISMfoldvisual(data, opt, idx_fold);
                    
                    idx_sel_DISM_local = find(m.theta~=0)';
                    idx_sel_DISM_ = [idx_sel_DISM_, idx_sel_DISM_local];
                    idx_sel_DISM_size(idx_fold) = length(idx_sel_DISM_local);
                    idx_sel_DISM_fold{idx_fold} = idx_sel_DISM_local;
                    idx_sel_DISM_theta{idx_fold} = m.theta;
                    if idx_fold == 1
                        idx_sel_DISM_com = idx_sel_DISM_local;
                    else
                        idx_sel_DISM_com = unique(intersect(idx_sel_DISM_com, idx_sel_DISM_local));
                    end

                    if isfield(data, 'gndFea')
                        gndFeaSelenNum = sum(data.gndFea(idx_sel_DISM_local));
                        gndFeaNum = sum(data.gndFea);
                        idx_sel_DISM_prec(idx_fold) = gndFeaSelenNum / length(idx_sel_DISM_local);
                        idx_sel_DISM_rec(idx_fold) = gndFeaSelenNum / gndFeaNum;
                    end

                    accDISM(idx_fold) = mean(m.acc);

                end
                idx_sel_DISM_ = unique(idx_sel_DISM_);
                
                % union precision and recall
                if isfield(data, 'gndFea')
                        gndFeaSelenNum = sum(data.gndFea(idx_sel_DISM_));
                        gndFeaNum = sum(data.gndFea);
                        model.idx_sel_DISM_prec_union = gndFeaSelenNum / length(idx_sel_DISM_);
                        model.idx_sel_DISM_rec_union = gndFeaSelenNum / gndFeaNum;
                end

                model.lambda1 = opt.lambda1s;
                model.lambda2 = opt.lambda2s;
                model.k = opt.k;
                model.acc = accDISM;
                model.accuCV = mean(accDISM);
                model.stdCV = std(accDISM);
                model.idx_sel_DISM_ = idx_sel_DISM_;
                model.idx_sel_DISM_fold = idx_sel_DISM_fold;
                model.idx_sel_DISM_theta = idx_sel_DISM_theta;
                model.idx_sel_DISM_size = idx_sel_DISM_size;
                model.idx_sel_DISM_com = idx_sel_DISM_com;
                model.feaNum = mean(idx_sel_DISM_size);

                if isfield(data, 'gndFea')
                    model.idx_sel_DISM_prec = idx_sel_DISM_prec;
                    model.idx_sel_DISM_rec = idx_sel_DISM_rec;
                end
                smodel = [smodel model];

                if (isfield(opt, 'output') && opt.output) || ~isfield(opt, 'output')
                    fprintf('%6.3f %6.3f %d       %6.3f %d %6.3f %6.3f      %6.3f %6.3f %6.3f %6.3f\n', ...
                            opt.lambda1s, opt.lambda2s, opt.k, ...
                            mean(idx_sel_DISM_size), length(idx_sel_DISM_), mean(accDISM), std(accDISM),...
                            mean(model.idx_sel_DISM_prec), mean(model.idx_sel_DISM_rec), ...
                            model.idx_sel_DISM_prec_union, model.idx_sel_DISM_rec_union);
                end

            end
        end
    end

end

