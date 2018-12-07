function DISMauc(theta, gndFea, nfold)

	thetaAll = zeros(size(gndFea,1),1);
	for i = 1:nfold,
		hit = sortrows([abs(theta{i}), gndFea], -1);
		auc = computeAUC([(size(hit,1):-1:1)', hit(:,2)]);
		disp(['=-=-=-=-= SINGLE FOLD ',num2str(i),' =-=-=-=-=']);
		disp(sprintf('auc values: %6.4f', auc));
		thetaAll = max(abs(thetaAll), abs(theta{i}));
	end
	hit = sortrows([abs(thetaAll), gndFea], -1);
	auc = computeAUC([(size(hit,1):-1:1)', hit(:,2)]);
	disp(['=-=-=-=-= COMBINED FOLDS =-=-=-=-=']);
	disp(sprintf('auc values: %6.4f', auc));

end



function auc = computeAUC(dROC)

	smodel = {};
	model.dROC = dROC;
	model.algo = legend;
	smodel = [smodel model];
	auc = AUC_plot(smodel);

end