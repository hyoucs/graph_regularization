function auc = ROC(id)

	X =load(['roc/roc_data_', num2str(id),'.txt']);
    %{
	Y = X(:, 2:end);
	Y = int8(Y);
	Y(Y > 0) = 1;
	X = [X(:,1), double(Y)];
    %}
    smodel={}; model.dROC=[X(:,1), X(:,2)]; 
 	model.algo='ROC-MIR&FBF';  smodel=[smodel model];
 	figure;  
    auc = AUC_plot(smodel);

% 	auc = zeros(1,3);
% 
% 	smodel={}; model.dROC=[X(:,1), X(:,2)]; 
% 	model.algo='ROC-MIR&FBF';  smodel=[smodel model];
% 	figure;  
% 	
% 	model.dROC=[X(:,1), X(:, 3)]; 
% 	model.algo='ROC-78+genes'; smodel=[smodel model];
% 	
% 	model.dROC=[X(:,1), X(:, 4)]; 
% 	model.algo='ROC-all-categories'; smodel=[smodel model];
% 	AUC_plot(smodel(1:end)); 
% 
% 	model.dROC=[X(:,1), X(:, 5)]; 
% 	model.algo='ROC-sum-of-four-categories'; smodel=[smodel model];
% 	figure; AUC_plot(smodel); box on;

	print(gcf,'-dpng', ['roc/', num2str(id),'.png']);

end