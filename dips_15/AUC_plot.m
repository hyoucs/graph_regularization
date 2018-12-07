function auc_value = AUC_plot(modelArray)
%-- Plot a ROC curve for each model m in smodel 
%-- m.dROC =[col1 col2] where:
%      + col1: ranked subjects (real numbers) by algo; 
%      + col2: ground truth (0/1 numbers) according to that rank
%  http://www.mathworks.se/help/matlab/ref/linespec.html   

patt_id={'-.k*';':bs';'--ro';'-.mx';':kd';'-gx';'--k';'-r';'-b';'-.r*';'--mo';':bs';'-.k*';':bs';'--ro';'-.mx';':kd';'-gx';'--k';'-r';'-b';'-.r*';'--mo';':bs'};    
cstr={};
for i=1:length(modelArray)
    model = modelArray{i}; 
    pattForm = patt_id{i};
    strLegend = '';

    [auc] = prec_rec(model.dROC(:,1),model.dROC(:,2),'numThresh',25,'plotBaseline',0,'style',pattForm);
    auc_value = auc;  
    
    if isfield(model,'algo')
        strLegend= strcat(strLegend,model.algo);
    end
    
%     if isfield(model,'k')
%         strLegend= strcat(strLegend,' k= ',num2str(model.k,'%3d'));
%     end
%     
%     if isfield(model,'alpha')
%         strLegend= strcat(strLegend,' alp= ',num2str(model.alpha,'%5.1f'));
%     end
   
    strLegend= strcat(strLegend,' auc= ',num2str(auc,'%.3f'));
    cstr=[cstr strLegend];
end
legend(cstr,'Location','SouthEast'); grid on;

end

