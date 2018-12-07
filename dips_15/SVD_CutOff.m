function [U, S, V]=SVD_CutOff(X,SVDpct)
%-- Perform SVD(X) and cutoff small singluars (or keep top singluars)
%-- SVDpct: 
%             > 1: for number of top eVals
%             =< 1: for percentage of top eVals

    
    [U, S, V]=svd(X);
%     diag(S)
    singluar_SVD = full(diag(S));
    if SVDpct > 1
        idx = SVDpct;
        if idx < length(singluar_SVD)
            U(:,idx+1:end)=[];
            V(:,idx+1:end)=[];
            S=diag(singluar_SVD(1:idx));
        end
    else
        sumEig = sum(singluar_SVD)*SVDpct;
        sumNow = 0;
        for idx = 1:length(singluar_SVD)
            sumNow = sumNow + singluar_SVD(idx);
            if sumNow >= sumEig
                break;
            end
        end
        U(:,idx+1:end)=[];
        V(:,idx+1:end)=[];
        S=diag(singluar_SVD(1:idx));
    end
end