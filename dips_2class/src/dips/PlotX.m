function PlotX(X,c,pos,ttxt,lbldisp,fdsp)
%------------------------------------------------------------------------%
% Plot data points in a specific screen's position
% Input:
%     + X: Data set in [nFea nSmp]
%     + c: cluster label (diff colors for diff clusters), 
%          if [] then all same colors!
%     + pos: position in col. with 3 figs each
%     + ttxt: title text to display on the figure:
%           - use [] for null text
%     + lbldisp: display a label-number next to a point:
%           - use [] for no label
%           - use 0 for label all instances
%           - use i>0 for label points in i-cluster 
%     + fdsp: text feature to display
%           - if provided (e.g., disp outlier factor)
%           - otherwise, display point id as label
%   Example:
%        PlotX(X,c,1,'Label points in 2nd cluster',2,'');
%------------------------------------------------------------------------%

%-- set figure position based on pos param
%-- scrsz:[left, bottom, width, height] with (1,1) is the screen's left,bot
scrsz = get(0,'ScreenSize'); 
nfigperrow=4;
nfigpercol=3;
totfig=nfigperrow*nfigpercol; 
figW=(scrsz(3)-0.1*scrsz(3))/nfigperrow;       %-- div 4 for 4 figures per row
figH=(scrsz(4)-0.1*scrsz(4))/nfigpercol;       %-- div 4 for 3 figures per col
%pattSet=['k.';'r*';'b+';'k+';'yo';'rx';'b*';'g*'];

%-- gray,red,blue,black,magenta,cyan,green
colSet={[.2 .2 .2],[1 0 0], [0 0 1], [0 0 0], [1 0 1], [0 1 1],[0 1 0]}; 
pattSet={'+','o','*','.','s','d','p','h'};
% c=COLORS(mod(i-1,length(COLORS))+1);

if isempty(pos), pos=0; end

if (nargin<6), fdsp=''; end

posArray=[];
for i=1:totfig
    if(mod(i,nfigpercol)~=0)
            pL=mod(i,nfigpercol)-1;
    else
            pL=nfigpercol-1;
    end
    pB=ceil(i/nfigpercol)-1;
    lef_bot=[10 50]+[pB*figW pL*figH];    
    posArray=[posArray' lef_bot']';
end


if(pos~=0) %-- new figure if pos is provided
    p=posArray(mod(pos-1,totfig)+1,:);
    figure('Position',[p figW figH]); 
end

[D N]=size(X);


if(D>3), X(4:end,:)=[];  end %-- max 3Dim for visualization

if isempty(c)||(~exist('c', 'var')),  c=ones(1,N); end

%-- set figure margin
minX=min(X');
maxX=max(X');
margin=0.1;
pct=((maxX-minX)+0.1)*margin; %-- add .1 to avoid max=min
K=unique(c);
% if K>length(colSet), error('Too many clusters...'); end;

for i=1:length(K) 
    ci=find(c==K(i));
    if ~isempty(ci)
        colId=mod(i,length(colSet))+1;
        if(D==1)
           plot(X([ci]),0,pattSet{colId},'Color',colSet{colId});    
        elseif(D==2)
           plot(X(1,ci),X(2,ci),pattSet{colId},'Color',colSet{colId});
           xlabel('feature 1'); ylabel('feature 2');
        else
           plot3(X(1,ci),X(2,ci),X(3,ci),pattSet{colId},'Color',colSet{colId}); 
%            xlabel('feature 1'); ylabel('feature 2');zlabel('feature 3');
        end
        hold on; 
    end
end

% if(pos~=0)
%     if(D==2), axis([minX(1)-pct(1) maxX(1)+pct(1) minX(2)-pct(2) maxX(2)+pct(2)]); end;
%     if(D==3), axis([minX(1)-pct(1) maxX(1)+pct(1) minX(2)-pct(2) maxX(2)+pct(2) minX(3)-pct(3) maxX(3)+pct(3)]); end;   
% end




if ~isempty(lbldisp)    %-- display point labels 
    cj=1:length(X); 
    dsp=cj;
    if (lbldisp~=0), cj=find(c==lbldisp); end %-- lable all x's in lbldisp class
    if ~isempty(fdsp), dsp=fdsp(cj);end
    for i=1:length(cj)    
        offset=(maxX(1)-minX(1))* (-0.02 + 0.04.*rand(1,1));
        if(D==1)
        	text(X(1,cj(i))+offset,0,num2str(dsp(i)));
        elseif(D==2)
            text(X(1,cj(i))+offset,X(2,cj(i))+offset,num2str(dsp(i)));
        else
        	text(X(1,cj(i))+offset,X(2,cj(i))+offset,X(3,cj(i))+offset,num2str(dsp(i)));
        end
    end
end

if ~isempty(ttxt),title(ttxt); end

end
