function hist1D (C)
%-- visualize cluster distribution stored in C
%-- divide the range into max 50 intervals 

if islogical(C) %-- convert to double to enable hist
    C = double(C);
end    


if (length(unique(C))>50)
    disp('50 intervals (clusters) used as too many unique clusters...');
    step= ceil(max(C)-min(C))/50;
    x_h = min(C)-1:step:max(C)+1;     
else
    x_h = unique(C); 
end
hist(C,x_h); %show cluster histogram
[n,x] = hist(C,x_h);
barstrings = num2str(n');
text(x,n,barstrings,'horizontalalignment','center','verticalalignment','bottom')


return
