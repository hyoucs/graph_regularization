function [ line_graph ] = genLineGraph(n)

%------------------------------------------------------------------------%
% Generate the line graph G of n nodes, which means 
% G(e(i,j),e(k,l)) = 1, if e(i,j) and e(k,l) are adjacent.
% Example: for n = 4, 
% then the adjacent matrix of line graph G is
%         (a,b)  (a,c) (a,d) (b,c) (b,d) (c,d)
% (a,b)     0     1     1     1     1     0
% (a,c)     1     0     1     1     0     1
% (a,d)     1     1     0     0     1     1
% (b,c)     1     1     0     0     1     1
% (b,d)     1     0     1     1     0     1
% (c,d)     0     1     1     1     1     0
%------------------------------------------------------------------------%


    % generate edge-index matrix
    map = gen_map(n);

    % generate 
    num_fea    = n*(n-1)/2;
    line_graph = zeros(num_fea);
    % edges on same column
    for i=1:n
        for j=1:n
            if j==i
                continue;
            end
            for k=j+1:n
                if k==i
                    continue;
                end
                line_graph(map(i,j),map(i,k)) = 1;
            end
        end
    end

    line_graph = max(line_graph, line_graph');

end



function [ map ] = gen_map( len )
%-- map indices between brain line_graphwork to edge-dual line_graphwork
%-- len is the number of nodes in brain line_graphwork
%-- example: gen_map(4) will create following edge-index matrix
%   A =
%        0     1     2     3
%        1     0     4     5
%        2     4     0     6
%        3     5     6     0

    num_fea = len*(len-1)/2;
    map = zeros(len);
    cnt = 0;
    for i=1:len-1
        for j=i+1:len
            cnt = cnt+1;
            map(j,i)=cnt;
        end
    end
            
    map = max(map,map');

end

