function [ map ] = gen_map( len )
%-- map indices between brain network to edge-dual network
%-- len is the number of nodes in brain network
%-- example: gen_map(4)

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

