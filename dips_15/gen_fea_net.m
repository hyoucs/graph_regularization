function [ net ] = gen_fea_net( map )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

len = size(map,1);
num_fea = len*(len-1)/2;
net = zeros(num_fea);
% edges on same column
for i=1:len
    for j=1:len
        if j==i
            continue;
        end
        for k=j+1:len
            if k==i
                continue;
            end
            net(map(i,j),map(i,k)) = 1;
        end
    end
end

net = max(net, net');

end

