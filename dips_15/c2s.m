function out=c2s(in)
% C2S Converts cell to structure array.
%
% Synopsis:
%  out = c2s({'item1',value1,'item2',value2,...})
%
% Description:
%  out = c2s({'item1',value1,'item2',value2,...}) is a shortcut 
%   for out = struct('item1',value1,'item2',value2,...).
%
% Example:
%  out = c2s({'a',1,'b',2,'c','Hello'})
%

if iscell(in),
  out = struct(in{:});
else
  out = in;
end

return;
% EOF