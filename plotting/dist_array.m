function varargout = dist_array(A,dim)
% distributes array to multiple output variables
%
% A: input array
% dim: dimension to distribute along (default: 1)

if nargin==1
    varargout = num2cell(A);
else
    varargout = num2cell(A,dim);
end

end


