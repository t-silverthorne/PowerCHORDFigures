function Z = Out2(FUN,varargin)
% Z = Out2(FUN,VARARGIN);
%
%	Provides the second output from the function
[~,Z] = FUN(varargin{:});
end