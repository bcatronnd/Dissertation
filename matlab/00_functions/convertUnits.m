function varargout = convertUnits(from2,varargin)

switch(from2)
    case 'inch2meter'
        scale = 0.0254;
end

for aa=1:length(varargin)
    varargout{aa} = varargin{aa}*scale;
end
end