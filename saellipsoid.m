function sa = saellipsoid(varargin) 


    switch nargin
        case 1
            
            r = abs(varargin{1});
            if size(r,2) == 1
                sa = 4*pi*r^2;
            elseif size(r,2) == 3
                a = r(1);
                b = r(2);
                c = r(3);
                sa = 4*pi*(((a*b)^1.6+(a*c)^1.6+(b*c)^1.6)/3)^(1/1.6);
            end
        case 2
            a = abs(varargin{1});
            b = abs(varargin{2});
            sa = 4*pi*(((a*b)^1.6+(a*b)^1.6+(b*b)^1.6)/3)^(1/1.6);

        case 3
            a = abs(varargin{1});
            b = abs(varargin{2});
            c = abs(varargin{3});
            sa = 4*pi*(((a*b)^1.6+(a*c)^1.6+(b*c)^1.6)/3)^(1/1.6);
    end
end