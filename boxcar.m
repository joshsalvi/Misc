%boxcar is a function for smoothing 1D data by
%imposing a sliding boxcar function over the data.  You set the width of
%box.
%m is the 1D data
%c is the width of the box
%Jonathan A. N. Fisher 2008

function G = boxcar(m, c)

%G will be the smoothed pixel trace, so let's initialize it
G = zeros(1,max(size(m)));  

for i=1:max(size(m));

        window_min = i-round((c-1)/2);
        window_max = i+round((c-1)/2);
        
        if window_min<1
            window_min = 1;
        elseif window_max>max(size(m))
            window_max = max(size(m));
        else
        end

        G(i) = sum(m(window_min:window_max))/c;
end




