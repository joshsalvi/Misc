
%wrapphase.m is a function that "unwraps" phase traces that are wrapped at
%+/- 2pi
%
%JANF 2010

function corrected = unwrapphase(tracein);

corrected = tracein;

phaseadd = 0;

for i = 2:length(tracein)
   
    corrected(i) = tracein(i)+phaseadd;
    
    if tracein(i-1)>tracein(i)+3
        corrected(i) = tracein(i)+2*pi;
        phaseadd = phaseadd+2*pi;
    elseif tracein(i-1)<tracein(i)-3;
        corrected(i) = tracein(i)-2*pi;
        phaseadd = phaseadd-2*pi;
    else
    end
end