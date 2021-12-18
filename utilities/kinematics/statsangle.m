function [m,s] = statsangle(angle);
angle = angle(:);
n = length(angle);
a0 = angle(1);
for i=1:n
    angle(i) = addangle(angle(i),-a0);
end
m = a0+mean(angle);
s = std(angle);
return;


    
