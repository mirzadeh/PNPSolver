function [v] = inter(x1,x2,v1,v2,x)
theta = (x-x1)/(x2-x1);
v = (1-theta)*v1+theta*v2;
end