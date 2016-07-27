function [C, Ceq] = rnotzero(r)
C = double((r==0)-1); 
Ceq = double(r==0);