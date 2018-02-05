function [result] = increasingArray(StartValue, endValue, number)
%INCREASINGARRAY Summary of this function goes here
%   Detailed explanation goes here
stepsize = (endValue - StartValue) / (number - 1);
result = StartValue : stepsize : endValue;
end

