function [ Pmp,Vmp,Imp ] = SolveMpp( ipv_func, vpv_range, errTolerance )
%SOLVEMPP Summary of this function goes here
%   Detailed explanation goes here
    func = @(vpv) -vpv*SolveIpv(ipv_func,vpv);
    [min,Vmp] = golden_search(func,vpv_range,errTolerance);
    Imp = SolveIpv(ipv_func, Vmp);
    Pmp=-min;
end

