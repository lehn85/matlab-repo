function Ipv = SolveIpv( ipv_func, Vpv )
%IPVFUNC Solve Ipv equation
%   ipv_func is a function of (Vpv,Ipv)
%   Detailed explanation goes here
% Iph_c = @(G,Tc) pv.Isc_c*(1+pv.Ki*(Tc-pv.Tstc))*G;
% I0_c = @(Tc) pv.Isc_c*(1+pv.Ki*(Tc-pv.Tstc))/(exp(pv.q*(pv.Voc_c*(1+pv.Kv*(Tc-pv.Tstc)))/(pv.kb*Tc*pv.A))-1);
% Ipv_func = @(G,Tc,Vpv,Ipv,Rs,Rp) pv.Np*Iph_c(G,Tc)-pv.Np*I0_c(Tc)*(exp(pv.q*(Vpv+Ipv*Rs)/(pv.Ns*pv.kb*Tc*pv.A))-1)-(Vpv+Ipv*Rs)/Rp;

Ipv=ipv_func(Vpv,0);
for i3=0:1:10
    Ipv=ipv_func(Vpv,Ipv);
end

%fix out of bound values
for i3=1:1:size(Vpv')
    if Ipv(i3)<0
        Ipv(i3)=0;
    end
end

end

