%% var
clear;% clear workspace
init_var();

%% define functions
pv.Isc_c=Isc/Np;
pv.Voc_c=Voc/Ns;
T0=273;%0 C degree
Tstc=T0+25;
pv.Tstc=Tstc;%standard test condition temp
pv.kb=1.381e-23;%Boltzmann const
pv.A=1.3;% ideal coeff
pv.q=1.6e-19;%electron charge
pv.Eg=1.12;%bang-gap energy
pv.Ki=Ki;
pv.Kv=Kv;
pv.Np=Np;
pv.Ns=Ns;
pv.Vmpp=Vmpp;
pv.Impp=Impp;
pv.Pmaxe=Pmaxe;

%Functions: G (kW/m2), Tc *K
% G-kW/m2;Tc-K;Isc(A)
pv.Iph_c = @(pv,G,Tc) pv.Isc_c*(1+pv.Ki*(Tc-pv.Tstc))*G;
% I reverse saturation
%Irs_c = @(Tc) Isc_c/(exp(q*Voc_c/(kb*A*Tc))-1);
% I saturation
%I0_c = @(Tc) Irs_c(Tc)*((Tc/Tstc)^3)*exp(q*Eg*(1/Tstc-1/Tc)/(kb*A));
% VILLALVA et al.: COMPREHENSIVE APPROACH TO MODELING AND SIMULATION OF PHOTOVOLTAIC ARRAYS
%Vth = @(Tc) (k*Tc)/q;
pv.I0_c = @(pv,Tc) pv.Isc_c*(1+pv.Ki*(Tc-pv.Tstc))/(exp(pv.q*(pv.Voc_c*(1+pv.Kv*(Tc-pv.Tstc)))/(pv.kb*Tc*pv.A))-1);
% output pv current
% Ipv_func = @(G,Tc,Vpv,Ipv) Np*Iph_c(G,Tc)-Np*I0_c(Tc)*(exp(q*(Vpv./Ns+Ipv.*Rs/Np)/(kb*Tc*A))-1)-(Np*Vpv./Ns+Ipv.*Rs)/Rp;
% use Rs Rp equivlent
pv.Ipv_func_Rspx = @(pv,G,Tc,Vpv,Ipv,rsx,rpx) pv.Np*pv.Iph_c(pv,G,Tc)-pv.Np*pv.I0_c(pv,Tc)*(exp(pv.q*(Vpv+Ipv*rsx)/(pv.Ns*pv.kb*Tc*pv.A))-1)-(Vpv+Ipv*rsx)/rpx;
% use Rs Rp of cell
% Ipv_func_Rspx = @(G,Tc,Vpv,Ipv,rsx,rpx) Np*Iph_c(G,Tc)-Np*I0_c(Tc)*(exp(q*(Vpv./Ns+Ipv.*rsx/Np)/(kb*Tc*A))-1)-(Np*Vpv./Ns+Ipv.*rsx)/rpx;

%% Rs, Rshunt estimation
GG=1;
TT=Tstc;
pv.GG=GG;
pv.TT=TT;
% Rp0=Vmpp/(Isc-Impp)-(Voc-Vmpp)/Impp;
% Rs Rp equivlent
pv.Rp_func = @(pv,rs) (pv.Vmpp+pv.Impp*rs)/(pv.Np*pv.Iph_c(pv,GG,TT)-pv.Np*pv.I0_c(pv,TT)*(exp(pv.q*(pv.Vmpp+pv.Impp*rs)/(pv.Ns*pv.kb*TT*pv.A))-1)-pv.Pmaxe/pv.Vmpp);
% Rs Rp of cell
% Rp_func = @(rs) (Vmpp*Np/Ns+Impp*rs)/(Np*Iph_c(GG,TT)-Np*I0_c(TT)*(exp(q*(Vmpp/Ns+Impp*rs/Np)/(kb*TT*A))-1)-Pmaxe/Vmpp);
% Rpx=Rp0;
Rsx=0;
% Ipv_func_Rspx = @(G,Tc,Vpv,Ipv,rsx,rpx) Np*Iph_c(G,Tc)-Np*I0_c(Tc)*(exp(q*(Vpv./Ns+Ipv.*rsx/Np)/(kb*Tc*A))-1)-(Np*Vpv./Ns+Ipv.*rsx)/rpx;
errPmax=1000;
errTolerance=1e-6;
Rsmin=0;
Rsmax=1;
iter=0;
maxIter=50;
tau=double((sqrt(5)-1)/2);
f_x=[0;0];
Rsbound=[Rsmin+(1-tau)*(Rsmax-Rsmin);Rsmin+tau*(Rsmax-Rsmin)];
while(iter<maxIter && abs(Rsmin-Rsmax)>errTolerance)
    iter=iter+1;    
    for i2=1:1:2
        Rsx=Rsbound(i2);
        % calc Rp from Rs
        Rpx=pv.Rp_func(pv, Rsx);
        % solve Ipv=Ipv_func()
%         I_pv_m=SolveIpv(@(Vpv,Ipv) pv.Ipv_func_Rspx(pv,GG,TT,V_pv_m,I_pv_m,Rsx,Rpx),V_pv_m);        
        % find max point
%         mppP = 0;
%         mppV = 0;
%         mppI = 0;
%         for i3=1:1:size(V_pv_m')
%             if (mppP<V_pv_m(i3)*I_pv_m(i3))
%                 mppP = V_pv_m(i3)*I_pv_m(i3);
%                 mppV = V_pv_m(i3);
%                 mppI = I_pv_m(i3);
%             end
%         end
        [mppP,mppV,mppI]=SolveMpp(@(v,i) pv.Ipv_func_Rspx(pv,GG,TT,v,i,Rsx,Rpx),[0 Voc],1e-5);
        f_x(i2)=mppP;
    end
    %compare
    if (f_x(1)<f_x(2))
        Rsmax=Rsbound(2);
        Rsbound(2)=Rsbound(1);
        Rsbound(1)=Rsmin+(1-tau)*(Rsmax-Rsmin);        
    else
        Rsmin=Rsbound(1);
        Rsbound(1)=Rsbound(2);
        Rsbound(2)=Rsmin+tau*(Rsmax-Rsmin);
    end
    
    fprintf('iter=%d Pm1=%0.3f Pm2=%0.3f Imp=%0.4f Vmp=%0.4f Rsmin=%0.5f Rsmax=%0.5f Rp=%0.5f\n',iter, f_x(1),f_x(2),mppI, mppV, Rsmin,Rsmax, Rpx);
    %     PmaxList(i2)=mppP;
%     errPmax = abs(mppP-Pmaxe);
%     if (errPmax<errTolerance || Rsx>1)
%         break;
%     else
%         Rsx=Rsx+1e-5;
%     end
%     fprintf('Pm=%0.3f Imp=%0.4f Vmp=%0.4f Rs=%0.5f Rp=%0.5f\n',mppP,mppI, mppV, Rsx,Rpx);
end
Rs=(Rsmin+Rsmax)/2;
Rp=pv.Rp_func(pv,Rs);
fprintf('Calculated Pmax=%0.3f Impp=%0.4f Vmpp=%0.4f Rs=%0.5f Rp=%0.5f I0=%s\n ', mppP, mppI, mppV, Rs, Rp, Np*pv.I0_c(pv,Tstc));
% output pv current
Ipv_func = @(G,Tc,Vpv,Ipv) pv.Ipv_func_Rspx(pv,G,Tc,Vpv,Ipv,Rs,Rp);
%% build Pmpp(rs) graph
Rsmin=0;
Rsmax=0.5;
RsList=(Rsmin:1e-3:Rsmax)';
PmaxList(size(RsList),1)=0;
Pmaxmin=10000;
RsForPmaxmin=0;
for i=1:1:size(RsList)
    Rsx=RsList(i);
    Rpx=pv.Rp_func(pv,Rsx);
    [mppP,mppV,mppI]=SolveMpp(@(v,i) pv.Ipv_func_Rspx(pv,GG,TT,v,i,Rsx,Rpx),[0 Voc],1e-5);
    PmaxList(i)=mppP;
    if (Pmaxmin>mppP)
        Pmaxmin=mppP;
        RsForPmaxmin=Rsx;
    end
end
fprintf('Calculated Pmaxmin=%0.4f Rsx=%0.5f\r\n',Pmaxmin,RsForPmaxmin);
plot(RsList,PmaxList);
%% build lookup table
index1=0;
index2=0;
V_pv_m=0:0.1:Voc;
I_pv_m(size(V_pv_m))=0;
G_m=0:10:1000;%W/m2
T_m=0:1:75;%*C

tbl_pv(size(V_pv_m),size(G_m),size(T_m))=0;
tbl_pv_mppt_P(size(G_m),size(T_m))=0;
tbl_pv_mppt_V(size(G_m),size(T_m))=0;
tbl_pv_mppt_I(size(G_m),size(T_m))=0;
% make P-V table
for Gi=G_m
    index1 = index1+1;
    index2 = 0;
    for Ti=T_m
        index2 = index2+1;
        I_pv_m=Ipv_func(Gi/1000,Ti+T0,V_pv_m,0);
        for i3=0:1:10
            I_pv_m=Ipv_func(Gi/1000,Ti+T0,V_pv_m,I_pv_m);
        end
        %fix outside
        for i3=1:1:size(V_pv_m')
            if V_pv_m(i3)>(Voc*(1+Kv*(Ti+T0-Tstc))) || I_pv_m(i3)<0
                break;
            end
        end
        for i4=i3:1:size(V_pv_m')
            I_pv_m(i4)=0;
        end
        tbl_pv(:,index1,index2)=I_pv_m;
        mppP = 0;
        mppV = 0;
        mppI = 0;
        for i3=1:1:size(V_pv_m')
            if (mppP<V_pv_m(i3)*I_pv_m(i3))
                mppP = V_pv_m(i3)*I_pv_m(i3);
                mppV = V_pv_m(i3);
                mppI = I_pv_m(i3);
            end
        end
        tbl_pv_mppt_P(index1,index2)=mppP;
        tbl_pv_mppt_V(index1,index2)=mppV;
        tbl_pv_mppt_I(index1,index2)=mppI;
    end
end
%% plot (by G)
TstcIndex=26;
Glist=[101;81;61;41;21];
PVmppt_G(size(Glist),3)=0;
IlistG(size(V_pv_m'),size(Glist))=0;
PlistG(size(V_pv_m'),size(Glist))=0;
for i=1:1:size(Glist)
    PVmppt_G(i,1)=tbl_pv_mppt_V(Glist(i),TstcIndex);
    PVmppt_G(i,2)=tbl_pv_mppt_P(Glist(i),TstcIndex);
    PVmppt_G(i,3)=tbl_pv_mppt_I(Glist(i),TstcIndex);
    IlistG(:,i)=tbl_pv(:,Glist(i),TstcIndex);
    PlistG(:,i)=tbl_pv(:,Glist(i),TstcIndex).*V_pv_m';
end
createfigure_I_V_byG(V_pv_m,IlistG,PVmppt_G(:,1),PVmppt_G(:,3));
set(gcf,'Position',[50   100   800   600]);
createfigure_P_V_byG(V_pv_m,PlistG,PVmppt_G(:,1),PVmppt_G(:,2));
set(gcf,'Position',[800   100   800   600]);

%% plot (by T)
GstcIndex=101;
Tlist=[26;36;46;56;66];
PVmppt_T(size(Tlist),3)=0;
IlistT(size(V_pv_m'),size(Tlist))=0;
PlistT(size(V_pv_m'),size(Tlist))=0;
for i=1:1:size(Tlist)
    PVmppt_T(i,1)=tbl_pv_mppt_V(GstcIndex,Tlist(i));
    PVmppt_T(i,2)=tbl_pv_mppt_P(GstcIndex,Tlist(i));
    PVmppt_T(i,3)=tbl_pv_mppt_I(GstcIndex,Tlist(i));
    IlistT(:,i)=tbl_pv(:,GstcIndex,Tlist(i));
    PlistT(:,i)=tbl_pv(:,GstcIndex,Tlist(i)).*V_pv_m';
end
createfigure_I_V_byT(V_pv_m,IlistT,PVmppt_T(:,1),PVmppt_T(:,3));
set(gcf,'Position',[50   100   750   640]);
createfigure_P_V_byT(V_pv_m,PlistT,PVmppt_T(:,1),PVmppt_T(:,2));
set(gcf,'Position',[800   100   750   640]);
