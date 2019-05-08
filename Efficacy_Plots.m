%% Efficacy plots
% Main comparison plot
dfigure('Position',[754 185 321 566]);
a1 = subplot(4,1,1);
p1=semilogy(oK.fssDc/1e6,oK.Psn_fDc);  hold all;
p2=semilogy(oK.fss/1e6,oK.Psn);  hold all;
uistack(p2,'bottom');
% perfectly-corrected signal without noise
ot = Correct_Signal(Vn(1:length(oK.p0))-wn(1:length(oK.p0)),dt,p0(1:length(oK.p0)),pD(1:length(oK.p0))); 
Pa = abs(As).^2;                                                                  
Pp = ot.Psn_fDc(ot.pklocs)/(ot.NDc^2*ot.dtDc/ot.NDc) - 0*2*1^2/ot.NDc;            
[xc,lgs] = xcorr(Pa,Pp); [~,ml]=max(xc); sa = lgs(ml); % amt to shift Pe by
if sa<0, Pp = Pp(1-sa:end);
else,    Pp = [zeros(sa,1);Pp]; end
if length(Pp)<length(Pa), Pp = [Pp;zeros(length(Pa)-length(Pp),1)];
else,                     Pp = Pp(1:length(Pa)); end
semilogy((mean(f0)+ns*mean(fD))/1e6,Pp*oK.NDc^2*oK.dtDc/oK.NDc,'Color',dColor(7));        % true amplitudes
xlim([-290000000 290000000]/1e6); ylim([1.184e-09 0.004362]);
xyt('Frequency (MHz)','Power (a.u.)','');

% Zoomed comparison plot
a2=subplot(4,1,2);
p1=semilogy(oK.fssDc/1e6,oK.Psn_fDc);  hold all;
p2=semilogy(oK.fss/1e6,oK.Psn);  hold all;
uistack(p2,'bottom');
semilogy((mean(f0)+ns*mean(fD))/1e6,Pp*oK.NDc^2*oK.dtDc/oK.NDc,'Color',dColor(7));        % true amplitudes
xyt('Frequency (MHz)','Power (a.u.)','');
xlim([-19.83 21.47]); ylim([3.806e-09 0.004024]);

% Frequency plots. Note that f0 can be off by an integer number of fr's
ts=[0:N-1]'*dt;
a3=subplot(4,1,3);
Nd=round(mean((oK.f0-f0)./oK.fD)); % which line we're locked to is essentially arbitrary, determine number of rep rates off for plot
plot(ts(1:length(oK.f0))/1e-6,(oK.f0-Nd*oK.fD)/1e6); hold all;
p2 = plot(ts/1e-6,f0/1e6,'Color',dColor(7)); uistack(p2,'bottom');
xyt('','f_0 (MHz)','');
a3.XTickLabel=[];
a4=subplot(4,1,4);
plot(ts(1:length(oK.f0))/1e-6,(oK.fD)/1e6); hold all;
p2 = plot(ts/1e-6,fD/1e6,'Color',dColor(7)); uistack(p2,'bottom');

xyt('Time (\mus)','f_r (MHz)','');
set(a2,'Position',[0.2461 0.8055 0.553 0.1574]);
set(a1,'Position',[0.13 0.4274 0.775 0.3041]);
set(a3,'Position',[0.13 0.2297 0.775 0.1158]);
set(a4,'Position',[0.13 0.11 0.775 0.1214]);
%% Residual plots
Pe = (oK.Psn_fDc(oK.pklocs)-0*oK.npsd)/(oK.NDc^2*oK.dtDc/oK.NDc) - 2*1^2/ot.NDc;
[xc,lgs] = xcorr(Pa,Pe); [~,ml]=max(xc); sa = lgs(ml); % amt to shift Pe by
if sa<0, Pe = Pe(1-sa:end);
else,    Pe = [zeros(sa,1);Pe]; end
if length(Pe)<length(Pa), Pe = [Pe;zeros(length(Pa)-length(Pe),1)];
else,                     Pe = Pe(1:length(Pa)); end

dfigure; loglog(Pp,abs(Pe-Pp)./Pp,'Marker','.','LineStyle','none');
[~,so]=sort(Pa); hold all; plot(Pp(so),sqrt(4*1^2/oK.NDc./Pp(so).*(1+1./Pp(so))));
[~,so]=sort(Pa); hold all; plot(Pp(so),sqrt(4*1^2/oK.NDc./Pp(so).*(1+1./Pp(so)/oK.NDc)),'Color',dColor(7));

xlim([3.41e-06 67.98]); ylim([0.0001374 100]);
legend('Power residuals','Incoherent limit','Coherent limit'); legend boxoff;
xyt('Power (a.u.)','Residual','');