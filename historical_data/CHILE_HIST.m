

c=wod_rd('chile/ocldb1671559124.5072.OSD');%,'params',[1 2 3]);
 
ctd=wod_rd('chile/ocldb1671559124.5072.CTD','params',[1 2 3]);
 
%%
ii=c.latitude>-43 & c.longitude>-72.6 & c.latitude<-41.9;

clf;orient landscape;wysiwyg;
subplot(1,2,1);
plot(c.longitude,c.mtime,'o');
line(c.longitude(ii),c.mtime(ii),'marker','o','linest','none','linewi',2,'color','m');
line(ctd.longitude,ctd.mtime,'marker','x','color','r','linest','none');
axdate('y',20)
axdeg('xW');xlabel('Longitude');
ylabel('Year');
legend('Hydrographic Station','...in Comau','CTD');
title('Chile Fjord area data in NODC archive');

subplot(1,2,2);
m_proj('lambert','lon',[-78 -72],'lat',[-48 -40]);
m_gshhs_i('patch',[.7 .7 .7],'edgecolor','none');
m_grid;
m_line(c.longitude,c.latitude,'marker','o','color','b','linest','none');
m_line(ctd.longitude,ctd.latitude,'marker','x','color','r','linest','none');
m_line(c.longitude(ii),c.latitude(ii),'marker','o','linest','none','linewi',2,'color','m');

%%
c
print -dpng  ChileFjord


%%

%ii=c.latitude>-43 & c.longitude>-72.6 & c.latitude<-41.9 & c.mtime<;

% Hudson data
ii=c.NODCInstitution==152 & c.latitude>-43 & c.longitude>-72.6 & c.latitude<-41.9;

clf;
m_proj('lambert','lon',[-78 -72],'lat',[-48 -40]);
m_gshhs_i('patch',[.7 .7 .7],'edgecolor','none');
m_grid;
m_line(c.longitude(ii),c.latitude(ii),'marker','o','color','b','linest','-');

%%

clf;
subplot(2,3,1);
plot(c.temp(:,ii),c.depth(:,ii));
set(gca,'ydir','reverse');
xlabel('Temperature/^oC');
ylabel('Depth/dbar');

subplot(2,3,2);
plot(c.sal(:,ii),c.depth(:,ii));
line(c.BucketSal(ii),0,'marker','o');
set(gca,'ydir','reverse');
xlabel('Salinity/ppt');

subplot(2,3,3);
plot(c.ox(:,ii),c.depth(:,ii));
set(gca,'ydir','reverse');
xlabel('O_2/(\mumol/kg)');

subplot(2,3,4);
plot(c.sal(:,ii),c.temp(:,ii));
xlabel('Salinity/ppt');
ylabel('Temperature/^oC');

subplot(2,3,[5 6]);
m_proj('lambert','lon',[-75 -72],'lat',[-43 -41]);
m_gshhs_f('patch',[.7 .7 .7],'edgecolor','none');
m_grid;
m_line(c.longitude(ii),c.latitude(ii),'marker','o','color','b','linest','-');
title('CCGS Hudson (1970)');

%%
print -dpng Hudson

%%
plot(c.mtime(ii),c.ox(:,ii),'ob'); %%,ctd.mtime(ic),ctd.ox(:,ic),'.r');
axdate;

 


