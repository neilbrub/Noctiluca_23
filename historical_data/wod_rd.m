function [c,dat,profindex]=wod_rd(filename,arg,iparams);
% WOD_RD reads NODC 'WOD' format file (WOD2001, WOD2005, and WOD13)
%        WOD_RD(FILENAME) returns all data.
%        WOD_RD(FILENAME,'info') gives stats
%        WOD_RD(FILENAME,'params',DATA) returns parameters
%        specified by DATA only, where DATA is a vector  containing
%        one or more of:
%                                   	  WOD13
%           1-Temperature            C 
%           2-Salinity               -
%           3-Oxygen                 ml/l
%           4-Phosphate              uM
%           5-Total Phosphorus      	   No
%           6-Silicate               uM
%           7-Nitrite               	   No
%           8-Nitrate (+Nitrite)     uM
%           9-pH                      - 
%           10-Ammonia              	   No
%           11-Chlorophyll           ug/L
%           12-Phaeophytin          	   No
%           13-Primary Prod         	   No
%           14-Biochem              	   No
%           15-LightC14             	   No 
%           16-DarkC14              	   No
%           17-Alkalinity           meq/l
%           18-POC                  	   No
%           19-DOC                  	   No
%           20-pCO2                 uatm
%           21-tCO2 (DIC)           mM
%           22-XCO2sea              	   No
%           23-NO2NO3               	   No
%           24-Transmissivity        /m
%           25-Pressure             dbar
%           26-Air temp               C         <- Conflict between wodASC (which sez 'Conductivity') and Documentation
%           27    Co2warming          C
%           28    xCO2-atmos         ppm
%           29    Air pressure       mbar
%           30    Lat
%           31    Long
%           32    Julian Year-day
%           33-Tritium               TU
%           34-Helium'               nM
%           35-DeltaHe                %
%           36-DeltaC14             permille
%           37-DeltaC13             permille
%           38-Argon                 nM
%           39-Neon                  nM
%           40-CFC11                 pM
%           41-CFC12                 pM
%           42-CFC113                pM   
%           43-O18                  permille
% 
%       e.g. DATA=[1 2 9] returns all profiles with temp, salinity, or pH.
% 
% Currently the code will decode the biological data that may also be
% stored but none of this data is returned to the calling function.
%

% Rich pawlowicz (rich@eos.ubc.ca) 18/Dec/04
%   OCt/2013 - WOD13 format handling added
%   June/2014  - Better handling of some format wrinkles
%                added code to read headers first to get stats
%                before reading data (good for CTD data)

% Standard depths
stdz=[ 0., 10., 20., 30., 50., 75., 100., 125., 150., ...
       200., 250., 300., 400., 500., 600., 700., 800., 900., ...
       1000., 1100., 1200., 1300., 1400., 1500., 1750., 2000., ...
       2500., 3000., 3500., 4000., 4500., 5000., 5500., 6000., ...
       6500., 7000., 7500., 8000., 8500., 9000. ];

stdnam={'temp','sal','ox','po4','Tpo4','si','no2','no3','pH','nh4',...
        'chl','phaeo','pprod','biochem','lightc14','darkc14','alk',...
	'POC','DOC','pco2','dic','xco2sea','no2no3','xmiss','press','air_temp',...
	'co2warm','xco2_atm','air_press','lat','lon','julday','tritium','helium','deltaHe','deltaC14','deltaC13',...
	'Ar','Ne','cfc11','cfc12','cfc113','o18'};
	
stdunit={'C','ppt/PSU','ml/l','uM','uM','uM','uM','uM','-','uM',...
         'ug/l','','','','','','meq/l','','','uatm','mM','','','/m','dbar','C',...
	 'C','ppm','mbar','deg','deg','day','TU','nM','%','permille','permille',...
	 'nM','nM','pM','pM','pM','permille'};	
info=0;
params=-1; %[1:26];  % -1 means all

if nargin>1,
 switch arg,
   case 'info',
     info=1;
   case 'params',
     params=iparams;
 end;
end;
     
     
if info,        

  % Suck in all the data
  fd=fopen(filename,'r');
  dat=fread(fd,'char');
  fclose(fd);

  % Get rid of line separators
  dat(dat==10 | dat==13 )=[];

  fmt=dat(1);
  if fmt=='A',
    disp('WOD01 format');
  elseif fmt=='B',
    disp('WOD05 format');
  elseif fmt=='C',
    disp('WOD13 format');    
  else
    error('WOD 1998 format?');
  end;

  % Find beginning of all casts.
  % This works a lot until a lot of 'C's started showing up in a WOD13 format file
  %profindex=find(dat(1:end-1)==fmt & (dat(2:end)>='0' & dat(2:end)<='9') & (['0';dat(1:end-2)]<'A' | ['0';dat(1:end-2)]>'Z') );

  % So try this:
  % format is always first character in the line, but sometimes we can have a 'C' there for other reasons
  profindex=find(dat(1:80:end)==fmt & (dat(2:80:end)>='0' & dat(2:80:end)<='9') )*80-79;


  nprof=length(profindex);
  profindex=[profindex;length(dat)+1];

  fprintf('%d profiles in file\n',nprof);

  c=struct('info',filename,...
         'format',char(fmt),...
	 'nprofiles',nprof,....
         'stationID',repmat(NaN,1,nprof),...
         'country',repmat(NaN,1,nprof),...
	 'countrycode',repmat(' ',nprof,2),...
	 'cruise',repmat(NaN,1,nprof),...
	 'mtime',repmat(NaN,1,nprof),...
	 'gtime',repmat(NaN,6,nprof),...
	 'latitude',repmat(NaN,1,nprof),...
	 'longitude',repmat(NaN,1,nprof),...
	 'nlevels',repmat(NaN,1,nprof),...
	 'paramnames',cell(1,1),...
	 'params',repmat(0,length(stdnam),nprof)   );
  c.paramnames=stdnam';	 
  c.paramunits=stdunit';	 
  
  
else

  % First read file once to get info
  fprintf(' First, scan for info on available data\n');
  [finfo,dat,profindex]=wod_rd(filename,'info');
  
  nprof=length(profindex)-1;
  fmt=finfo.format;
  maxlevels=max(finfo.nlevels);
  if params<0,
    params=find(any(finfo.params,2));
  end;  
  fprintf('\n#params = %d, maxlevels=%d \n param#---name----#profiles--units\n',length(params),maxlevels);

  c=struct('data',filename,...
         'format',char(fmt),...
	 'nprofiles',nprof+1,...
	 'stationID',repmat(NaN,1,nprof),...
         'country',repmat(NaN,1,nprof),...
	 'countrycode',repmat(' ',nprof,2),...
	 'WODplatform',repmat(NaN,1,nprof),...
	 'NODCInstitution',repmat(NaN,1,nprof),...
	 'cruise',repmat(NaN,1,nprof),...
	 'mtime',repmat(NaN,1,nprof),...
	 'gtime',repmat(NaN,6,nprof),...
	 'latitude',repmat(NaN,1,nprof),...
	 'longitude',repmat(NaN,1,nprof),...
	 'bottomdepth',repmat(NaN,1,nprof),...
	 'SECCHI',repmat(NaN,1,nprof),...
	 'BucketSal',repmat(NaN,1,nprof),...
	 'biodata_available',repmat(0,1,nprof),...
	 'biodata_ntaxa',repmat(NaN,1,nprof) ,...
	 'nlevels',repmat(NaN,1,nprof),...
	 'depth',repmat(NaN,maxlevels,nprof),...
	 'params',finfo.params(params,:),...
	 'paramnames',cell(1,1) );
  c.paramnames=finfo.paramnames(params)';
  for k=1:length(params),
    eval(['c.' stdnam{params(k)} '=repmat(NaN,maxlevels,nprof);']);
    fprintf(' %2d %10s  -> %5d     (%s)\n',params(k), stdnam{params(k)}, sum(finfo.params(params(k),:))  ,  stdunit{params(k)});
  end;  
  fprintf('Reading data\n');

end;
	 
for k=1:nprof,
  if rem(k,round(nprof/10))==0, fprintf('.'); end;
  
  ichar=char(dat(profindex(k)+1:profindex(k+1)-1)');
  [nbytes,ptr] =read_int(ichar,1);
%%fprintf('length - %d, %d\n',length(ichar),nbytes);
  [station,ptr]=read_int(ichar,ptr);
  if fmt=='A' | fmt=='B',
    [country,ptr]=read_int(ichar,ptr,2);
  elseif fmt=='C',
    countrycode=ichar(ptr+[0:1]);
    ptr=ptr+2;
    country=countrymatch(countrycode);
    if isempty(country),country=NaN; end;
  end;  
  [cruise,ptr] =read_int(ichar,ptr);
  [year,ptr]   =read_int(ichar,ptr,4);
  [mon,ptr]    =read_int(ichar,ptr,2);
  [dy,ptr]     =read_int(ichar,ptr,2);
  [tim,ptr]    =read_float(ichar,ptr);
  [lat,ptr]    =read_float(ichar,ptr);
  [long,ptr]   =read_float(ichar,ptr);
 
  c.stationID(k)=station;
  c.country(k)=country;
  c.countrycode(k,:)=countrymatch(country);
  c.cruise(k)=cruise;
  c.gtime(:,k)=[year mon dy fix(tim) (tim-fix(tim))*60 0]';
  if isfinite(tim), c.mtime(k)=datenum(c.gtime(:,k)');
  else            c.mtime(k)=datenum(year,mon,dy); end;
  c.latitude(k)=lat;
  c.longitude(k)=long;

  [nlevels,ptr]=read_int(ichar,ptr);
  [isoor,ptr]  =read_int(ichar,ptr,1);
  [nvar,ptr]   =read_int(ichar,ptr,2);

  c.nlevels(k)=nlevels;
%%fprintf('Meta\n');

  % Read meta-data for parameters in file (currently only the name of the variables
  % stored in this profile are saved in 'c')
  
  for l=1:nvar,
    [varcode,ptr]=read_int(ichar,ptr);  % Variable code (Table 3)
    if varcode>0,
     if info, c.params(varcode,k)=1; end;
     colnum(l)=varcode;
    end; 
    [qcode,ptr]  =read_int(ichar,ptr,1);   % Quality control flag (Table 12):
                                           % Various tests against climatology
					   % - not saved here!
    [nmeta,ptr]  =read_int(ichar,ptr);     % Number of variable-specific metadata
    for m=1:nmeta,
      [varcode,ptr]=read_int(ichar,ptr);   %Variable-specific code (Table 5)
                                           % - not saved here!
               % 1 Accession #
	       % 2 Project
	       % 3 Scale  (ITS90 or IPTS68, etc.)
	       % 4 Institution
	       % 5 Instrument (brand name) 
               % 6 methods
               % 8 Originators Units
	       % 10 Equilibrator type (for CO2 measurements)
	       % 11 Filter type and size (for filtered samples)
	       % 12 Incubation time (for primary productuvuty measurements)
	       % 13 CO2 Sea warming (when sample analysed in lab)
	       % 15 Analysis temperature (for CO2 measurement)
	       % 16 Uncalibrated  (set to 1 if instrument uncalibrated)
	       % 17 Contains nitrite (set to 1 if nitrate+nitrite)
	       % 18 Standard Seawater Batch number
	       % 19 Adjustment (in ARGO profile data)
      [varval,ptr] =read_float(ichar,ptr); % ...and its value
    end;  
    %end;
  end;
  
  %% This is all we need to read to see included variables.
  %% Keep going only to get actual data!
  
  if ~info,
  
  % nbyte for each header appears to be the number of bytes
  % not including the amount used to store nbyte - i.e. just add
  % that much to skip reas of block.
  
  % Read character data for PI block
  %  - currently none of this is saved to 'c'
  
  [nbyte,ptr]=read_int(ichar,ptr);
  if nbyte>0,
    [nentry,ptr]=read_int(ichar,ptr,1);
    for l=1:nentry,
      [ctype,ptr]=read_int(ichar,ptr,1);
      if ctype<3,  % ctype==1 or ==2
        % originators cruise code if ==1 or station code if ==2
        [ndat,ptr]=read_int(ichar,ptr,2);
	[chardat,ptr]=read_int(ichar,ptr,ndat);
      else  % ctype==3 means this is PI info
        [ndat,ptr]=read_int(ichar,ptr,2);  % Number of PI names
	for  m=1:ndat,
	 [varcode,ptr]=read_int(ichar,ptr);  % Variable code
	 [PIcode,ptr]=read_int(ichar,ptr);   % PI code for that variable
	end;
      end;
    end;
  end;
  
  
  % Read secondary header:
  % This contains much info that would be in the header of the
  % individual data files - local weather conditions, platform
  % codes, digitization methods (for really old data), bottom
  % depth
  %%fprintf(' sec start %d ',ptr);  
  [nbyte ,ptr]=read_int(ichar,ptr);
  if nbyte>0,
    [nentry,ptr]=read_int(ichar,ptr);
    for l=1:nentry,
      [varcode,ptr]=read_int(ichar,ptr);
      [headval,ptr]=read_float(ichar,ptr);
      switch varcode,
        case 3,
	  c.WODplatform(k)=headval;
	case 4,
	  c.NODCInstitution(k)=headval;
        case 10,
	  c.bottomdepth(k)=headval;
	case 15,
	  c.SECCHI(k)=headval;  
	case 47,
	  c.BucketSal(k)=headval;  
       end;	  
    end;
  end;
  %%fprintf(' end %d .... %d\n',ptr,nbyte);
  
  % Read biological header
  % This is completely ignored!
  
  optr=ptr;  
  [nbyte ,ptr]=read_int(ichar,ptr);
  
  if nbyte>0,
     % fprintf(' bio start %d ',optr);  
     [nentry,ptr]=read_int(ichar,ptr);
     c.biodata_available(k)=nentry;
     for l=1:nentry,
      [varcode,ptr]=read_int(ichar,ptr);   % Bio header variables (Table 6)
                                           % not saved here!
	% 1 Volume filtered
	% 2 sampling duration (mins)
	% 3 mesh size (um)
	% 4 Type of tow (vert/horiz/oblique)
	% 5 Large removed volume (ml)
	% 6 Large plankters removed flag
	% 7 Gear code
	% 8 Sampler volume (leters)
	% 9 Net mouth area (m^2)
	% 10 Preservative
	% 11 Weight method
	% 12 large removed length criteria (cm)
	% 13 Count method
	% 14 Tow distance (m)
	% 15 Average tow speed (knots)
	% 16 Sampling start time (GMT)
	% 18 Flowmeter type
	% 19 Flowmeter calibration
	% 20 Counting INStitution
	% 21 Vouncher Institution
	% 22 Wire angle start (degrees)
	% 23 Wire angle end (degrees)
	% 24 Depth determination method
	% 25 Volume method
	% 30 Accession number				   
      [biohead,ptr]=read_float(ichar,ptr);
     end; 
          	 	      
     % Read taxonomic data
     % This is completely ignored!
  
    [ntaxa ,ptr]=read_int(ichar,ptr);
    %fprintf(' taxo start %d %s',ptr,ichar(ptr+[0:10]));
    c.biodata_ntaxa(k)=ntaxa;
    for l=1:ntaxa,
       [nentry,ptr]=read_int(ichar,ptr); % Number of entrieds per taxa
       for m=1:nentry,
        [varcode,ptr]=read_int(ichar,ptr);  %  Taxa code (Table 7)
        [biohead,ptr]=read_float(ichar,ptr);
        [bioflags,ptr]=read_int(ichar,ptr,2);
       end; 
   end;        	 	
   %fprintf(' end %d....%d\n',ptr,ntaxa);
  end; 
    

  for l=1:nlevels,
   if isoor==0 || fmt=='C', % not at standard depths
     %fprintf('dep? %s\n',ichar(ptr+[0:10]))
     [depval,ptr]=read_float(ichar,ptr);
     c.depth(l,k)=depval;
   else
     c.depth(l,k)=stdz(l);      
   end;  
   %fprintf('nlevel=%d, dep=%f',l,c.depth(l,k));    
   [deperror,ptr]=read_int(ichar,ptr,1);
   [origdeperror,ptr]=read_int(ichar,ptr,1);
   for m=1:nvar,
     [val,ptr]=read_float(ichar,ptr);
     %fprintf('(%d,%f)',m,val);
     if isfinite(val),
       [QCflag,ptr]=read_int(ichar,ptr,1);
       [Oflag,ptr]=read_int(ichar,ptr,1);
     end;
     if any(find(colnum(m)==params)),
       eval(['c.' stdnam{colnum(m)} '(l,k)=val;']);
     end;	 
   end;  % variables in each level 
   %fprintf('\n');   
  end; % of levels
    
  end;% reading full profile data instead of just info alone
   
end; % of this profile

  
  
    
return;

function [val,ptr]=read_int(ichar,ptr,len);
% READS an integer

if nargin<3,
  field_bytes=sscanf(ichar(ptr),'%d');
%fprintf('v%5s %10s\n',ichar(ptr+[0:field_bytes]),ichar(ptr:min(end,ptr+10)));    
  if field_bytes>0,
    val=sscanf(ichar(ptr+[1:field_bytes]),'%d');
  else
    val=0;
  end;  
  ptr=ptr+field_bytes+1;
else
%fprintf('f%5s %10s\n',ichar(ptr+[0:len-1]),ichar(ptr:min(end,ptr+10)));    
  val=sscanf(ichar(ptr+[0:len-1]),'%d');
  ptr=ptr+len;
end;
  
return;

function [val,ptr]=read_float(ichar,ptr);
% READS a float
if ichar(ptr)=='-',
%fprintf('%s\n',ichar(ptr+[0:10]));
   val=NaN;
   ptr=ptr+1;
else
   sig=sscanf(ichar(ptr),'%d');
   tot=sscanf(ichar(ptr+1),'%d');  
   prec=sscanf(ichar(ptr+2),'%d');
   val=sscanf(ichar(ptr+2+[1:tot]),'%d')/(10^prec);
   ptr=ptr+tot+3;
end;
      
function c=countrymatch(country);

countrycodes=[...	 
'DE';% GERMANY
'DU';%  EAST GERMANY
'AR';%  ARGENTINA
'AU';%  AUSTRALIA
'AT';%  AUSTRIA
'BE';%  BELGIUM
'BR';%  BRAZIL
'BG';%  BULGARIA
'CA';%  CANADA
'CL';%  CHILE
'TW';%  TAIWAN
'CO';%  COLOMBIA
'KR';%  KOREA; REPUBLIC OF
'DK';%  DENMARK
'EG';%  EGYPT
'EC';%  ECUADOR
'ES';%  SPAIN
'US';%  UNITED STATES
'FI';%  FINLAND
'FR';%  FRANCE
'GR';%  GREECE
'IN';%  INDIA
'ID';%  INDONESIA
'IE';%  IRELAND
'IS';%  ICELAND
'IL';%  ISRAEL
'IT';%  ITALY
'JP';%  JAPAN
'LB';%  LEBANON
'LR';%  LIBERIA
'MG';%  MADAGASCAR
'MA';%  MOROCCO
'MX';%  MEXICO
'NO';%  NORWAY
'NC';%  NEW CALEDONIA
'NZ';%  NEW ZEALAND
'PK';%  PAKISTAN
'NL';%  NETHERLANDS
'PE';%  PERU
'PH';%  PHILIPPINES
'PL';%  POLAND
'PT';%  PORTUGAL
'RO';%  ROMANIA
'GB';%  GREAT BRITAIN
'CN';%  CHINA
'SE';%  SWEDEN
'TH';%  THAILAND
'TN';%  TUNISIA
'TR';%  TURKEY
'SU';%  SOVIET UNION
'ZA';%  SOUTH AFRICA
'UY';%  URUGUAY
'VE';%  VENEZUELA
'YU';%  YUGOSLAVIA
'99';%  UNKNOWN
'AG';%  ANTIGUA
'DZ';%  ALGERIA
'AO';%  ANGOLA
'BB';%  BARBADOS
'BS';%  BAHAMAS
'CR';%  COSTA RICA
'CU';%  CUBA
'CY';%  CYPRUS
'EE';%  ESTONIA
'FJ';%  FIJI
'GH';%  GHANA
'HN';%  HONDURAS
'HK';%  HONG KONG
'CI';%  COTE D'IVOIRE
'KW';%  KUWAIT
'LV';%  LATVIA
'LT';%  LITHUANIA
'MU';%  MAURITIUS
'MT';%  MALTA
'MC';%  MONACO
'MY';%  MALAYSIA
'MR';%  MAURITANIA
'NG';%  NIGERIA
'PA';%  PANAMA
'CD';%  CONGO; THE DEMOCRATIC REPUBLIC OF THE
'RU';%  RUSSIAN FEDERATION
'SA';%  SAUDI ARABIA
'SC';%  SEYCHELLES
'SN';%  SENEGAL
'SG';%  SINGAPORE
'SL';%  SIERRA LEONE
'VC';%  SAINT VINCENT AND THEN GRENADINES
'TO';%  TONGA
'TT';%  TRINIDAD AND TOBAGO
'UA';%  UKRAINE
'WS';%  SAMOA; WESTERN
'YE';%  YEMEN
'ZZ';%  MISCELLANEOUS ORGANIZATION
'MH';%  MARSHALL ISLANDS
'HR';%  CROATIA
'EU'];%  EUROPEAN UNION

if isstr(country), % if we pass in a string, get the number
  c=strmatch(country,countrycodes);

else % if we pass in a number, get the string
  if isfinite(country) & country<=size(countrycodes,1),
    c=countrycodes(country,:);
  else
    c='  ';
  end;  
end;
  

