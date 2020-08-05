function [EEM_correct,EEM_Cutted_rrr2]=eemscat(X,MissRayleh,MissRaman,MissRayleh2)

%EEMSCAT for removing Rayleigh and Raman scattering from fluorescence
%        EEM data and interpolating the excised areas.
%
% INPUT:
%    X      X-array of EEMs. X is size IxJxK, where I is number of
%           samples, J emissions and K excitations. X has to be a
%           dataset object where the axisscales contain wavelengths.
%           You can convert an array to a dataset doing
%           X = dataset(X);
%           X.axisscale{2} = EmAx; % The emission wavelengths (nm)
%           X.axisscale{3} = ExAx; % The excitation wavelengths (nm)
%
% OPTIONAL INPUT
% RayeleighWitdth
%           RayeleighWitdth is a two-element vector defining how many nanometers 
%           to the left and right of the Rayleigh center is removed. Use
%           [0 0] to avoid removing. Default if not given is [25 25];
%
%
% RamanWidth 
%           As above but for Raman. Default is [0 0]. Assumes water samples
%           (for the position of the Raman scatter)
%
% Rayleigh2Width
%           As above but for second order Rayleigh scatter. Default is none 
%           [0 0].
%
%                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%					Example:
% NewEEM = eemscat( X,  [20 20],   [10 10],   [10 10]);
%                  Data Rayleigh   RamanWidth 2.orderRayl
%
%                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:
%
% NewEEM     EEM data with interpolated areas.
%
% Additional optional output
% EEMNaN    EEM data with missing data rather than interpolated values
%
% 2013, Sept, Version 3, minor update to new matlab version


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Morteza Bahram -Ph.D in Analytical Chemistry 
% Department of chemistry, 
% Faculty of science, 
% Urmia university- 
% Urmia - Iran
%
% morteza.bahram@gmail.com
% 
% &
% Rasmus Bro - rb@life.ku.dk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LowZero = 10;

if nargin<1
  error('At least one input must be given in EEMSCAT');
end


if nargin<2
  MissRayleh =[25 25];
end
if nargin<3
  MissRaman =[0 0];
end
if nargin<4
  MissRayleh2 =[0 0];
end


tic

r1=MissRayleh(1);
r2=MissRayleh(2);
r3=LowZero;
r4=MissRaman(1);
r5=MissRaman(2);
r6=MissRayleh2(1);
r7=MissRayleh2(2);
ax=X.axisscale{2};    %%%%  emission axcisscale
ax2=X.axisscale{3};  %%%%  excitation axcisscale

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q=isnan(X.data);  %%% NaN element at the original data due to emission out
%%% of spectroscopic instrumental range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xmiss1=X;   %%  initialize Xmiss1 for Rayleh scatter
%      Rayleh removal
for i=1:size(X.data,3);
  j = find((ax<(ax2(i)+r2))& (ax>(ax2(i)-r1)));
  Xmiss1.data(:,j,i)=NaN;
end
EEM_correct_Rayleh=Xmiss1.data;  %%INITIALIZING  EEM_correct

%%%%%%%%%%%%%%   Raman removal
wavenumber_forraman = 3600; % Possibly correct for other solvents than water
for k = 1:size(X.data,3) % EVERY EXCITATION
  Emission(k)=1e7/( (1e7/ax2(k))-wavenumber_forraman);
end
Xmiss2=X;    %%initialize Xmiss1 for Raman scatter
for k=1:length(Emission);
  j=find((ax>Emission(k)-r5)& (ax<Emission(k)+r4));
  Xmiss2.data(:,j,k)=NaN;
end

%      Second order Rayleh removal
Xmiss3=X;   %%initialize Xmiss1 for Rayleh2 scatter
for i=1:size(X.data,3);
  j = find((ax>((2*ax2(i))-r7))& (ax<((2*ax2(i))+r6)));
  Xmiss3.data(:,j,i)=NaN;
end
EEM_correct_Rayleh2=Xmiss3.data;   %%%% only Rayleh2 scatering was removed here

%%%%%%%%%%%%     Rayleh and Raman removal

for k=1:length(Emission);
  j=find((ax>Emission(k)-r5)& (ax<Emission(k)+r4));
  Xmiss1.data(:,j,k)=NaN;
end
EEM_Cutted_rr=Xmiss1.data;     %%%%%%%%%  Only  rayle 1 and raman was removed here
for i=1:size(X.data,3);
  j = find((ax>((2*ax2(i))-r7))& (ax<((2*ax2(i))+r6)));
  Xmiss1.data(:,j,i)=NaN;
end

EEM_Cutted_rrr2=Xmiss1.data;   %%% Rayleh 1 and rayleh2 and raman was removed



%%INITIALIZING  EEM_correct
EEM_correct=EEM_Cutted_rrr2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%  Interpolating the latest line of emission at every sample

for j=1:size(X.data,1);  % LOPE FOR EVERY SAMPLE
  ppp=squeeze(EEM_correct(j,end,:));  %%%%%%% the last EMISSION vector in every sample
  w= ~isnan(ppp);
  ww=find(w==1);
  Eendcut=ppp(ww);  %%% the last CUTTED EMISSION vector in every sample

  ax2cut=ax2(ww);      %%% Cutted EMISSION AXIS
  mm=interp1(ax2cut,Eendcut,ax2,'pchip','extrap'); %%% interpolation using cubic option

  ffff=find(mm<0);
  mm(ffff)=0;   %% Forcing Non negativity at the interpolation

  EEM_correct(j,end,:)=mm;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Interpolation  %%%%%%%%%%%


for j=1:size(X.data,1);  % LOOP FOR EVERY SAMPLE
  for i=1:length(ax2);  % LOOP FOR EVRY EXCITATION WAVELENGHT
    dd=EEM_correct(j,:,i);
    if ~isnan(dd(1));
      gg= ~isnan(dd);
      ff=find(gg==1);
      Ecut=dd(ff);       %%% CUTTED EMISSION
      axcut=ax(ff);      %%% Cutted EMISSION AXIS

      p=interp1(axcut,Ecut,ax,'pchip'); %%% interpolation using cubic option
      EEM_correct(j,:,i)=p;

    elseif ~isnan(dd(1));
      gg= ~isnan(dd);
      ff=find(gg==1);
      Ecut=dd(ff);       %%% CUTTED EMISSION
      axcut=ax(ff);      %%% Cutted EMISSION AXIS

      a2=ax2(i)-r3; a3=a2-5;%% calculation the emisions
      %%%%wavelengh r3 nm below excitation wavelengh
      Ecut2=[0,Ecut];    %%% forcing zero for emission 30 nm below excitation wavelenght
      axcut2=[a3,axcut];  %%% axis scale scale of emission vector
      p=interp1(axcut2,Ecut2,ax,'pchip'); %%% interpolation using spline option
      EEM_correct(j,:,i)=p;


    else

      a1=ax2(i)-r3-5;a2=ax2(i)-r3; %% calculation the emisions
      %%%%wavelengh 30 nm below excitation wavelengh
      gg= ~isnan(dd);
      ff=find(gg==1);
      Ecut=dd(ff);       %%% CUTTED EMISSION
      axcut=ax(ff);      %%% Cutted EMISSION AXIS

      Ecut2=[0,0,Ecut];    %%% forcing zero for emission 30 nm below excitation wavelenght
      axcut2=[a2,a1,axcut];  %%% axis scale scale of emission vector
      p=interp1(axcut2,Ecut2,ax,'pchip'); %%% interpolation using cubic option
      EEM_correct(j,:,i)=p;

    end
  end
end

Data_Correct=EEM_correct;
Data_Correct(Q)=NaN;
data = Data_Correct;
Data_Correct=dataset(data);
Data_Correct.axisscale{2}=ax;
Data_Correct.axisscale{3}=ax2;
Data_Correct.name=('INTRRR2(Cubic)');

% inc=X.include;
% Data_Correct.include{1}=inc{1};
% Data_Correct.include{2}=inc{2};
% Data_Correct.include{3}=inc{3};


EEM_Cut=EEM_Cutted_rrr2;
EEM_Cut(Q)=NaN;
EEM_Cut=dataset(EEM_Cut);
EEM_Cut.axisscale{2}=ax;
EEM_Cut.axisscale{3}=ax2;
EEM_Cut.name=('Missed Data');

% inc=X.include;
% EEM_Cut.include{1}=inc{1};
% EEM_Cut.include{2}=inc{2};
% EEM_Cut.include{3}=inc{3};

t=toc;

subplot(1,3,1)
mesh(ax2,ax,squeeze(X.data(1,:,:)))
title ('first sample - raw data')
axis tight

subplot(1,3,2)
mesh(ax2,ax,squeeze(EEM_Cutted_rrr2(1,:,:)))
title ('first sample - scattering removed')
axis tight

subplot(1,3,3)
mesh(ax2,ax,squeeze(Data_Correct.data(1,:,:)))
title ('first sample - interpolated')
axis tight
