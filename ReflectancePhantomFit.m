function [q,RefFit,Wav,AbsFit,ScatFit] = ReflectancePhantomFit(InputRef,InputWav,varargin)
%  Copyright (C) 2022, A.Walter (personal alec.b.walter@vanderbilt.edu) MIT license
% Summary of this function:
%   varargin takes in two options, and PlotResults and MinScattering, as
%   name-value pairs. It can also take in specified bands and weights if
%   desired. If so, they must come before the options.
%
%   Bands should be notated by an array of wavelengths in nm (eg [540:610])
%   and its corresponding weight, usually as a power of 10 (but not
%   required).
%
%   'PlotResults' value of 1 will plot the target reflectance and the
%   reflectance of the fit, while highlighting the selected bands. A value
%   of 2 will also plot the absorption and scattering properties of the
%   phantom fit. A value of zero will supress all plots and is default.

%   'MinScattering' indicates a desired minimum that the reduced
%   scattering of the diffuse reflectance phantom at 950nm must have in 1/mm. This
%   is used to prevent unrealsitic Absorption-Scattering pairs that result
%   in the same reflectance. The default value is 0.5 1/mm.


%% Check for options
optioncount=0;
if any(strcmp(varargin,'MinScattering'))
    temploc=find(strcmp(varargin,'MinScattering')==1);
    SolveSeperate=varargin{temploc+1};
    optioncount=optioncount+2;
else
    MinScattering=0.5;
end
if any(strcmp(varargin,'PlotResults'))
    temploc=find(strcmp(varargin,'PlotResults')==1);
    PlotResults=varargin{temploc+1};
    optioncount=optioncount+2;
else
    PlotResults=0;
end

NumVars=length(varargin);
if rem(NumVars,2)>0
    error('An input is missing its partner');
end

if NumVars>optioncount
    for i=1:(NumVars-optioncount)/2
    Bands.(['a' num2str(i)])=varargin{2*i-1};
    Weights.(['a' num2str(i)])=varargin{2*i};
    end
    UseBands=1;
else
    UseBands=0;
end

%% Load Nomalized Pigment Properties
temp_a=readtable([pwd, '\NormalizedPigmentProperties\Absorption.csv']); %open normalized absorption file located in current directory
temp_s=readtable([pwd, '\NormalizedPigmentProperties\ReducedScattering.csv']); %open normalized reduced scattering file located in current directory

wav=temp_a{:,1}; %seperate wavelengths

if min(wav)<min(InputWav)
    BlueDiff=min(InputWav)-min(wav);
else
    BlueDiff=0;
end

if max(wav)>max(InputWav)
    RedDiff=max(wav)-max(InputWav);
else
    RedDiff=0;
end

wav=temp_a{1+RedDiff:end-BlueDiff,1};
Absorption=temp_a{1+RedDiff:end-BlueDiff,2:end}; %seperate absorption data
Scattering=temp_s{1+RedDiff:end-BlueDiff,2:end}; %seperate scattering data

%% Bring target reflectance into correct wavelength-space
TargetRef=interp1(InputWav,InputRef,wav);

%% Set weights
if UseBands==1
    w=zeros(length(wav),1);
    for k=1:length(fieldnames(Bands))
        [i,j]=find(wav==Bands.(['a' num2str(k)]));
        w(i)=Weights.(['a' num2str(k)])(1);
    end
    w2=diag(sqrt(w));
else
    w=ones(length(wav),1);
    w2=diag(sqrt(w));
end
%% Solve
fun = @(q)double(w2*TargetRef-w2*WMC_model([Absorption(:,:)*q,Scattering(:,:)*q],1.556,[0])); %Minimization function
fun2= @(q)double(WMC_model([Absorption(:,:)*q,Scattering(:,:)*q],1.556,[0])); %Gives Reflectance given pig. conc.
options = optimoptions(@lsqnonlin,'Display','iter');
EndScat=MinScattering./2;
q0=zeros(size(Absorption,2),1)+1;
qmax=q0*10000;
qmax(end)=20; %limit the amount of Al2O3 used
qmin=zeros(size(Absorption,2),1);
q=qmin;
while EndScat<MinScattering
    qmin(end-2:end-1)=q(end-2:end-1);
    q0=q.*MinScattering./EndScat;q0(end)=10;
    disp('This may take a bit...')
    q=lsqnonlin(fun,q0,qmin,qmax,options);

    AbsFit=Absorption*q;
    ScatFit=Scattering*q;
    EndScat=min(ScatFit)
end
Wav=wav;
RefFit=fun2(q);

%% Plot
if PlotResults>0
    clr=[0.6667	0.2	0.4667
        0.9333	0.4	0.4667
        0.8	0.7333	0.2667
        0.1333	0.5333	0.2
        0.4	0.8	0.9333
        0.2667	0.4667	0.6667];
    figure();
    colororder(clr)
    plot(wav,TargetRef,'-');
    hold on;plot(wav,RefFit,'--');
    if UseBands==1
    for k=1:length(fieldnames(Bands))
        r1=rectangle('Position',[min(Bands.(['a' num2str(k)])) 0 max(Bands.(['a' num2str(k)]))-min(Bands.(['a' num2str(k)])) 1],'FaceColor',[0,0,0,.2],'EdgeColor','none');
        uistack(r1,'bottom')
    end
    end
    set(gca, 'Layer', 'top')
    xlim([370 950])
    ylim([0 1])
    ax=gca;
    ax.FontSize = 12;
    xlabel('Wavelength (nm)','Fontsize',15)
    ylabel('Diffuse Reflectance','Fontsize',15)
    legend('Target','Phantom Fit')
    if PlotResults>1
    figure();
    colororder(clr)
    plot(wav,AbsFit,'-');
    hold on;plot(wav,ScatFit,'-');
    set(gca, 'YScale', 'log')
    if UseBands==1
    for k=1:length(fieldnames(Bands))
        r1=rectangle('Position',[min(Bands.(['a' num2str(k)])) 0 max(Bands.(['a' num2str(k)]))-min(Bands.(['a' num2str(k)])) 1],'FaceColor',[0,0,0,.2],'EdgeColor','none');
        uistack(r1,'bottom')
    end
    end
    set(gca, 'Layer', 'top')
    xlim([370 950])
    ylim([min(min(AbsFit),min(ScatFit))./2 max(max(AbsFit),max(ScatFit)).*2])
    ax=gca;
    ax.FontSize = 12;
    xlabel('Wavelength (nm)','Fontsize',15)
    ylabel('Absorption / Reduced Scattering Coefficient (mm^-^1)','Fontsize',15)
    legend('Phantom Fit Absorption','Phantom Fit Scattering')

end
end

end

function [R_of_fx]= WMC_model(op,n,fx)
% white MC model for light transport in tissue  
load([pwd, '\WhiteMonteCarlo\WMC_10S_DoubleLog.mat']);musr=10;DvT2=DvT2';
v=300/n;

R=DvT2'./1e7;
R_of_rho_and_t=R;
t=b;
delta_t=t(2)-t(1);
rho=a;
delta_rho=rho(2)-rho(1);
 
fresn=((1-n)/(1+n))^2;
 
% for each pair of optical properties (i.e. each wavelength)
mua=op(:,1);
musp=op(:,2);
 

for i=1:length(mua) %loop over wavelengths
R_of_rho = sum(R_of_rho_and_t.*repmat(exp(-mua(i)*v*t*(musr/musp(i))).* ... 
    (musr/musp(i)),length(rho),1), 2)*(musr/musp(i))^-2; %reflectance per m^2 per s, 

f=length(fx); %number of spatial frequencies used in measurement
r=length(rho); % number of spatial points used in Monte Carlo simulation
 
R_of_fx(i,:) = sum( repmat(R_of_rho',[f,1]) .* besselj(0,repmat(2*pi*fx',[1,r]).*repmat((rho).*(musr/musp(i)),[f,1])) .*  repmat((musr/musp(i)),[f,r]), 2)';
end
end