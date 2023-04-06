function [q,AbsFit,ScatFit,Wav] = BroadbandPhantomFit(InputAbs,InputScat,InputWav,varargin)
%  Copyright (C) 2022, A.Walter (personal alec.b.walter@vanderbilt.edu) MIT license
%Summary of this function goes here
%   Detailed explanation goes here
%   varargin takes in three options, SolveSeperate, AbsPigScattering and PlotResults, as
%   name-value pairs.
%
%   'SolveSeperate' value of 1 if Absoption sould be solved for first using
%   only the absorbing pigments (the first 17 included pigments) and 0 if
%   they should be solved together (accounts for white pigment absorption
%   in the UV). Default value is 1, as I find it more useful to ignore the
%   UV absorption unless that wavelength range is specifically needed.
%
%   'AbsPigScattering' value of 1 takes into account the scattering
%   provided by the absorbing pigments, 0 only uses the scattering of the
%   white pigments. Default value is 1, though ideally the contribution
%   should be small. At high aborption values, increased variation from fit
%   to phantom should be expected.
%
%   'PlotResults' value of 1 will plot the target properties and the
%   properties of the fit, while highlighting the selected bands. A value
%   of 0 will supress the plot and is default.

%% Check for options
optioncount=0;
if any(strcmp(varargin,'SolveSeperate'))
    temploc=find(strcmp(varargin,'SolveSeperate')==1);
    SolveSeperate=varargin{temploc+1};
    optioncount=optioncount+2;
else
    SolveSeperate=1;
end
if any(strcmp(varargin,'AbsPigScattering'))
    temploc=find(strcmp(varargin,'AbsPigScattering')==1);
    AbsPigScattering=varargin{temploc+1};
    optioncount=optioncount+2;
else
    AbsPigScattering=1;
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

%% Bring target properties into correct wavelength-space
TargetAbs=interp1(InputWav,InputAbs,wav);
TargetScat=interp1(InputWav,InputScat,wav);

%% Weight and Solve
w=max(TargetAbs)./TargetAbs; %create weight-function to correct for log-range of absorption
w2=diag(sqrt(w));
w3=TargetScat./TargetScat; %placeholder scattering weight function
% w3=max(TargetScat)./TargetScat; %uncomment to create weight-function to correct for range of scattering
w4=diag(sqrt(w3));
if SolveSeperate == 1 %Check if absorption should be solved before scattering
    x=lsqnonneg(w2*Absorption(:,1:end-3),w2*TargetAbs); %solve for concentrations of absorbing pigments

    if AbsPigScattering == 1 %Check if the scattering of the absorping pigments should be taken into account
        z=lsqnonneg(w4*Scattering(:,end-2:end),w4*(TargetScat-Scattering(:,1:end-3)*x));
        q=[x;z];
    else
        z=lsqnonneg(w4*Scattering(:,end-2:end),w4*TargetScat);
        q=[x;z];
    end
else
    q0=zeros(size(Absorption,2),1)+.1;
    qmin=zeros(size(Absorption,2),1);
    qmax=q0*10000;
    w=max(TargetAbs)./TargetAbs; %create weight-function to correct for range of absorption
    w2=diag(sqrt(w));
    w3=TargetScat./TargetScat; %placeholder scattering weight function
%     w3=max(TargetScat)./TargetScat; %uncomment to create weight-function to correct for range of scattering
    w4=diag(sqrt(w3));
    if AbsPigScattering == 1 %Check if the scattering of the absorping pigments should be taken into account
        fun = @(q)sum(((w2*TargetAbs-w2*Absorption*q)).^2+(w4*TargetScat-w4*Scattering*q).^2);
    else
        fun = @(q)sum(((w2*TargetAbs-w2*Absorption*q)).^2+(w4*TargetScat-w4*Scattering(:,end-2:end)*q(end-2:end)).^2);
    end
    options = optimoptions(@fmincon,'Display','none');
    q0=zeros(size(Absorption,2),1)+.1;
    qmin=zeros(size(Absorption,2),1);
    qmax=q0*10000;
%     qmax(end)=20; %limit the amount of Al2O3 used
    q=fmincon(fun,q0,[],[],[],[],qmin,qmax,[],options);

end

AbsFit=Absorption*q;
ScatFit=Scattering*q;
Wav=wav;

%% Plot
if PlotResults==1
    clr=[0.6667	0.2	0.4667
        0.9333	0.4	0.4667
        0.8	0.7333	0.2667
        0.1333	0.5333	0.2
        0.4	0.8	0.9333
        0.2667	0.4667	0.6667];
    figure();
    colororder(clr)
    plot(wav,TargetAbs,'-');
    hold on;plot(wav,AbsFit,'--');
    set(gca, 'YScale', 'log')
    set(gca, 'Layer', 'top')
    xlim([370 950])
%     ylim([min(min(AbsFit),min(TargetAbs))./2 max(max(AbsFit),max(TargetAbs)).*2])
    ax=gca;
    ax.FontSize = 12;
    xlabel('Wavelength (nm)','Fontsize',15)
    ylabel('Absorption Coefficient (mm^-^1)','Fontsize',15)
    legend('Target','Phantom Fit')

    figure();
    colororder(clr)
    plot(wav,TargetScat,'-');
    hold on;plot(wav,ScatFit,'--');
    set(gca, 'Layer', 'top')
    xlim([370 950])
%     ylim([0 max(max(ScatFit),max(TargetScat)).*1.2])
    ax=gca;
    ax.FontSize = 12;
    xlabel('Wavelength (nm)','Fontsize',15)
    ylabel('Reduced Scattering Coefficient (mm^-^1)','Fontsize',15)
    legend('Target','Phantom Fit')
end

end