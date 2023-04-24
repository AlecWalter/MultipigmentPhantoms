function [q,AbsFit,ScatFit,Wav] = MultibandPhantomFit(InputAbs,InputScat,InputWav,varargin)
%  Copyright (C) 2022, A.Walter (personal alec.b.walter@vanderbilt.edu) MIT license
%Summary of this function:
%   varargin takes in the specified bands (with paired weight) and three
%   options, SolveSeperate, AbsPigScattering and PlotResults, as
%   name-value pairs. The options must come last.
%
%   Bands should be notated by an array of wavelengths in nm (eg [540:610])
%   and its corresponding set of weights in the format
%   [W_abs,W_sct],usually as a power of 10 (but not required). Weights
%   should be adjusted to acheive the desired result. Your mileage may
%   vary, but I have found having a scattering weight around 10^4 to 10^6
%   less than the corresponding absorption weight to yield good results.
%
%   'SolveSeperate' value of 1 if Absoption sould be solved for first using
%   only the absorbing pigments (the first 17 included pigments) and 0 if
%   they should be solved together (accounts for white pigment absorption
%   in the UV). Default value is here is 0, as UV absorption only comes
%   into play if you need that band, in which case solving together will
%   yield better results.
%
%   'AbsPigScattering' value of 1 takes into account the scattering
%   provided by the absorbing pigments, 0 only uses the scattering of the
%   white pigments. Default value is 1, though ideally the contribution
%   should be small. At high aborption values, increased variation from fit
%   to phantom should be expected.
%   
%   'PlotResults' value of 1 will plot the target properties and the
%   properties of the fit, while highlighting the selected bands. A value
%   of 0 will supress the plot and is default. A value of 2 will plot the
%   target and predicted properties while not highlighting the selected
%   bands.
%% Check varargin, grab band information
optioncount=0;
if any(strcmp(varargin,'SolveSeperate'))
    temploc=find(strcmp(varargin,'SolveSeperate')==1);
    SolveSeperate=varargin{temploc+1};
    optioncount=optioncount+2;
else
    SolveSeperate=0;
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

for i=1:(NumVars-optioncount)/2
    Bands.(['a' num2str(i)])=varargin{2*i-1};
    Weights.(['a' num2str(i)])=varargin{2*i};
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

%% Set Band Weights
w=zeros(length(wav),1); %Absorption weights
w3=zeros(length(wav),1); %Scattering weights
for k=1:length(fieldnames(Bands))
    [i,j]=find(wav==Bands.(['a' num2str(k)]));
    w(i)=Weights.(['a' num2str(k)])(1);
    w3(i)=Weights.(['a' num2str(k)])(2);
end
w2=diag(sqrt(w));
w4=diag(sqrt(w3));
%% Solve
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
    if AbsPigScattering == 1 %Check if the scattering of the absorping pigments should be taken into account
        fun = @(q)sum(((w2*TargetAbs-w2*Absorption*q)).^2+(w4*TargetScat-w4*Scattering*q).^2);
    else
        fun = @(q)sum(((w2*TargetAbs-w2*Absorption*q)).^2+(w4*TargetScat-w4*Scattering(:,end-2:end)*q(end-2:end)).^2);
    end
    options = optimoptions(@fmincon,'Display','none');
    q0=zeros(size(Absorption,2),1)+.1;
    qmin=zeros(size(Absorption,2),1);
    qmax=q0*10000;
    qmax(end)=5; %limit the amount of Al2O3 used
    q=fmincon(fun,q0,[],[],[],[],qmin,qmax,[],options);

end

AbsFit=Absorption*q;
ScatFit=Scattering*q;
Wav=wav;
%% Plots
if PlotResults>0
    clr=[0.6667	0.2	0.4667
        0.1333	0.5333	0.2
        0.9333	0.4	0.4667
        0.8	0.7333	0.2667
        0.4	0.8	0.9333
        0.2667	0.4667	0.6667];
    figure();
    colororder(clr)
    plot(wav,TargetAbs,'-','LineWidth',2.25);
    hold on;plot(wav,AbsFit,':','LineWidth',2.25);
    set(gca, 'YScale', 'log')
    if PlotResults==1
    for k=1:length(fieldnames(Bands))
        r1=rectangle('Position',[min(Bands.(['a' num2str(k)])) min(min(AbsFit),min(TargetAbs))./2 max(Bands.(['a' num2str(k)]))-min(Bands.(['a' num2str(k)])) max(max(AbsFit),max(TargetAbs)).*2-min(min(AbsFit),min(TargetAbs))./2],'FaceColor',[0,0,0,.2],'EdgeColor','none');
        uistack(r1,'bottom')
    end
    end
    set(gca, 'Layer', 'top')
    xlim([min(wav) max(wav)])
    ylim([min(min(AbsFit),min(TargetAbs))./2 max(max(AbsFit),max(TargetAbs)).*2])
    ax=gca;
    ax.FontSize = 12;
    xlabel('Wavelength (nm)','Fontsize',15)
    ylabel('Absorption Coefficient (mm^-^1)','Fontsize',15)
    legend('Target','Phantom Fit')

    figure();
    colororder(clr)
    plot(wav,TargetScat,'-','LineWidth',2.25);
    hold on;plot(wav,ScatFit,':','LineWidth',2.25);
    if PlotResults==1
    for k=1:length(fieldnames(Bands))
        r1=rectangle('Position',[min(Bands.(['a' num2str(k)])) 0 max(Bands.(['a' num2str(k)]))-min(Bands.(['a' num2str(k)])) max(max(ScatFit),max(TargetScat)).*1.2-0],'FaceColor',[0,0,0,.2],'EdgeColor','none');
        uistack(r1,'bottom')
    end
    end
    set(gca, 'Layer', 'top')
    xlim([min(wav) max(wav)])
    ylim([0 max(max(ScatFit),max(TargetScat)).*1.2])
    ax=gca;
    ax.FontSize = 12;
    xlabel('Wavelength (nm)','Fontsize',15)
    ylabel('Reduced Scattering Coefficient (mm^-^1)','Fontsize',15)
    legend('Target','Phantom Fit')
end

end