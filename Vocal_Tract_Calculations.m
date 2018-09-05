% first specify name = 'filename.mat';
% load(name);

name = 'Example_MRI.mat';
load(name)

% for i=1:48
%     SG.radius(i) = 0.0075
% end

VT.length = 0.17;
VT.radius = 0.013;
VT.alpha_multiplier = 5;
glottis.length = 0.01;
glottis.radius = 0.00000001;
glottis.alpha_multiplier = 5;
SG.length = 0.17;

SG.radius = 0.013;
SG.alpha_multiplier = 5;


%% run the combined processing

% result=Vocal_Tract_General_Model_2017_SVT(tube);
% if glottis and SG are specified then use
% result=Vocal_Tract_General_Model_2016(tube, glottis, SG);

result=Vocal_Tract_General_Model_2017_SVT(VT, glottis, SG);

%% Figure 1 Impedance
fg(1) = figure('Name',strcat(name, ' Rigid'),'NumberTitle','off');

% set(gcf,'outerposition',get(0,' screensize'));
ax(1) = subplot(4,1,1);
% plot([-result.tube.radius_fn(1,1); ((result.tube.radius_fn(:,1)...
%     -result.tube.radius_fn(1,1))-result.tube.glottis.length)],...
%     [result.tube.radius_fn(1,2); result.tube.radius_fn(:,2)],...
%     '-', 'LineWidth',2); 

% % adjust profile to glottal 0
% distance = result.tube.radius_fn(:,1)-sum(result.tube.SG.length)...
%     -sum(result.tube.glottis.length);

plot([-result.tube.radius_fn(1,1); result.tube.distance],...
    [result.tube.radius_fn(1,2); result.tube.radius_fn(:,2)],...
    '-', 'LineWidth',2); 
% if SG and glottis tube segments are > 1
% plot(distance, result.tube.radius_fn(:,2), '-', 'LineWidth',2);
xlabel('Distance from glottis (m)');
ylabel('Radius (m)');
title('Airway profile');
hold all;

% Z at lips
ax(2) =subplot(4,2,5);
% semilogy(result.freq,abs(result.tube.Z_in(end,:)), '-', 'LineWidth',2);
semilogy(result.freq,abs(result.tube.in.Z_flex), '-', 'LineWidth',2);
axis([0 4000 10^4 10^8])
% xlabel('frequency(Hz)');
ylabel('|Z| (Pa\bullets\bulletm^{-3})');
title('Z at the lips');
set(gca,'xticklabel',[]);
xlabel('Frequency (Hz)')
hold all;

% Phase at lips
ax(3) = subplot(4,2,6);
% plot(result.freq,angle(result.tube.Z_in(end,:))*180/pi, 'LineWidth',2);
plot(result.freq,angle(result.tube.in.Z)*180/pi, 'LineWidth',2);
% xlabel('frequency(Hz)');
ylabel('Angle(Deg)');
title('Phase at the lips');
set(gca,'xticklabel',[]);
hold all;
 
%Z at glottis
ax(4) =subplot(4,2,3);
% semilogy(result.freq,abs(result.tube.Z(1,:)) ,'-', 'LineWidth',2);
semilogy(result.freq,abs(result.tube.out.Z) ,'-', 'LineWidth',2);
% semilogy(result.freq,abs(result.tube.Z(2,:)) ,'-', 'LineWidth',2);
% xlabel('frequency(Hz)');
ylabel('Z (Pa\bullets\bulletm^{-3})');
title('Z at the glottis');
set(gca,'xticklabel',[]);
hold all;

%Phase at glottis
ax(5) =subplot(4,2,4);
% plot(result.freq,angle(result.tube.Z(1,:))*180/pi, 'LineWidth',2);
plot(result.freq,angle(result.tube.out.Z)*180/pi, 'LineWidth',2);
% xlabel('frequency(Hz)');
ylabel('Angle (Deg)');
title('Phase at the glottis');
set(gca,'xticklabel',[]);
hold all;

ax(6) =subplot(4,2,7);
plot(result.freq, 20.*log10(abs(result.tube.in.Z_ratio_1)),...
    '-', 'LineWidth',2);  hold all;
% plot(result.freq, 20.*log10(abs(result.tube.in.Z_ratio_2)),...
%     '-', 'LineWidth',2); 
% plot(result.freq, 20.*log10(abs(result.tube.in.Z_ratio_3)),...
%     '-', 'LineWidth',2); 
% plot(result.freq, 20.*log10(abs(result.tube.in.Z_ratio_4)),...
%     '-', 'LineWidth',2); 
xlabel('Frequency (Hz)');
ylabel('dB');
title('Z ratio at lips');
hold all;

ax(7) =subplot(4,2,8);
plot(result.freq,angle(result.tube.in.Z_ratio_1)*180/pi,...
    'LineWidth',2); hold all;
% plot(result.freq,angle(result.tube.in.Z_ratio_2)*180/pi,...
%     'LineWidth',2);
% plot(result.freq,angle(result.tube.in.Z_ratio_3)*180/pi,...
%     'LineWidth',2);
% plot(result.freq,angle(result.tube.in.Z_ratio_4)*180/pi,...
%     'LineWidth',2);
xlabel('frequency (Hz)');
ylabel('Angle (Deg)');
title('Phase of Z ratio at lips');
hold all;

% link all x axes to enable synchronised zooming
linkaxes(ax(2:7),'x');


%% Figure 2 Transfer Functions
fg(2) = figure('Name',strcat(name, ' Transfer Functions'),'NumberTitle','off');

% set(gcf,'outerposition',get(0,' screensize'));
bx(1) = subplot(4,1,1);
plot([-result.tube.radius_fn(1,1); result.tube.distance],...
    [result.tube.radius_fn(1,2); result.tube.radius_fn(:,2)],...
    '-', 'LineWidth',2); 
xlabel('Distance from glottis (m)');
ylabel('Radius (m)');
title('Airway profile');
hold all;

% T_PP
bx(2) =subplot(4,2,5);
plot(result.freq, 20.*log10(abs(result.tube.out.T_PP)),...
    '-', 'LineWidth',2);
% xlabel('frequency(Hz)');
ylabel('(dB)');
title('Power out/in');
set(gca,'xticklabel',[]);
hold all;

% T_Pp phase
bx(3) = subplot(4,2,6);
plot(result.freq, 20.*log10(abs(result.tube.out.T_Pp)),...
    'LineWidth',2);
% xlabel('frequency(Hz)');
ylabel('(dB)');
title('Power out/ Pressure in');
set(gca,'xticklabel',[]);
hold all;
 
%T_pp
bx(4) =subplot(4,2,3);
plot(result.freq, 20.*log10(abs(result.tube.out.T_pp)),...
    '-', 'LineWidth',2);
% xlabel('frequency(Hz)');
ylabel('(dB)');
title('Pressure out/in');
set(gca,'xticklabel',[]);
hold all;

%T_uu phase
bx(5) =subplot(4,2,4);
% plot(result.freq,angle(result.tube.Z(1,:))*180/pi, 'LineWidth',2);
plot(result.freq, 20.*log10(abs(result.tube.out.T_uu)),...
    '-', 'LineWidth',2);
% xlabel('frequency(Hz)');
ylabel('(dB)');
title('Flow out/in');
set(gca,'xticklabel',[]);
hold all;

% T_pu
bx(6) =subplot(4,2,7);
plot(result.freq, 20.*log10(abs(result.tube.out.T_pu)),...
    '-', 'LineWidth',2); 
xlabel('Frequency (Hz)');
ylabel('dB');
title('Transpedance pressure out / flow in');
hold all;

% T_Pu
bx(7) =subplot(4,2,8);
plot(result.freq, 20.*log10(abs(result.tube.out.T_Pu)),...
    '-', 'LineWidth',2);
xlabel('frequency (Hz)');
ylabel('(dB)');
title('Power out / flow in');
hold all;

% link all x axes to enable synchronised zooming
linkaxes(bx(2:7),'x');


%%
%% Figure 3 Impedance with non-rigid walls
fg(3) = figure('Name',strcat(name, ' Non-rigid'),'NumberTitle','off');

% set(gcf,'outerposition',get(0,' screensize'));
cx(1) = subplot(4,1,1);
plot([-result.tube.radius_fn(1,1); result.tube.distance],...
    [result.tube.radius_fn(1,2); result.tube.radius_fn(:,2)],...
    '-', 'LineWidth',2); 
xlabel('Distance from glottis (m)');
ylabel('Radius (m)');
title('Airway profile');
hold all;

% Z at lips
cx(2) =subplot(4,2,5);
% semilogy(result.freq,abs(result.tube.Z_in(end,:)), '-', 'LineWidth',2);
semilogy(result.freq,abs(result.tube.in.Z_flex), '-', 'LineWidth',2);
% xlabel('frequency(Hz)');
ylabel('Z (Pa\bullets\bulletm^{-3})');
title('Z at the lips');
set(gca,'xticklabel',[]);
hold all;

% Phase at lips
cx(3) = subplot(4,2,6);
% plot(result.freq,angle(result.tube.Z_in(end,:))*180/pi, 'LineWidth',2);
plot(result.freq,angle(result.tube.in.Z_flex)*180/pi, 'LineWidth',2);
% xlabel('frequency(Hz)');
ylabel('Angle(Deg)');
title('Phase at the lips');
set(gca,'xticklabel',[]);
hold all;
 
%Z at glottis
cx(4) =subplot(4,2,3);
% semilogy(result.freq,abs(result.tube.Z(1,:)) ,'-', 'LineWidth',2);
semilogy(result.freq,abs(result.tube.out.Z_flex) ,'-', 'LineWidth',2);
% semilogy(result.freq,abs(result.tube.Z(2,:)) ,'-', 'LineWidth',2);
% xlabel('frequency(Hz)');
ylabel('Z (Pa\bullets\bulletm^{-3})');
title('Z at the glottis');
set(gca,'xticklabel',[]);
hold all;

%Phase at glottis
cx(5) =subplot(4,2,4);
% plot(result.freq,angle(result.tube.Z(1,:))*180/pi, 'LineWidth',2);
plot(result.freq,angle(result.tube.out.Z_flex)*180/pi, 'LineWidth',2);
% xlabel('frequency(Hz)');
ylabel('Angle (Deg)');
title('Phase at the glottis');
set(gca,'xticklabel',[]);
hold all;

%ratio
Zparallel_3 = (result.tube.Z_rad(end,:).*result.tube.in.Z_flex)./...
    (result.tube.Z_rad(end,:) + result.tube.in.Z_flex);

z_flex_ratio_at_the_lips = ...
    Zparallel_3./(result.tube.Z_rad(end,:));

cx(6) =subplot(4,2,7);
plot(result.freq, 20.*log10(abs(z_flex_ratio_at_the_lips)), '-', 'LineWidth',2); 
xlabel('Frequency (Hz)');
ylabel('dB');
title('Z ratio at lips');
hold all;

cx(7) =subplot(4,2,8);
plot(result.freq,angle(z_flex_ratio_at_the_lips)*180/pi, 'LineWidth',2);
xlabel('frequency (Hz)');
ylabel('Angle (Deg)');
title('Phase of Z ratio at lips');
hold all;

% link all x axes to enable synchronised zooming
linkaxes(cx(2:7),'x');


%% Figure 4 Transfer Functions Non-rigid
fg(4) = figure('Name',strcat(name, ' Non-rigid Transfer Functions'),...
    'NumberTitle','off');

% set(gcf,'outerposition',get(0,' screensize'));
dx(1) = subplot(4,1,1);
plot([-result.tube.radius_fn(1,1); result.tube.distance],...
    [result.tube.radius_fn(1,2); result.tube.radius_fn(:,2)],...
    '-', 'LineWidth',2); 
xlabel('Distance from glottis (m)');
ylabel('Radius (m)');
title('Airway profile');
hold all;

% T_PP
dx(2) =subplot(4,2,5);
plot(result.freq, 20.*log10(abs(result.tube.out.T_PP_flex)),...
    '-', 'LineWidth',2);
% xlabel('frequency(Hz)');
ylabel('(dB)');
title('Power out/in');
set(gca,'xticklabel',[]);
hold all;

% T_Pp phase
dx(3) = subplot(4,2,6);
plot(result.freq, 20.*log10(abs(result.tube.out.T_Pp_flex)),...
    'LineWidth',2);
% xlabel('frequency(Hz)');
ylabel('(dB)');
title('Power out/ Pressure in');
set(gca,'xticklabel',[]);
hold all;
 
%T_pp
dx(4) =subplot(4,2,3);
plot(result.freq, 20.*log10(abs(result.tube.out.T_pp_flex)),...
    '-', 'LineWidth',2);
% xlabel('frequency(Hz)');
ylabel('(dB)');
title('Pressure out/in');
set(gca,'xticklabel',[]);
hold all;

%T_uu phase
dx(5) =subplot(4,2,4);
% plot(result.freq,angle(result.tube.Z(1,:))*180/pi, 'LineWidth',2);
plot(result.freq, 20.*log10(abs(result.tube.out.T_uu_flex)),...
    'LineWidth',2);
% xlabel('frequency(Hz)');
ylabel('Angle (Deg)');
title('Flow out/in');
set(gca,'xticklabel',[]);
hold all;

% T_pu
dx(6) =subplot(4,2,7);
plot(result.freq, 20.*log10(abs(result.tube.out.T_pu_flex)),...
    '-', 'LineWidth',2); 
xlabel('Frequency (Hz)');
ylabel('dB');
title('Transpedance pressure out / flow in');
hold all;

% T_Pu
dx(7) =subplot(4,2,8);
plot(result.freq, 20.*log10(abs(result.tube.out.T_Pu_flex)),...
    'LineWidth',2);
xlabel('frequency (Hz)');
ylabel('(dB)');
title('Power out / flow in');
hold all;

% link all x axes to enable synchronised zooming
linkaxes(dx(2:7),'x');

%% specify axis properties
for fignum = 1:length(fg)

h=get(fg(fignum),'children'); % gets a handle for all subplots

% instead of loop select individual axis and use gca
for l = 1:length(h)
    
    %   'Box'         , 'off'     , ...
    %   'YTick'       , 0:100:500, ...
    %     'XTick'     , 0:1:10    , ...
% 'FontName',   'Helvetica' , ...
% 'YTickLabel'  , [0; 2; 4; 6; 8; 10];% change the label

    set(h(l), ...
        'layer'       , 'top'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.01 .01] , ...
        'XMinorTick'  , 'on'      , ...
        'XGrid'       , 'on'      , ...
        'YGrid'       , 'off'      , ...
        'XColor'      , [.3 .3 .3], ...
        'YColor'      , [.3 .3 .3], ...
        'ZColor'      , [.3 .3 .3], ...
        'LineWidth'   , 1.5       , ...
        'FontSize',   12          );
    
    
    
end
set(h(7), 'XGrid', 'off');

%% specify overall figure properties
set(gcf, 'Color', 'w'); % sets the background colour
set(gcf,'units','inches')
set(gcf,'papersize',[8,12])
% set(gcf,'papertype', 'usletter')
set(gcf, 'PaperOrientation', 'portrait')

%% choose size
% set(gcf,'paperposition',[1,1,6.69,5]) % two columns
% set(gcf,'paperposition',[1,1,3.37,3]) % one column - need to adjust font sizes


 %% save examplename.eps at 300dpi with a tiff preview
 % this is cropped to the size of the figure
% print( gcf,'-depsc2','-tiff','-r300','-painters', strcat(name, '.eps') )
% savefig(strcat(name, '.fig'));


%% for surface and images use tiff
% print(gcf, '-dtiff', strcat(name, '.tiff'))

%% make a pdf with the figure on the page defined earlier
% print(gcf, '-dpdf', '-r300', '-painters', strcat(name, '.pdf'))


end