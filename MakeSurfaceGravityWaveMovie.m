
file = '/Users/jearly/Desktop/BoussinesqWave.nc';
FramesFolder ='/Users/jearly/Desktop/BoussinesqWaveFrames';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Make the frames folder
%
if exist(FramesFolder,'dir') == 0
	mkdir(FramesFolder);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Read in the problem dimensions
%
x = ncread(file, 'x');
t = ncread(file, 'time');

ssh = ncread(file, 'SSH');
bathymetry = ncread(file, 'bathymetry');
sponge = ncread(file, 'sponge');
wave_generator = ncread(file, 'wave_generator');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Setup the figure
%
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 1280 720])

for iTime=1:length(t)
    clf
	% 'area' fills everything between the line and the zero axis,
	harea = area(x/1000, -bathymetry);
	% so we need to color that area 'blue'
	set(harea(1),'FaceColor',[0.46 0.58 0.96], 'EdgeColor', [0.46 0.58 0.96],  'LineWidth', 1); % color the water blue
	% but then color the plot background brown, to represent the bathymetry.
	set(gca,'Color',[.60 .40 .20]);
	hold on
	% We instructed matlab not to draw a line at the bathymetry, so let's do that now.
	plot(x/1000, -bathymetry, 'LineWidth', 2, 'Color', 'black')
    
	wave = ssh(:,iTime);
    
	% adjust the axis to show a bit more earth and wave surface.
    xlim([min(x/1000) max(x/1000)])
	ylim([-max(bathymetry)*1.1 max(max(ssh))*2.0])
    set( gca, 'FontSize', 18);
	xlabel('distance (km)','FontName', 'Helvetica','FontSize',24)
	ylabel('depth (m)','FontName', 'Helvetica','FontSize',24)
	title('Surface gravity waves','FontName', 'Helvetica','FontSize',36);
    text(0,double(-max(bathymetry)),sprintf('t=%d seconds',round(t(iTime))),'FontName', 'Helvetica','FontSize',20, 'BackgroundColor', 'white')


	
		% Find the regions where the wave is expected to break.
	[mins maxes] = FindCriticalSlopeIndices(x, wave);

	% now we plot the surface waves
	% Because harea fills between the line and the axis at zero, we have to fill the upper portion and the lower portion separately.
	waveTop = wave;
	waveTop(find( waveTop < 0)) = 0;
	waveBottom = wave;
	waveBottom(find( waveTop > 0)) = 0;

	harea = area(x/1000, 50*ones(size(x)));
	set(harea(1),'FaceColor',[1.0 1.0 1.0], 'EdgeColor', [1.0 1.0 1.0],  'LineWidth', .1); % color the surface white

	harea = area(x/1000, wave);
	set(harea(1),'FaceColor',[1.0 1.0 1.0], 'EdgeColor', [1.0 1.0 1.0],  'LineWidth', .001); % color the surface white
	
	if 1
		for i=1:length(mins)
			range = mins(i):maxes(i);
			harea = area(x(range)/1000, 50*ones(size(x(range))) );
			set(harea(1),'FaceColor',[1.0 0.0 0.0], 'EdgeColor', [1.0 1.0 1.0],  'LineWidth', .1); % color the surface white
			
			harea = area(x(range)/1000, wave(range));
			set(harea(1),'FaceColor',[1.0 1.0 1.0], 'EdgeColor', [1.0 1.0 1.0],  'LineWidth', .001); % color the surface white
			
		end
	end
	
	harea = area(x/1000, waveTop);
	set(harea(1),'FaceColor',[0.46 0.58 0.96], 'EdgeColor', 'none'); % color the water blue

	if 0
		for i=1:length(mins)
			range = mins(i):maxes(i);
			harea = area(x(range)/1000, waveTop(range));
			set(harea(1),'FaceColor',[1.0 0.0 0.0], 'EdgeColor', 'none'); % color the water red
		end
	end
	
	plot(x/1000, wave, 'LineWidth', 1.0, 'Color', 'black')
    
	% tighten the plot
	T = get(gca,'tightinset');
	set(gca,'position',[T(1) T(2) 1-T(1)-T(3) 1-T(2)-T(4)-.05]);
    
    output = sprintf('%s/t_%03d.png', FramesFolder,iTime-1);
    export_fig(output,'-r72')
end

