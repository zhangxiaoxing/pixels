%% Point attractor
% Define the mesh surface
[X, Y] = meshgrid(-2.3:0.1:2.3, -1.3:0.1:1.3);
zr=0.5*sqrt((abs(X)-1).^2+(Y).^2);
gk = fspecial('gaussian', [7,7], 1);
Z=conv2(zr,gk,'valid');
Z=Z+rand(size(Z))*0.1;

% Plot the mesh surface
fh=figure('Position',[100,100,640,320]);
sh=surf(X(4:end-3,4:end-3), Y(4:end-3,4:end-3), Z);
sh.EdgeAlpha=0.15;
zlim([0,1])
view(-15,60)
colormap('turbo')
grid OFF


%% synfire chain

[X, Y] = meshgrid(-2.3:0.1:2.3, -1.3:0.1:1.3);
zr=zeros(size(X));
zr(:,[12:16,32:36])=-1;
zr=zr+repmat(linspace(-1,1,27).',1,47);

gk = fspecial('gaussian', [7,7], 1);
Z=conv2(zr,gk,'valid');
Z=Z+rand(size(Z))*0.2;

fh=figure('Position',[100,100,640,320]);
sh=surf(X(4:end-3,4:end-3), Y(4:end-3,4:end-3), Z);
sh.EdgeAlpha=0.15;
zlim([-2,1])
view(-15,60)
colormap('turbo')
grid OFF


%% Cell assembly
zr=(double(imread(fullfile('binary','cellasm_illustration_S.png')))-255)./255;
gk = fspecial('gaussian', [5,5],2);
Z=conv2(zr,gk,'same');
Z=Z+rand(size(Z))*0.1;

fh=figure('Position',[100,100,640,320]);
sh=surf(Z);
sh.EdgeAlpha=0.15;
zlim([-2,1])
view(-15,60)
colormap('turbo')
grid OFF




%% Trajectory
zr=(double(imread(fullfile('binary','traj.png')))-255)./240;
gk = fspecial('gaussian', [5,5],2);
Z=flip(conv2(zr,gk,'same'))+rand(size(Z))*0.1;
Z=Z+rand(size(Z))*0.1;

fh=figure('Position',[100,100,640,320]);
sh=surf(Z);
sh.EdgeAlpha=0.15;
zlim([-2,1])
view(-15,60)
colormap('turbo')
grid OFF
