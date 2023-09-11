function rcv = construct_box(x0,y0,L_x,L_y,dx)
% INPUTS
% x0,y0 - starting point to construct mesh

rcv = [];

xval = linspace(x0,x0 + L_x,round(L_x/dx));
yval = linspace(y0,y0 + L_y,round(L_y/dx));
nxv = length(xval);
nyv = length(yval);

x1 = [xval(1:end-1),...
      xval(end)*ones(1,nyv-1),...
      xval(end:-1:2),...
      xval(1)*ones(1,nyv-1)];
x2 = [xval(2:end),...
      xval(end)*ones(1,nyv-1),...
      xval(end-1:-1:1),...
      xval(1)*ones(1,nyv-1)];

y1 = [yval(1)*ones(1,nxv-1),...
      yval(1:end-1),...
      yval(end)*ones(1,nxv-1),...
      yval(end:-1:2)];
y2 = [yval(1)*ones(1,nxv-1),...
      yval(2:end),...
      yval(end)*ones(1,nxv-1),...
      yval(end-1:-1:1)];

rcv.x1 = x1';
rcv.x2 = x2';
rcv.y1 = y1';
rcv.y2 = y2';

rcv.xc = [(rcv.x1+rcv.x2)./2 , (rcv.y1+rcv.y2)./2];
rcv.W = sqrt((rcv.y2-rcv.y1).^2 + (rcv.x2-rcv.x1).^2);
rcv.N = length(rcv.W);

rcv.dip = atan2d(rcv.y1-rcv.y2,rcv.x1-rcv.x2);
rcv.nv = [-(rcv.y2-rcv.y1)./rcv.W,(rcv.x2-rcv.x1)./rcv.W];

end