function rcv = construct_triangle(A,B,C,dx,Ntopo)

pts = [A;B;C;A];

rcv = [];

% construct triangle
r1 = sqrt((A(1)-B(1)).^2 + (A(2)-B(2)).^2);
r2 = sqrt((B(1)-C(1)).^2 + (B(2)-C(2)).^2);
r3 = sqrt((C(1)-A(1)).^2 + (C(2)-A(2)).^2);

disp(['Length of AB = ' num2str(r1) ', dip = ' num2str(atan2d(A(2)-B(2),A(1)-B(2)))])
disp(['Length of BC = ' num2str(r2) ', dip = ' num2str(atan2d(B(2)-C(2),B(1)-C(2)))])
disp(['Length of CA = ' num2str(r3) ', dip = ' num2str(atan2d(C(2)-A(2),C(1)-A(2)))])

Nvec = [round(r1/dx);round(r2/dx);round(r3/dx)];
N = sum(Nvec);

x1 = zeros(N,1);
x2 = zeros(N,1);
y1 = zeros(N,1);
y2 = zeros(N,1);

for i = 1:3
    xval = linspace(pts(i,1),pts(i+1,1),Nvec(i)+1);
    yval = linspace(pts(i,2),pts(i+1,2),Nvec(i)+1);
    if i == 1
        st = 1;
        en = Nvec(i);
    else
        st = sum(Nvec(1:(i-1)))+1;
        en = sum(Nvec(1:i));
    end

    x1(st:en) = xval(1:end-1);
    x2(st:en) = xval(2:end);
    y1(st:en) = yval(1:end-1);
    y2(st:en) = yval(2:end);
    
end
% construct free surface
topox = linspace(-50,50,Ntopo+1)';
topoy = 0.*topox;

rcv.x1 = [x1;topox(1:end-1)];
rcv.x2 = [x2;topox(2:end)];
rcv.y1 = [y1;topoy(1:end-1)];
rcv.y2 = [y2;topoy(2:end)];


rcv.xc = [(rcv.x1+rcv.x2)./2 , (rcv.y1+rcv.y2)./2];
rcv.W = sqrt((rcv.y2-rcv.y1).^2 + (rcv.x2-rcv.x1).^2);
rcv.N = length(rcv.W);

rcv.dip = atan2d(rcv.y1-rcv.y2,rcv.x1-rcv.x2);
rcv.nv = [-(rcv.y2-rcv.y1)./rcv.W,(rcv.x2-rcv.x1)./rcv.W];


end