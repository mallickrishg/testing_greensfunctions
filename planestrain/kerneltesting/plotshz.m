function plotshz(shz,val,kmt)

% dx = shz.T(:)/2;
% dy = shz.W(:)/2;
% 
% x = [shz.x2(:)-dx,shz.x2(:)-dx,shz.x2(:)+dx,shz.x2(:)+dx,shz.x2(:)-dx]';
% y = [shz.x3(:)-dy,shz.x3(:)+dy,shz.x3(:)+dy,shz.x3(:)-dy,shz.x3(:)-dy]';
% 
% patch(x,y,val)

if exist('kmt','var')==1 && kmt==1
    scf = 1/1e3;
else 
    scf = 1;
end

if isfield(shz,'W')
    dx = shz.T(:)/2;
    dy = shz.W(:)/2;
    
    x = [shz.x2(:)-dx,shz.x2(:)-dx,shz.x2(:)+dx,shz.x2(:)+dx,shz.x2(:)-dx]';
    y = [shz.x3(:)-dy,shz.x3(:)+dy,shz.x3(:)+dy,shz.x3(:)-dy,shz.x3(:)-dy]';
    
    patch(x.*scf,y.*scf,val)

elseif isfield(shz,'tri')
    x = [shz.A(:,1),shz.B(:,1),shz.C(:,1)]';
    y = [shz.A(:,2),shz.B(:,2),shz.C(:,2)]';
    
    patch(x.*scf,y.*scf,val)
    
else
    disp('Nothing to plot')
end
end