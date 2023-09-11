function [Disp,Stress] = LTkernelFS(x,y,xe,ye,a,Beta,Fx,Fy,nu,mu)

f = constkernel(x - xe, y - ye, a, nu);
[ux, uy, sxx, syy, sxy] = trac2dispstress(Fx, Fy, f, y - ye, mu, nu);

Disp = [ux(:),uy(:)];
Stress = [sxx(:),sxy(:),syy(:)];

end

function f = constkernel(x, y, a, nu)
    f = zeros(length(x), 7);
    leadingconst = 1/(4*pi*(1-nu));

    for i=1:length(x)
        % f(x) = f_1
        f(i, 1) = -leadingconst * (y(i) * (atan2(y(i), (x(i)-a)) - atan2(y(i), (x(i)+a))) - (x(i) - a) * log(sqrt((x(i) - a)^2 + y(i)^2)) + (x(i)+a) * log(sqrt((x(i)+a)^2 + y(i)^2)));

        % df/dy = f_2
        f(i, 2) = -leadingconst * ((atan2(y(i), (x(i)-a)) - atan2(y(i), (x(i)+a))));

        % df/dx = f_3
        f(i, 3) = leadingconst * (log(sqrt((x(i)-a)^2 + y(i)^2)) - log(sqrt((x(i)+a)^2 + y(i)^2)));

        % d2f/dxdy = f_4
        f(i, 4) = leadingconst * (y(i) / ((x(i)-a)^2 + y(i)^2) - y(i) / ((x(i)+a)^2 + y(i)^2));

        % d2f/dxdx = -d2f/dydy = f_5
        f(i, 5) = leadingconst * ((x(i)-a) / ((x(i)-a)^2 + y(i)^2) - (x(i)+a) / ((x(i)+a)^2 + y(i)^2));

        % d3f/dxdydy = -d3f/dxdxdx = f_6
        f(i, 6) = leadingconst * (((x(i)-a)^2 - y(i)^2) / ((x(i)-a)^2 + y(i)^2)^2 - ((x(i)+a)^2 - y(i)^2) / ((x(i)+a)^2 + y(i)^2)^2);

        % d3f/dydydy = -d3f/dxdxdy = f7
        f(i, 7) = 2 * y(i) / (4 * pi * (1 - nu)) * ((x(i)-a) / ((x(i)-a)^2 + y(i)^2)^2 - (x(i)+a) / ((x(i)+a)^2 + y(i)^2)^2);
    end
end


function [ux, uy, sxx, syy, sxy] = trac2dispstress(xcomp, ycomp, f, y, mu, nu)
    ux = zeros(size(y));
    uy = zeros(size(y));
    sxx = zeros(size(y));
    syy = zeros(size(y));
    sxy = zeros(size(y));
    for i=1:length(y)
        ux(i) = xcomp / (2.0 * mu) * ((3.0 - 4.0 * nu) * f(i, 1) + y(i) * f(i, 2)) + ycomp / (2.0 * mu) * (-y(i) * f(i, 3));
        uy(i) = xcomp / (2.0 * mu) * (-y(i) * f(i, 3)) + ycomp / (2.0 * mu) * ((3.0 - 4.0 * nu) * f(i, 1) - y(i) * f(i, 2));
        sxx(i) = xcomp * ((3.0 - 2.0 * nu) * f(i, 3) + y(i) * f(i, 4)) + ycomp * (2.0 * nu * f(i, 2) + y(i) * -f(i, 5));
        syy(i) = xcomp * (-1.0 * (1.0 - 2.0 * nu) * f(i, 3) - y(i) * f(i, 4)) + ycomp * (2.0 * (1.0 - nu) * f(i, 2) - y(i) * -f(i, 5));
        sxy(i) = xcomp * (2.0 * (1.0 - nu) * f(i, 2) + y(i) * -f(i, 5)) + ycomp * ((1.0 - 2.0 * nu) * f(i, 3) - y(i) * f(i, 4));
    end
end