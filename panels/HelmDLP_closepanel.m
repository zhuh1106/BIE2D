function A = HelmDLP_closepanel(omega,t,s,a,b,side, meth)
% special quadrature for 2d helmholtz
% via kernel split Laplace DLP special quadrature

if nargin == 0, test_HelmDLP_closepanel(); return; end

p = numel(s.x); % GL nodes order
if isfield(s,'be'), be = s.be; else, be = 1; end  % upsample factor of GL nodes
if be > 1
  [x, w] = gauss(p); w = w';
  [xf, wf] = gauss(ceil(be*p)); wf = wf';
  Imn = interpmat(x, xf);
  hl = mean(s.w./(w.*s.sp));
  s.w = wf; % overwrite source panel
  s.x = Imn*s.x; s.xp = Imn*s.xp;
  s.sp = abs(s.xp); s.tang = s.xp./s.sp; s.nx = -1i*s.tang;    % outward unit normals
  s.wxp = hl*s.w.*s.xp; % complex speed weights (Helsing's wzp)
  s.w = abs(s.wxp); % speed weights
end

% Variants of an explicit kernel-split panel-based Nystrom discretization scheme for Helmholtz boundary value problems.
% Johan Helsing & Anders Holst. p 696, eq(26), pay attention to possible typo in laplace DLP coefficient 1/(2*pi) instead of 1/pi
lapslpmat = -2*pi*LapSLP_closepanel(t,s,a,b,side, meth); % close eval for log kernel (Laplace SLP)
lapdlpmat = -2*pi*LapDLP_closepanel(t,s,a,b,side, meth); % close eval for Cauchy kernel (Laplace DLP)
u = 1i/4*(omega*abs(s.x(:).'-t.x(:))).*besselh(1,omega*abs(s.x(:).'-t.x(:))).*imag(s.wxp(:).'./(t.x(:)-s.x(:).'));
A = ( u + 2/pi*log(abs(s.x(:).'-t.x(:))).*imag(u) + 1/(2*pi)*real(s.nx(:).'./(s.x(:).'-t.x(:))).*s.w' ) ...
         - 2/pi*(lapslpmat./s.w').*imag(u) - 1/(2*pi)*real(lapdlpmat); 
if be > 1, A = A*Imn; end

end

function test_HelmDLP_closepanel
% test a bvp

k = 10;                       % wave number           
a = .3; w = 5;                % smooth wobbly radial shape params
verb = 1;                     % 0 text only; 1 makes a figure
s = wobblycurve(1,a,w,100);   % parametrix descrip of curve (dummy # pts)
Np = 16*6;                      % # panels (not a convergence study, just a test)
p = 16;                       % panel order
[pa tpan s] = quadr_uniform_panels(s,Np,p);
zpan = s.Z(tpan);             % panel endpoint locations, in C plane

%%% self mat
Smat = 1i/4*besselh(0,k*abs(s.x(:).'-s.x(:))).*s.w(:)'; % naive eval mat
Dmat = 1i/4*(k*abs(s.x(:).'-s.x(:))).*besselh(1,k*abs(s.x(:).'-s.x(:))).*imag(s.wxp(:).'./(s.x(:)-s.x(:).'));
for i=1:Np % close eval
  pa{i}.be = 2; % upsample factor
  cflag=(min(abs(s.x((i-1)*p+1:i*p).'-s.x(:)),[],2)/sum(s.w((i-1)*p+1:i*p))<=0.7); % close targ rule
  Smat(cflag,(i-1)*p+1:i*p) = HelmSLP_closepanel(k,struct('x',s.x(cflag)),pa{i},zpan(i),zpan(i+1),'e','h');
  Dmat(cflag,(i-1)*p+1:i*p) = HelmDLP_closepanel(k,struct('x',s.x(cflag)),pa{i},zpan(i),zpan(i+1),'e','h');
end

%%% solve
z0=0.3+0.1i; % interior source
rhs=besselh(0,k*abs(z0-s.x));
rho = (Smat+Dmat)\rhs; % Helmholtz SLP

%%% eval mat
x0 = 1.2; nx = 100; gx = max(abs(s.x))*linspace(-x0,x0,nx); % grid for soln errors...
[xx yy] = meshgrid(gx); zz = xx+1i*yy; ii = ~s.inside(zz);  % outside
g.x = zz(ii(:));
Stmat = 1i/4*besselh(0,k*abs(s.x(:).'-g.x(:))).*s.w(:)'; % naive eval mat
Dtmat = 1i/4*(k*abs(s.x(:).'-g.x(:))).*besselh(1,k*abs(s.x(:).'-g.x(:))).*imag(s.wxp(:).'./(g.x(:)-s.x(:).'));
for i=1:Np % close eval
  pa{i}.be = 2; % upsample factor
  cflag=(min(abs(s.x((i-1)*p+1:i*p).'-g.x(:)),[],2)/sum(s.w((i-1)*p+1:i*p))<=0.7); % close targ rule
  Stmat(cflag,(i-1)*p+1:i*p) = HelmSLP_closepanel(k,struct('x',g.x(cflag)),pa{i},zpan(i),zpan(i+1),'e','h');
  Dtmat(cflag,(i-1)*p+1:i*p) = HelmDLP_closepanel(k,struct('x',g.x(cflag)),pa{i},zpan(i),zpan(i+1),'e','h');
end
Unum = (Stmat+Dtmat)*rho;
Uexa = besselh(0,k*abs(z0-g.x));

%%% error
err = nan*(1+1i)*zz;
err(ii(:)) = abs(Unum - Uexa);

figure(1),clf,imagesc(log10(abs(err)));colorbar
axis equal
fprintf('S+D max err = %.3g\n',max(err(:)))
end