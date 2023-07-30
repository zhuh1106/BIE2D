function A = HelmSLP_closepanel(omega,t,s,a,b,side, meth)
% special quadrature for 2d helmholtz
% via kernel split Laplace SLP special quadrature,

if nargin == 0, test_HelmSLP_closepanel(); return; end

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
% Johan Helsing & Anders Holst. p 696, eq(20)
lapslpmat = -2*pi*LapSLP_closepanel(t,s,a,b,side, meth); % close eval for log kernel (Laplace SLP)
u = 1i/4*besselh(0,omega*abs(s.x(:).'-t.x(:))).*s.w'; % naive eval
A = ( u + 2/pi*log(abs(s.x(:).'-t.x(:))).*imag(u) ) ... % S_0k(r,r') = S_k(r,r')+2/pi*log|r-r'|*Im{S_k{r,r'}} smooth part
         - 2/pi*(lapslpmat./s.w').*imag(u); % -2/pi*log|r-r'|*Im{S_k(r,r')}
if be > 1, A = A*Imn; end

end

function test_HelmSLP_closepanel
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
for i=1:Np % close eval
  pa{i}.be = 2; % upsample factor
  cflag=(min(abs(s.x((i-1)*p+1:i*p).'-s.x(:)),[],2)/sum(s.w((i-1)*p+1:i*p))<=0.7); % close targ rule
  Smat(cflag,(i-1)*p+1:i*p) = HelmSLP_closepanel(k,struct('x',s.x(cflag)),pa{i},zpan(i),zpan(i+1),'e','h');
end

%%% solve
z0=0.3+0.1i; % interior source
rhs=besselh(0,k*abs(z0-s.x));
rho = Smat\rhs; % Helmholtz SLP

%%% eval mat
x0 = 1.2; nx = 100; gx = max(abs(s.x))*linspace(-x0,x0,nx); % grid for soln errors...
[xx yy] = meshgrid(gx); zz = xx+1i*yy; ii = ~s.inside(zz);  % outside
g.x = zz(ii(:));
Stmat = 1i/4*besselh(0,k*abs(s.x(:).'-g.x(:))).*s.w(:)'; % naive eval mat
for i=1:Np % close eval
  pa{i}.be = 2; % upsample factor
  cflag=(min(abs(s.x((i-1)*p+1:i*p).'-g.x(:)),[],2)/sum(s.w((i-1)*p+1:i*p))<=0.7); % close targ rule
  Stmat(cflag,(i-1)*p+1:i*p) = HelmSLP_closepanel(k,struct('x',g.x(cflag)),pa{i},zpan(i),zpan(i+1),'e','h');
end
Unum = Stmat*rho;
Uexa = besselh(0,k*abs(z0-g.x));

%%% error
err = nan*(1+1i)*zz;
err(ii(:)) = abs(Unum - Uexa);

figure(1),clf,imagesc(log10(abs(err)));colorbar
axis equal
fprintf('S max err = %.3g\n',max(err(:)))
end