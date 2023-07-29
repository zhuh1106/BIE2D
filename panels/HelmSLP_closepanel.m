function A = HelmSLP_closepanel(omega,t,s,a,b,side, meth)
% special quadrature for 2d helmholtz
% via kernel split Laplace SLP special quadrature,

p = numel(s.x); % GL nodes order
if isfield('be',s), be = s.be; else, be = 1; end  % upsample factor of GL nodes
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