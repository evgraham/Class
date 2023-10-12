syms z1 p1;
Dz1 = EG_C2D_matched(RR_tf(z1, 10, 1), .01)
function [Dz] = EG_C2D_matched(Ds, h, omega_bar, causality)
    syms z1 p1;
    % function [Dz]=EG_C2D_matched(Ds,h,omega_bar,causality)
    % Convert a continuous D(s) (RR_tf class) to a discrete D(z) using 
    % the matched z-transform method. If an omega_bar is specififed as an 
    % input it will be used to caluclate the gain, otherwise the default 
    % value of 0 is used. The causality input determines whether a 
    % semi-causal or strictly causal D(z) is determined. Please input 
    % 'semi' or 'strict' accordingly or a default value of 'semi' will be 
    % used. 
    % Determine user inputs and switch to defaults if necessary.
    %syms z1 p1;
    if nargin == 2
        omega_bar=0;
        causality='semi';
    elseif nargin < 4 || isempty(omega_bar)
        omega_bar=0;
    elseif nargin < 4 || isempty(causality)
        causality='semi';
    end
    % disp(omega_bar); disp(causality);
    % Determine gain @ s=i*(omega_bar)
    gain = RR_evaluate(Ds, omega_bar*1i);
    % Initialize poles and zeros
    zero = Ds.z; pole = Ds.p;
    % Determine tf order
    m = length(Ds.num.n); n = length(Ds.den.n);
    % Find infinte zeros
    infinite_zeros = m-n;
    % Initialize discrete zeros
    zs=zeros(1, length(zero)+infinite_zeros);
    % Add discrete zeros for infinte zeros
    for k=length(zero):length(zero)+infinite_zeros
        zs(k)=-1;
    end
    % When a strict tf is needed, add an zero at infinity
    if strcmp(causality, 'strict')
        zs(end) = inf;
    end
    
    % Match continuous zeros to discrete zeros
    for j=1:length(zero)
        zs(j) = exp(subs(zero(j))*h);
    end
    
    % Match continuous poles to discrete poles
    ps=[];
    for j=1:length(pole)
        % match pole according to z=e^(sh)        
        ps(j) = exp(eval(subs(pole(j)))*h);
    end
    Dz = RR_tf(zs, ps, gain); Dz.h=h;
end % function EG_C2D_matched
