% Problem 3 of HW1
% Code to use the built-in MATLAB c2d function
% Dz_builtin = c2d(tf([1 1], [1 10 0]), 0.01, 'matched');
sympref('FloatingPointOutput', true);
% TEST CASES BELOW
% Dz1 = EG_C2D_matched(RR_tf([1 1], [1 10 0]), .01, [], []);
% syms z1 p1;
% Dz1 = EG_C2D_matched(RR_tf([1 z1], [1 p1 0]), .01, .1, 'strict');
% FUNCTION DECLARATION
function [Dz] = EG_C2D_matched(Ds, h, omega_bar, causality)
    % function [Dz]=EG_C2D_matched(Ds,h,omega_bar,causality)
    % Convert a continuous D(s) (RR_tf class) to a discrete D(z) using 
    % the matched z-transform method. If an omega_bar is specififed as an 
    % input it will be used to caluclate the gain, otherwise the default 
    % value of 0 is used. The causality input determines whether a 
    % semi-causal or strictly causal D(z) is determined. Please input 
    % 'semi' or 'strict' accordingly or a default value of 'semi' will be 
    % used. IF LEAVING OPTION BLANK, USE []
    % Determine user inputs and switch to defaults if necessary.
    if nargin == 2, omega_bar=0; causality='semi';
    elseif nargin == 3 && isempty(causality), causality='semi';
    elseif nargin == 3 && isempty(omega_bar), omega_bar=0;
    elseif nargin == 4 && isempty(causality) && ~isempty(omega_bar) 
        causality='semi';
    elseif nargin == 4 && isempty(omega_bar) && ~isempty(causality) 
        omega_bar=0;
    elseif nargin == 4 && isempty(omega_bar) && isempty(causality), 
        omega_bar=0; causality='semi';
    end
    % Initialize poles and zeros
    zero = Ds.z; pole = Ds.p;
    % Test if pole at 0 and omega_bar==0
    for p=1:length(pole)
        if pole(p)==0 && omega_bar==0
            error(['If a pole at zero is present, omega_bar cannot be ' ...
                'equal to 0. Please input a different omega_bar.'])
        end
    end
    % Determine tf order
    m = Ds.num.n; n = Ds.den.n;
    % Find infinte zeros
    infinite_zeros = n-m;
    % Initialize discrete zeros
    if infinite_zeros > 0
        zs=sym(zeros(1, length(zero)+infinite_zeros));
        % Add discrete zeros at -1 for infinite zeros
        for k=length(zero)+1:length(zero)+infinite_zeros
            zs(k)=-1;
        end
        % When a strict tf is needed, add a zero at infinity
        if strcmp(causality, 'strict')
            zs(end) = [];
        end
    else
        zs=sym(zeros(1, length(zero)));
    end
    % Match continuous zeros to discrete zeros
    for j=1:length(zero), zs(j) = exp(subs(zero(j))*h);
    end
    % Match continuous poles to discrete poles according to z=e^(sh)
    ps=sym([]);
    for j=1:length(pole), ps(j) = exp(subs(pole(j))*h);
    end
    % Determine correct gain value
    Dz = RR_tf(zs, ps, 1); Dz.h=h;
    % Determine gain @ s=i*(omega_bar)
    sgain = RR_evaluate(Ds, omega_bar*1i);
    zgain = RR_evaluate(Dz, exp(omega_bar*1i*h));
    if sgain==zgain, gain = 1;
    else, gain=sgain/zgain;
    if isnan(gain), gain=1;
    end
    end
    % Final gain of Dz
    Dz = RR_tf(zs, ps, gain);
end % function EG_C2D_matched
