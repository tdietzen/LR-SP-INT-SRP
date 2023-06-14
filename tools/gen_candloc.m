function [loc, Delta_t_p] = gen_candloc(micPos, c, varargin) 
% [loc, Delta_t_p] = gen_candloc(micPos c, varargin)
% converts grid points to candidate locations and computes TDOAs.
%
% IN:
% micPos         microphone positions - channels x coordinates
% c              speed of sound
% varargin       arrays of grid points; if 3 arrays are provided, the grid
%                points are interpreted as Cartesian (x, y, z), and if 2 
%                arrays are provided, the grid points are interpeted as 
%                angles (polar, azimuth)
%
% OUT:
% loc            location coordinates or DOA vectors - candidate locations 
%                x microphone pairs
% Delta_t_p      TDOAs - candidate locations x microphone pairs


% interpret grid points as NF or FF
switch length(varargin)
    
    case 2
        mode = 'FF';
        pol = varargin{1};
        azi = varargin{2};
        
    case 3
        mode = 'NF';
        x = varargin{1};
        y = varargin{2};
        z = varargin{3};
                
    otherwise
        error('search grid generation requires two (FF) or three (NF) dimension vectors');
        
end

% compute locations and TDOAs
switch mode
    
    case 'FF'
        
        i = 0;
        
        for idx_pol = 1:length(pol)
            for idx_azi = 1:length(azi)
                
                i = i + 1;
                % location in far field spherical coordinates
                [Delta_t_p(i,:), loc(i,:)] = spher2TDOA(micPos, pol(idx_pol),  azi(idx_azi), c);
                if pol(idx_pol) == 0 || pol(idx_pol) == 180
                    break  % we do not need to loop further over azi because DOAs won't change
                end
                
            end
        end
        
    case 'NF'
        
      i = 0;
        
      for idx_z = 1:length(z)
          for idx_y = 1:length(y)
              for idx_x = 1:length(x)
                  
                  i = i + 1;
                  % location in cartesian coordinates
                  loc(i, :) = [x(idx_x) y(idx_y) z(idx_z)];
                  % TDOA
                  Delta_t_p(i,:) = cart2TDOA(micPos, [x(idx_x) y(idx_y) z(idx_z)], c);
                  
                  
              end
          end
      end
end 

end

function Delta_t = cart2TDOA(micPos, loc, c)
% Cartesian coordinates to TDOAs

M = size(micPos,1);
P = M*(M-1)/2;

% propagation time
dist = sqrt(sum((micPos - repmat(loc, size(micPos, 1), 1)).^2, 2));
propTime = dist/c;

% TDOA
Delta_t = zeros(P,1);
p = 0;
for mprime = 1:M
    for m = mprime+1:M
        p = p+1;
        Delta_t(p) = propTime(m) - propTime(mprime);
    end
end

end


function [Delta_t_p, incident_vec] = spher2TDOA(micPos, pol, azi, c)
% spherical coordinates to TDOAs

M = size(micPos,1);
P = M*(M-1)/2;

% incident vector
incident_vec = -[sin(deg2rad(pol))*cos(deg2rad(azi)),...
    sin(deg2rad(pol))*sin(deg2rad(azi)),...
    cos(deg2rad(pol))];

b = transp(incident_vec);

% TDOAs
Delta_t_p = zeros(P,1);
p = 0;
for mprime = 1:M
    for m = mprime+1:M
        p = p+1;
  
        a = micPos(m,:).' - micPos(mprime,:).';
        
        Delta_t_p(p) = (norm(a)/c)*(a.'*b/(norm(a)*norm(b)));

    end
end

end
