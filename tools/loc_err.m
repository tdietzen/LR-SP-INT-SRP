function [e_loc, src_loc_hat] = loc_err(z, cand_loc, src_loc, scenario)
% [e_loc, src_loc_hat] = loc_err(z, cand_loc, src_loc, scenario) computes
% the localization error.
%
% IN:
% z            SRP map
% cand_loc     candidate locations
% src_loc      true source locations
% scenario     scenario ('NF' or 'FF')
%
% OUT:
% e_loc        localization error
% src_loc_hat  estimated source location

[~, max_idx] = max(z); % find maximum

q_s = src_loc.';
q_s_hat = cand_loc(max_idx,:).';  % peak picking

% localization error
switch scenario
    case 'NF' 
        e_loc = norm(q_s_hat - q_s);  
    case 'FF'   
        e_loc = rad2deg(acos(q_s_hat.'*q_s));
end

src_loc_hat = q_s_hat.'; % estimated source location

end
