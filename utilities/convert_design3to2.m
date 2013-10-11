function new_design = convert_design3to2(design, data, tr)
%CONVERT_DESIGN3TO2 Convert design from type 3 (cell arrays with onsets) to type 2 (arrays with marked onset volumes)
%
% INPUTS:
%
% OUTPUTS:
%
%  Yaroslav Halchenko                                            Dartmouth
%  web:     http://www.onerussian.com                              College
%  e-mail:  yoh@onerussian.com                              ICQ#: 60653192
%

% yoh: disclaimer -- yarik used matlab more than 10 years ago last time
new_design = {};
for chunk=1:length(design)
    dimtime = length(size(data{chunk}));
    d = zeros(size(data{chunk}, dimtime), size(design{chunk}, 1));
    for ev=1:length(design{chunk})
        onsets = design{chunk}{ev};
        for i=1:length(onsets)
            vol = round(onsets(i)/tr);
            d(vol+1, ev) = 1;
        end
    end
    new_design{chunk} = d;
end
assert(length(new_design) == length(design));
