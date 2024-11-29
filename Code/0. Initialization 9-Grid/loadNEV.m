function NEV = loadNEV(fn, input)

    if ~exist('input', 'var')
        NEV = openNEV(fn, 'nomat', 'nosave');

    elseif strcmpi(input, 'noread')
        NEV = openNEV(fn, 'nomat', 'nosave');

    elseif ~strcmpi(input, 'read')
        disp('Invalid input.');

    else
        NEV = openNEV(fn, 'nomat', 'nosave');
    end
end