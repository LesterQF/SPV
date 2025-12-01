function [Jr_ini, s_tildeT, tildeT] = physical_enhancement(img, omega)

    imgvec = reshape(img, size(img, 1) * size(img, 2), 3);

    x_RGB(1, 1, 1:3) = mean(imgvec, 1); % unified spectrum
    x_mean = repmat(x_RGB, [size(img, 1) size(img, 2) 1]); % unified spectrum in each pixel

    scat_basis = x_mean ./ max(sqrt(sum(x_mean.^2, 3)), 0.001); % normalization
    fog_basis = img ./ max(sqrt(sum(img.^2, 3)), 0.001); % normalization
    cs_sim = repmat((((sum(scat_basis .* fog_basis, 3)))), [1 1 3]); % cos similarity
    scattering_light = (cs_sim) .* (sum(img, 3) ./ max(sum(x_mean, 3), 0.001)) .* x_mean;

    intial_img = img;
    [atmosphere, scattering_light] = get_atmosphere(intial_img, scattering_light);
    T = 1 - omega * scattering_light; % T = 1 - \tilde{t}
    T_ini = scattering_mask_sample(T);

    ShowR_d = (intial_img - atmosphere) ./ max(T_ini, 0.001) + atmosphere; % Winv.*
    mi = prctile2019(ShowR_d, 1, [1 2]);
    ma = prctile2019(ShowR_d, 99, [1 2]);
    Jr_ini = (ShowR_d - mi) ./ (ma - mi);

    Jr_ini = gamma0(uint8(Jr_ini * 255));

    tildeT = 1 - T;
    s_tildeT = 1 - T_ini;

    function [atmosphere, scatterlight] = get_atmosphere(image, scatterlight)
        for i = 1:1
            scatter_est = sum(scatterlight, 3);
            n_pixels = numel(scatter_est);
            n_search_pixels = floor(n_pixels * 0.001);
            image_vec = reshape(image, n_pixels, 3);
            [~, indices] = sort(scatter_est(:), 'descend');
            atmosphere = mean(image_vec(indices(1:n_search_pixels), :), 1);

            atmos(1, 1, :) = atmosphere;
            atmosphere = repmat(atmos, [size(scatter_est) 1]);

            sek = scatter_est(indices(n_search_pixels));

            scatterlight = scatterlight .* repmat(scatter_est <= sek, [1 1 3]) + ...
                (2/3 * sek - scatterlight) .* repmat(scatter_est > sek, [1 1 3]);
        end
    end
end

function img = gamma0(img)
    i = 0:255;
    f = ((i + 0.5) ./ 256).^(5/6);
    LUT(i + 1) = uint8(f * 256 - 0.5);
    img = rgb2ycbcr(img);
    img(:, :, 1) = LUT(img(:, :, 1) + 1);
    img = ycbcr2rgb(img);
end

function F = scattering_mask_sample(I)
    mSize = min(size(I, 1), size(I, 2));
    if mSize < 800
        r = 0.02;
    elseif mSize >= 800 && mSize < 1500
        r = 0.01;
    else
        r = 0.005;
    end
    I0 = imresize(I, r);
    F = imresize(I0, [size(I, 1), size(I, 2)], 'bicubic');
end

function y = prctile2019(varargin)
    par = inputParser();
    par.addRequired('x');
    par.addRequired('p');
    par.addOptional('dim', 1, @(x) isnumeric(x) || validateDimAll(x));
    par.addParameter('Delta', 1e3);
    par.addParameter('RandStream', []);

    par.parse(varargin{:});

    x = par.Results.x;
    p = par.Results.p;
    dim = par.Results.dim;
    delta = par.Results.Delta;
    rs = par.Results.RandStream;

    % Figure out which dimension prctile will work along.
    sz = size(x);

    % Permute the array so that the requested dimension is the first dim.
    if ~isequal(dim, 1)
        nDimsX = ndims(x);
        dim = sort(dim);
        perm = [dim setdiff(1:max(nDimsX, max(dim)), dim)];
        x = permute(x, perm);
    end
    sz = size(x);
    dimArgGiven = true;

    % Drop X's leading singleton dims, and combine its trailing dims.  This
    % leaves a matrix, and we can work along columns.
    work_dim = 1:numel(dim);

    work_dim = work_dim(work_dim <= numel(sz));
    nrows = prod(sz(work_dim));
    ncols = numel(x) ./ nrows;
    x = reshape(x, nrows, ncols);

    x = sort(x, 1);
    n = sum(~isnan(x), 1); % Number of non-NaN values in each column

    % For columns with no valid data, set n = 1 to get nan in the result
    n(n == 0) = 1;

    % If the number of non-nans in each column is the same, do all cols at once.
    if all(n == n(1))
        n = n(1);
        if isequal(p, 50) % make the median fast
            if rem(n, 2) % n is odd
                y = x((n + 1) / 2, :);
            else        % n is even
                y = (x(n / 2, :) + x(n / 2 + 1, :)) / 2;
            end
        else
            y = interpColsSame(x, p, n);
        end

    else
        % Get percentiles of the non-NaN values in each column.
        y = interpColsDiffer(x, p, n);
    end

    % Reshape Y to conform to X's original shape and size.
    szout = sz;
    szout(work_dim) = 1;
    szout(work_dim(1)) = numel(p);
    y = reshape(y, szout);

    % Undo the DIM permutation
    if dimArgGiven && ~isequal(dim, 1)
        y = ipermute(y, perm);
    end

    function y = interpColsSame(x, p, n)
        if isrow(p)
            p = p';
        end

        % Form the vector of index values (numel(p) x 1)
        r = (p / 100) * n;
        k = floor(r + 0.5); % K gives the index for the row just before r
        kp1 = k + 1;      % K+1 gives the index for the row just after r
        r = r - k;        % R is the ratio between the K and K+1 rows

        % Find indices that are out of the range 1 to n and cap them
        k(k < 1 | isnan(k)) = 1;
        kp1 = bsxfun(@min, kp1, n);

        % Use simple linear interpolation for the valid percentages
        y = (0.5 + r) .* x(kp1, :) + (0.5 - r) .* x(k, :);

        % Make sure that values we hit exactly are copied rather than interpolated
        exact = (r == -0.5);
        if any(exact)
            y(exact, :) = x(k(exact), :);
        end

        % Make sure that identical values are copied rather than interpolated
        same = (x(k, :) == x(kp1, :));
        if any(same(:))
            x = x(k, :); % expand x
            y(same) = x(same);
        end
    end

    function y = interpColsDiffer(x, p, n)
        % INTERPCOLSDIFFER A simple 1-D linear interpolation of columns that can
        % deal with columns with differing numbers of valid entries (vector n).

        [nrows, ncols] = size(x);

        % Make p a column vector. n is already a row vector with ncols columns.
        if isrow(p)
            p = p';
        end

        % Form the grid of index values (numel(p) x numel(n))
        r = (p / 100) * n;
        k = floor(r + 0.5); % K gives the index for the row just before r
        kp1 = k + 1;      % K+1 gives the index for the row just after r
        r = r - k;        % R is the ratio between the K and K+1 rows

        % Find indices that are out of the range 1 to n and cap them
        k(k < 1 | isnan(k)) = 1;
        kp1 = bsxfun(@min, kp1, n);

        offset = nrows * (0:ncols - 1);
        k = bsxfun(@plus, k, offset);
        kp1 = bsxfun(@plus, kp1, offset);

        y = (0.5 - r) .* x(k) + (0.5 + r) .* x(kp1);

        exact = (r == -0.5);
        if any(exact(:))
            y(exact) = x(k(exact));
        end

        same = (x(k) == x(kp1));
        if any(same(:))
            x = x(k); % expand x
            y(same) = x(same);
        end
    end

    function bool = validateDimAll(dim)
        bool = ((ischar(dim) && isrow(dim)) || ...
            (isstring(dim) && isscalar(dim) && (strlength(dim) > 0))) && ...
            strncmpi(dim, 'all', max(strlength(dim), 1));
    end
end