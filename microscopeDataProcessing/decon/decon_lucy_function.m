function [J_2, err_mat, k] = decon_lucy_function(I, PSF, NUMIT, varargin)
% adapted from matlab deconvlucy.m
% 
% xruan (05/18/2021): add support for early stop with stop criteria
% xruan (05/19/2021): add support for gpu computing
% xruan (06/02/2021): add support for fix iteration decon
% xruan (06/03/2021): add support for debug (evaluate error and save
% intermediate result).
% xruan (07/09/2021): add output for the number of iterations run; also
% output err mat for fixed iterations.
% xruan (11/11/2021): when using GPU, use single precision unless the data
% is in double
% xruan (11/12/2021): refactor code to save memory (J_4), and disable
% weighting or subsampling
% xruan (11/17/2022): add save16bit, bbox, background option for faster processing 
% on GPU, also remove err_thrsh and fixIter parameters
% xruan (08/31/2023): add scaleFactor, deconOffset, damper, and
%   EdgeErosion, and use input parser for parameters.


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('I', @isnumeric);
ip.addRequired('PSF', @isnumeric);
ip.addRequired('NUMIT', @isnumeric);
ip.addParameter('Background', [], @isnumeric);
ip.addParameter('useGPU', true, @islogical); % use GPU processing
ip.addParameter('save16bit', true, @islogical);
ip.addParameter('dampFactor', 1, @isnumeric); % damp factor for decon result
ip.addParameter('scaleFactor', 1, @isnumeric); % scale factor for result
ip.addParameter('deconOffset', 0, @isnumeric); % offset for decon result
ip.addParameter('EdgeErosion', 0, @isnumeric); % edge erosion for decon result
ip.addParameter('bbox', [], @isnumeric); % bounding box to crop data after decon
ip.addParameter('debug', false, @islogical);
ip.addParameter('debug_folder', './debug/', @ischar);
ip.addParameter('saveStep', 1, @isnumeric); % save intermediate results every given iterations

ip.parse(I, PSF, NUMIT, varargin{:});

pr = ip.Results;
Background = pr.Background;
useGPU = pr.useGPU;
save16bit = pr.save16bit;
dampFactor = pr.dampFactor;
scaleFactor = pr.scaleFactor;
deconOffset = pr.deconOffset;
EdgeErosion = pr.EdgeErosion;
bbox = pr.bbox;
debug = pr.debug;
debug_folder = pr.debug_folder;
saveStep = pr.saveStep;

if isempty(Background)
    Background = 0;
end

PSF = single(PSF);
sizeI = size(I);

% 1. Prepare PSFs
persistent H psz sq_p gpuNum last_useGPU
if isempty(psz) || ~(numel(psz) == numel(sizeI) && all(psz == sizeI) && last_useGPU == useGPU ...
        && sq_p == sum(PSF .^ 2, 'all')) || (useGPU && ~isempty(H) && gpuNum > 0 && ~existsOnGPU(H))
    psz = sizeI;
    
    gpuNum = 0;
    last_useGPU = useGPU;
    if useGPU
        gpuNum = gpuDeviceCount("available");
    end

    sq_p = sum(PSF .^ 2, 'all');
    if useGPU && gpuNum > 0
        PSF = gpuArray(PSF);
    end
    H = decon_psf2otf(PSF ./ sum(PSF(:)), double(sizeI));
end
clear PSF;

if gpuNum == 0
    useGPU = false;
end

% 2. Prepare parameters for iterations
if useGPU
    I = gpuArray(I);
    if deconOffset ~= 0 || EdgeErosion > 0
        I_mask = I ~= 0;
    end    
    if strcmp(underlyingType(I), 'uint16')
        I = I - uint16(Background);            
        I = single(I);
    else
        I = max(I - Background, 0);    
    end
    classI = underlyingType(I);

    J_2 = I;
    J_3 = zeros(sizeI, classI, 'gpuArray');
    J_4 = zeros(prod(sizeI), 1, classI, 'gpuArray');
else
    if ~isa(I, 'double')
        I = single(I);
    end
    if deconOffset ~= 0 || EdgeErosion > 0
        I_mask = I ~= 0;
    end    
    I = max(I - Background, 0);    
    classI = underlyingType(I);

    J_2 = I;
    J_3 = zeros(sizeI, classI);
    J_4 = zeros(prod(sizeI), 1, classI);
end

% 3. L_R Iterations
estep = saveStep;
err_mat = zeros(NUMIT, 4);

lambda = 2*any(J_4(:)~=0);
for k = 1 : NUMIT
    
    % 3.a Make an image predictions for the next iteration
    if k > 2
        % lambda = (J_4(:,1).'*J_4(:,2))/(J_4(:,2).'*J_4(:,2) +eps);        
        lambda = sum((J_2(:) - Y(:)) .* J_4) / (sum(J_4 .* J_4) + eps);
        lambda = max(min(lambda,1),0);% stability enforcement
        J_4 = J_2(:) - Y(:);
    elseif k == 2
        J_4 = J_2(:) - Y(:);
    end
    Y = max(J_2 + lambda*(J_2 - J_3), 0);% plus positivity constraint
        
    % 3.b  Make core for the LR estimation
    % directly compute CC within the same function to reduce overhead 
    ReBlurred = max(real(ifftn(H.*fftn(Y))), eps);
    % ReBlurred = ReBlurred + (ReBlurred == 0) * eps;
    ReBlurred = I./ReBlurred;
    
    % 3.c Determine next iteration image & apply positivity constraint
    J_3 = J_2;
    J_2 = max(Y.*real(ifftn(conj(H).*fftn(ReBlurred))),0);
        
    if debug
        err_mat(k, 1:2) = [k, mean((double(J_2) - double(I)) .^ 2, 'all')];
        err_mat(k, 3) = err_mat(k, 2) ./ err_mat(1, 2); 
        if rem(k, estep) == 0 && k > estep * 2
            istp = k/estep;
            err_mat(k, 4) = min(abs(err_mat(istp, 3) - err_mat(istp-1, 3)) / 10, abs(err_mat(istp, 3) + err_mat(istp-2, 3) - 2 * err_mat(istp-1, 3)) / 2);
        else
            if k < 3 * estep
                err_mat(k, 4) = 1;
            else
                err_mat(k, 4) = err_mat(istp * estep, 4);
            end
        end
        
        if k == 1
            mkdir(debug_folder);
        end
        if rem(k, estep) == 0
            if useGPU 
                writetiff(uint16(gather(J_2)), sprintf('%s/Iter_%04d.tif', debug_folder, k));
            else
                writetiff(uint16(J_2), sprintf('%s/Iter_%04d.tif', debug_folder, k));
            end
        end
    end
end

if debug
    err_mat = err_mat(1 : istp, :);
end

clear J_3 J_4 Y;

% post-processing of deconvolved result
if dampFactor > 1
    J_2 = J_2 .* (J_2 <= dampFactor * I) + min(J_2, dampFactor * I .* (J_2 >= dampFactor * I));
end
clear I;

if EdgeErosion > 0
    I_mask = decon_mask_edge_erosion(I_mask, EdgeErosion);
end

if scaleFactor ~= 1 || deconOffset ~= 0 || EdgeErosion > 0
    if deconOffset ~= 0 || EdgeErosion > 0
        J_2 = (J_2 * scaleFactor + deconOffset) .* I_mask;
    else
        J_2 = J_2 * scaleFactor;
    end
end

if save16bit
    J_2 = uint16(J_2);
end

if useGPU
    if ~isempty(bbox)
        J_2 = J_2(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6));
    end
    J_2 = gather(J_2);
else
    if ~isempty(bbox)
        try 
            J_2 = crop3d_mex(J_2, bbox);
        catch ME
            disp(ME);
            J_2 = J_2(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6));
        end
    end
end

end

