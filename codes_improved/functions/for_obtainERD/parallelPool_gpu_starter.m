function parallelPool_gpu_starter(config)
% Start parallel pool if not already running
if isempty(gcp('nocreate'))
    parpool('local', feature('numcores'));
    fprintf('Started parallel pool with %d workers\n', feature('numcores'));
end

if config.useGPU
    gpuDevice(1); % Select the first GPU
    fprintf('Using GPU for acceleration\n');
else
    fprintf('No GPU available, using CPU only\n');
end