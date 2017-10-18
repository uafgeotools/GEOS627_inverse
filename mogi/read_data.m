% unwrapped interferogram of Okmok volcano
% pixel values represent line-of-sight (LOS) range change in radians
% zeros in this interferogram represent data voids

% file must be in the local directory
filename = 'E451_20000818_20020719'; 

%file width (sample), length (line)
SAMPLE = 1100;
LINE = 980;
POSTING = 40.0;
HALF_WAVE = 28.3;

pha_file = [filename '.unw'];
obs_rngchg = read_real_be(pha_file,LINE,SAMPLE);
obs_rngchg = obs_rngchg*HALF_WAVE/2.0/pi;

% apply mask
obs_rngchg(obs_rngchg==0) = NaN;
data_mask = isnan(obs_rngchg);
data_mask = (data_mask==0);

% save mask (local run directory)
save data_mask
plot_model(obs_rngchg);
