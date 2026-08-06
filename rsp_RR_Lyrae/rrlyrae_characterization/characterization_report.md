# RR Lyrae RSP + Colors characterization report

## Run provenance

- Work directory: `/home/njm/MESA/MESA_Colors_Tests/rsp_RR_Lyrae`
- History file: `/home/njm/MESA/MESA_Colors_Tests/rsp_RR_Lyrae/LOGS_colors/history.data`
- Total history rows: 73162
- Analysis rows: 6518
- Model range: 131024--204185
- Completed periods in full history: 0--100
- Detected filter columns: F062, F087, F106, F129, F146, F158, F184, F213

## Main paper-ready numbers

- Median pulsation period: **0.571924 d** (window scatter 7.049e-06 d).
- Median effective temperature: **6324 K**, with peak-to-peak variation **1378 K**.
- Median luminosity: **log(L/Lsun) = 1.626**.
- Median photospheric radius: **5.337 Rsun**, with peak-to-peak variation **0.800 Rsun**. Radius source: `photosphere_r`.
- Final RSP growth/stability diagnostic: **|rsp_GREKM| = 3.797e-05** (log10=-4.421).
- Strongest/default color diagnostic: **F062-F213**, Pearson r(color, Teff) = **-0.999**.

## Filter amplitudes and phase landmarks

| filter   |   lambda_eff_A_assumed |   mag_min_bright |   mag_max_faint |   amplitude_mag |   median_mag |   phase_max_light |   phase_min_light |   rise_time_phase_min_to_max_light |   phase_lag_max_light_vs_first_filter |
|:---------|-----------------------:|-----------------:|----------------:|----------------:|-------------:|------------------:|------------------:|-----------------------------------:|--------------------------------------:|
| F062     |                   6200 |          0.8741  |          1.6672 |         0.79315 |       1.4709 |           0.54179 |           0.45253 |                           0.089263 |                               0       |
| F087     |                   8700 |          0.80713 |          1.3842 |         0.5771  |       1.1999 |           0.64726 |           0.45671 |                           0.19055  |                               0.10547 |
| F106     |                  10600 |          0.83099 |          1.3258 |         0.49478 |       1.1367 |           0.65327 |           0.45917 |                           0.1941   |                               0.11148 |
| F129     |                  12900 |          0.89662 |          1.2992 |         0.40256 |       1.1068 |           0.67607 |           0.46267 |                           0.2134   |                               0.13428 |
| F146     |                  14600 |          0.96204 |          1.3376 |         0.37561 |       1.1452 |           0.68585 |           0.46426 |                           0.22159  |                               0.14406 |
| F158     |                  15800 |          0.97174 |          1.2976 |         0.32591 |       1.1132 |           0.84952 |           0.46975 |                           0.37976  |                               0.30773 |
| F184     |                  18400 |          1.1452  |          1.4611 |         0.31585 |       1.2968 |           0.90198 |           0.47643 |                           0.42555  |                               0.36019 |
| F213     |                  21300 |          1.3846  |          1.6989 |         0.31426 |       1.5403 |           0.90198 |           0.47795 |                           0.42403  |                               0.36019 |

## Color diagnostics

| color     | blue_filter   | red_filter   |   median_color_mag |   min_color_mag |   max_color_mag |   amplitude_color_mag |   pearson_color_Teff |   pearson_color_radius |   loop_area_color_Teff |   loop_area_color_radius |
|:----------|:--------------|:-------------|-------------------:|----------------:|----------------:|----------------------:|---------------------:|-----------------------:|-----------------------:|-------------------------:|
| F062-F213 | F062          | F213         |          -0.083619 |        -0.65949 |        0.039627 |               0.69911 |             -0.999   |                0.59315 |                108.07  |                   2.7095 |
| F062-F184 | F062          | F184         |           0.15474  |        -0.41396 |        0.27544  |               0.6894  |             -0.99904 |                0.59345 |                107.37  |                   2.6664 |
| F087-F184 | F087          | F184         |          -0.10753  |        -0.46991 |       -0.029358 |               0.44055 |             -0.99896 |                0.59378 |                 69.597 |                   1.7066 |
| F106-F184 | F106          | F184         |          -0.15699  |        -0.43527 |       -0.09853  |               0.33674 |             -0.9991  |                0.59346 |                 50.262 |                   1.2978 |

## Physical summary statistics

| quantity         | unit     |    n |           min |           p05 |        median |          mean |           p95 |           max |           ptp |          std |
|:-----------------|:---------|-----:|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|-------------:|
| Teff             | K        | 6518 | 6135.7        | 6138.2        | 6324.3        | 6570.4        | 7393.9        | 7514.1        | 1378.4        | 467.75       |
| log_L            | dex Lsun | 6518 |    1.5499     |    1.5529     |    1.6263     |    1.6765     |    1.8447     |    1.8537     |    0.30381    |   0.10165    |
| L                | Lsun     | 6518 |   35.471      |   35.722      |   42.295      |   48.846      |   69.93       |   71.398      |   35.927      |  12.08       |
| radius           | Rsun     | 6518 |    4.985      |    4.9944     |    5.3375     |    5.3501     |    5.7645     |    5.7849     |    0.7999     |   0.26646    |
| log_g            | cgs      | 6518 |    2.7035     |    2.7077     |    2.7903     |    2.7843     |    2.8487     |    2.8509     |    0.14739    |   0.047994   |
| star_mass        | Msun     | 6518 |    0.65       |    0.65       |    0.65       |    0.65       |    0.65       |    0.65       |    0          |   1.1102e-16 |
| RSP period       | days     | 6518 |    0.57191    |    0.57191    |    0.57192    |    0.57192    |    0.57194    |    0.57194    |    2.4013e-05 |   7.0489e-06 |
| rsp_GREKM        |          | 6518 |    3.3828e-05 |    3.3828e-05 |    4.3026e-05 |    4.2744e-05 |    5.2477e-05 |    5.2477e-05 |    1.8649e-05 |   4.9594e-06 |
| abs_rsp_GREKM    |          | 6518 |    3.3828e-05 |    3.3828e-05 |    4.3026e-05 |    4.2744e-05 |    5.2477e-05 |    5.2477e-05 |    1.8649e-05 |   4.9594e-06 |
| max_abs_v_div_cs |          | 6518 |    0.38667    |    1.3858     |    4.1187     |    3.9342     |    6.0851     |    6.3173     |    5.9306     |   1.3283     |
| v_surf           | km/s     | 6518 |  -39.149      |  -37.058      |  -22.014      |   -8.0591     |   45.703      |   47.454      |   86.602      |  28.73       |
| rsp_DeltaR       | Rsun?    | 6518 |    0.92623    |    0.92636    |    0.92643    |    0.92643    |    0.92654    |    0.92654    |    0.00031208 |   5.0744e-05 |
| rsp_DeltaMag     | mag      | 6518 |    0.75893    |    0.75918    |    0.75939    |    0.75938    |    0.75952    |    0.75952    |    0.00059028 |   0.00010529 |
| SB_logL_residual | dex      | 6518 |   -1.032e-06  |   -1.032e-06  |   -1.032e-06  |   -1.032e-06  |   -1.032e-06  |   -1.032e-06  |    4.6629e-15 |   9.1808e-16 |

## Fourier morphology diagnostics

These are descriptive Fourier coefficients of phase-folded curves, useful for shape comparison; they are not meant as an observational calibration unless the model is tuned.

| series      |   fourier_order |   n_fit | status   |      mean |   rms_residual |        A1 |         A2 |         A3 |     R21 |      R31 |   phi21_rad |   phi31_rad |
|:------------|----------------:|--------:|:---------|----------:|---------------:|----------:|-----------:|-----------:|--------:|---------:|------------:|------------:|
| Teff        |               6 |    6518 | ok       | 6545.2    |     72.649     | 540.13    | 215.85     | 146.33     | 0.39963 | 0.27093  |      6.0737 |     5.3382  |
| radius_Rsun |               6 |    6518 | ok       |    5.4562 |      0.0015751 |   0.37349 |   0.06063  |   0.022469 | 0.16233 | 0.060159 |      3.4849 |     0.69274 |
| F062        |               6 |    6518 | ok       |    1.3129 |      0.046264  |   0.29765 |   0.13073  |   0.091737 | 0.43922 | 0.30821  |      2.349  |     4.2046  |
| F087        |               6 |    6518 | ok       |    1.0856 |      0.033955  |   0.21214 |   0.093809 |   0.065192 | 0.4422  | 0.30731  |      2.0142 |     3.6248  |
| F106        |               6 |    6518 | ok       |    1.0504 |      0.028973  |   0.17974 |   0.079068 |   0.054274 | 0.43991 | 0.30196  |      1.7717 |     3.2207  |
| F129        |               6 |    6518 | ok       |    1.0499 |      0.022838  |   0.14962 |   0.061366 |   0.040999 | 0.41014 | 0.27402  |      1.3954 |     2.5556  |

## Selected SED diagnostics

| file                  | path                                                                    |   n_wave |   wave_min_A |   wave_max_A |    flux_min |   flux_max |   flux_integral |   positive_area |   negative_area |   negative_area_fraction | has_negative_flux   | landmark       | filter_file_prefix   |   model_number |     phase |
|:----------------------|:------------------------------------------------------------------------|---------:|-------------:|-------------:|------------:|-----------:|----------------:|----------------:|----------------:|-------------------------:|:--------------------|:---------------|:---------------------|---------------:|----------:|
| F062_SED_00203205.csv | /home/njm/MESA/MESA_Colors_Tests/rsp_RR_Lyrae/SED/F062_SED_00203205.csv |     1199 |        147.2 |      1.6e+06 | -4.2832e-38 | 2.1159e-09 |      1.142e-05  |      1.142e-05  |      9.8221e-37 |               8.6008e-32 | True                | hottest_Teff   | F062                 |         203205 | 0.54131   |
| F213_SED_00203205.csv | /home/njm/MESA/MESA_Colors_Tests/rsp_RR_Lyrae/SED/F213_SED_00203205.csv |     1199 |        147.2 |      1.6e+06 | -4.2832e-38 | 2.1159e-09 |      1.142e-05  |      1.142e-05  |      9.8221e-37 |               8.6008e-32 | True                | hottest_Teff   | F213                 |         203205 | 0.54131   |
| F062_SED_00203546.csv | /home/njm/MESA/MESA_Colors_Tests/rsp_RR_Lyrae/SED/F062_SED_00203546.csv |     1199 |        147.2 |      1.6e+06 | -1.5754e-23 | 7.7382e-10 |      6.4919e-06 |      6.4919e-06 |      6.0979e-22 |               9.3931e-17 | True                | coolest_Teff   | F062                 |         203546 | 0.19217   |
| F213_SED_00203546.csv | /home/njm/MESA/MESA_Colors_Tests/rsp_RR_Lyrae/SED/F213_SED_00203546.csv |     1199 |        147.2 |      1.6e+06 | -1.5754e-23 | 7.7382e-10 |      6.4919e-06 |      6.4919e-06 |      6.0979e-22 |               9.3931e-17 | True                | coolest_Teff   | F213                 |         203546 | 0.19217   |
| F062_SED_00204185.csv | /home/njm/MESA/MESA_Colors_Tests/rsp_RR_Lyrae/SED/F062_SED_00204185.csv |     1199 |        147.2 |      1.6e+06 | -1.0459e-24 | 8.6887e-10 |      7.1944e-06 |      7.1944e-06 |      1.0971e-23 |               1.525e-18  | True                | max_radius     | F062                 |         204185 | 0.0077383 |
| F213_SED_00204185.csv | /home/njm/MESA/MESA_Colors_Tests/rsp_RR_Lyrae/SED/F213_SED_00204185.csv |     1199 |        147.2 |      1.6e+06 | -1.0459e-24 | 8.6887e-10 |      7.1944e-06 |      7.1944e-06 |      1.0971e-23 |               1.525e-18  | True                | max_radius     | F213                 |         204185 | 0.0077383 |
| F062_SED_00203927.csv | /home/njm/MESA/MESA_Colors_Tests/rsp_RR_Lyrae/SED/F062_SED_00203927.csv |     1199 |        147.2 |      1.6e+06 | -7.4566e-45 | 2.062e-09  |      1.1247e-05 |      1.1247e-05 |      1.0891e-43 |               9.6834e-39 | True                | min_radius     | F062                 |         203927 | 0.53489   |
| F213_SED_00203927.csv | /home/njm/MESA/MESA_Colors_Tests/rsp_RR_Lyrae/SED/F213_SED_00203927.csv |     1199 |        147.2 |      1.6e+06 | -7.4566e-45 | 2.062e-09  |      1.1247e-05 |      1.1247e-05 |      1.0891e-43 |               9.6834e-39 | True                | min_radius     | F213                 |         203927 | 0.53489   |
| F062_SED_00202481.csv | /home/njm/MESA/MESA_Colors_Tests/rsp_RR_Lyrae/SED/F062_SED_00202481.csv |     1199 |        147.2 |      1.6e+06 | -4.4189e-38 | 2.1161e-09 |      1.1421e-05 |      1.1421e-05 |      1.0412e-36 |               9.1169e-32 | True                | max_light_F062 | F062                 |         202481 | 0.54179   |
| F213_SED_00202481.csv | /home/njm/MESA/MESA_Colors_Tests/rsp_RR_Lyrae/SED/F213_SED_00202481.csv |     1199 |        147.2 |      1.6e+06 | -4.4189e-38 | 2.1161e-09 |      1.1421e-05 |      1.1421e-05 |      1.0412e-36 |               9.1169e-32 | True                | max_light_F062 | F213                 |         202481 | 0.54179   |
| F062_SED_00203810.csv | /home/njm/MESA/MESA_Colors_Tests/rsp_RR_Lyrae/SED/F062_SED_00203810.csv |     1199 |        147.2 |      1.6e+06 | -1.8999e-28 | 6.9323e-10 |      5.6765e-06 |      5.6765e-06 |      8.4407e-27 |               1.487e-21  | True                | min_light_F062 | F062                 |         203810 | 0.45253   |
| F213_SED_00203810.csv | /home/njm/MESA/MESA_Colors_Tests/rsp_RR_Lyrae/SED/F213_SED_00203810.csv |     1199 |        147.2 |      1.6e+06 | -1.8999e-28 | 6.9323e-10 |      5.6765e-06 |      5.6765e-06 |      8.4407e-27 |               1.487e-21  | True                | min_light_F062 | F213                 |         203810 | 0.45253   |

## SED integrity sample

- Sampled SED files: 2000
- Median negative-area fraction: 5.680e-24
- 95th percentile negative-area fraction: 9.136e-17
- Maximum negative-area fraction: 1.036e-16

## Inlist/control metadata extracted naively

### inlist
- `RSP_L` = `50d0`
- `RSP_Teff` = `6600d0`
- `RSP_X` = `0.75d0`
- `RSP_Z` = `0.0004d0`
- `RSP_kick_vsurf_km_per_sec` = `0.0d0`
- `RSP_mass` = `0.65d0`
- `RSP_max_dt` = `600d0`
- `RSP_target_steps_per_cycle` = `100`
- `colors_results_directory` = `'SED'`
- `create_RSP_model` = `.false.`
- `distance` = `3.0857d19`
- `history_columns_file` = `'history_columns_colors.list'`
- `history_interval` = `1`
- `instrument` = `'/home/njm/SED_Tools/data/filters/Roman/WFI'`
- `load_model_filename` = `'settled.mod'`
- `load_saved_model` = `.true.`
- `max_model_number` = `4000000`
- `profile_columns_file` = `'profile_columns.list'`
- `save_model_filename` = `'colors_final.mod'`
- `save_model_when_terminate` = `.true.`
- `sed_per_model` = `.true.`
- `stellar_atm` = `'/home/njm/SED_Tools/data/stellar_models/Kurucz2003all/'`
- `trace_history_value_name(1)` = `'rel_E_err'`
- `trace_history_value_name(2)` = `'log_rel_run_E_err'`
- `use_colors` = `.true.`
- `vega_sed` = `'data/colors_data/stellar_models/vega_flam.csv'`
- `x_integer_ctrl(1)` = `2`
- `x_integer_ctrl(5)` = `100`

### inlist_rsp_RR_Lyrae
- `RSP_L` = `50d0`
- `RSP_Teff` = `6600d0`
- `RSP_X` = `0.75d0`
- `RSP_Z` = `0.0004d0`
- `RSP_kick_vsurf_km_per_sec` = `20.0d0`
- `RSP_mass` = `0.65d0`
- `RSP_max_dt` = `600d0`
- `RSP_target_steps_per_cycle` = `150`
- `create_RSP_model` = `.true.`
- `history_columns_file` = `'history_columns_settle.list'`
- `history_interval` = `10`
- `initial_model_number` = `0`
- `max_model_number` = `4000000`
- `profile_columns_file` = `'profile_columns.list'`
- `save_model_filename` = `'settled.mod'`
- `save_model_when_terminate` = `.true.`
- `sed_per_model` = `.false.`
- `set_initial_model_number` = `.true.`
- `trace_history_value_name(1)` = `'rel_E_err'`
- `trace_history_value_name(2)` = `'log_rel_run_E_err'`
- `use_colors` = `.false.`
- `x_ctrl(2)` = `1d-3`
- `x_integer_ctrl(1)` = `1`
- `x_integer_ctrl(2)` = `80`
- `x_integer_ctrl(3)` = `8`
- `x_integer_ctrl(4)` = `600`

### inlist_rsp_RR_Lyrae_colors
- `RSP_L` = `50d0`
- `RSP_Teff` = `6600d0`
- `RSP_X` = `0.75d0`
- `RSP_Z` = `0.0004d0`
- `RSP_kick_vsurf_km_per_sec` = `0.0d0`
- `RSP_mass` = `0.65d0`
- `RSP_max_dt` = `600d0`
- `RSP_target_steps_per_cycle` = `100`
- `colors_results_directory` = `'SED'`
- `create_RSP_model` = `.false.`
- `distance` = `3.0857d19`
- `history_columns_file` = `'history_columns_colors.list'`
- `history_interval` = `1`
- `instrument` = `'/home/njm/SED_Tools/data/filters/Roman/WFI'`
- `load_model_filename` = `'settled.mod'`
- `load_saved_model` = `.true.`
- `max_model_number` = `4000000`
- `profile_columns_file` = `'profile_columns.list'`
- `save_model_filename` = `'colors_final.mod'`
- `save_model_when_terminate` = `.true.`
- `sed_per_model` = `.true.`
- `stellar_atm` = `'/home/njm/SED_Tools/data/stellar_models/Kurucz2003all/'`
- `trace_history_value_name(1)` = `'rel_E_err'`
- `trace_history_value_name(2)` = `'log_rel_run_E_err'`
- `use_colors` = `.true.`
- `vega_sed` = `'data/colors_data/stellar_models/vega_flam.csv'`
- `x_integer_ctrl(1)` = `2`
- `x_integer_ctrl(5)` = `100`

## Suggested wording/caveats

- Treat this as a time-dependent synthetic-photometry demonstration from a self-consistent RSP calculation, not as a calibrated fit to a particular observed RR Lyrae.
- If max |v/cs| is large, describe the model as shocky/supersonic in some phases rather than pretending the detailed morphology is tuned.
- The robust quantities for the section are the period, Teff/radius/luminosity ranges, color-temperature anti-correlation, wavelength-dependent amplitude trend, and phase-dependent SED/filter response.
