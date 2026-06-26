FN Lyr / KIC 6936115 RSP+Colors continuation run

This directory contains the large continuation run that loaded final_period100.mod
and continued toward the period-300 target, but stopped at max_model_number.

Important files:
  history.data
      Restored from archive_FN_Lyr/history_big/FN_Lyr_MESA_history_period300ish.data.gz

  final_period300.mod
      Saved model from the large continuation run.

  fn_lyr_continue_period100_to_300_MAXMODEL.log
      Full stdout log. The run stopped because model_number reached max_model_number,
      not because the model failed physically.

  FN_Lyr_KIC6936115_kepler_binned.csv
      Observed Kepler folded/binned light curve.

Known result from quick diagnostics:
  Late cycles reached photometric amplitudes around 0.17--0.24 mag,
  still below the observed FN Lyr Kepler amplitude of about 1.05 mag.
