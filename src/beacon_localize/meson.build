executable('beacon-localize',
           ['beacon_localize.c', 'beacon_calibration.c'],
           c_args: ['-Wfloat-conversion', '-Werror', '-O3', '-Wno-unused-result', '-march=native'],
           dependencies : [common_dep],
           include_directories : src_inc,
           install : true,
           install_dir : bin_dir)
