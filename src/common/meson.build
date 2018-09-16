common_src = ['getopt.c', 'string_util.c', 'time_util.c', 'zarray.c', 'zhash.c']

common_lib = library('common',
                     common_src,
                     dependencies : [math_dep],
                     include_directories : src_inc,
                     install : true,
                     install_dir : lib_dir)

common_dep = declare_dependency(link_with : common_lib,
                                include_directories : src_inc)
