project('non-parametric-uwb', 'c',
        default_options : ['c_std=gnu99',
                           'default_library=static',
                           'buildtype=release'])

compiler = meson.get_compiler('c')

math_dep = compiler.find_library('m', required: false)

src_inc = include_directories(['src'])
bin_dir = meson.source_root() + '/bin'
lib_dir = meson.source_root() + '/lib'

subdir('src')
