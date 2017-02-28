# vim: set syntax=python:

import os, platform, shutil, sys

from waflib import Configure
from waflib import Errors
from waflib import Logs
from waflib import Options
from waflib import Scripting
from waflib import Task
from waflib import Utils
from waflib import TaskGen

APPNAME = 'hstcal'
VERSION = '0.1.1'

top = '.'
out = 'build.' + platform.platform()

# A list of subdirectories to recurse into
SUBDIRS = [
    'cfitsio', # cfitsio needs to go first
    'applib',
    'cvos',
    'hstio',
    'hstio/test',
    'include',
    'pkg',
    'tables',
    ]

# Have 'waf dist' create tar.gz files, rather than tar.bz2 files
Scripting.g_gz = 'gz'

# Have gcc supersede clang
from waflib.Tools.compiler_c import c_compiler
c_compiler['darwin'] = ['gcc', 'clang']
c_compiler['default'] = ['gcc', 'clang']

option_parser = None
def options(opt):
    # We want to store
    # the option parser object so we can use it later (during the
    # configuration phase) to parse options stored in a file.
    global option_parser
    option_parser = opt.parser

    opt.load('compiler_c')
    opt.load('compiler_fc')

    opt.add_option(
        '--disable-openmp', action='store_true',
        help="Disable OpenMP")

    opt.add_option(
        '--debug', action='store_true',
        help="Create a debug build")

    opt.recurse('cfitsio')

def _setup_openmp(conf):
    conf.start_msg('Checking for OpenMP')

    if conf.options.disable_openmp:
        conf.end_msg('OpenMP disabled.', 'YELLOW')
        return

    try:
        conf.check_cc(
            header_name="omp.h", lib="gomp", cflags="-fopenmp",
            uselib_store="OPENMP")
    except Errors.ConfigurationError:
        conf.end_msg("OpenMP not found.", 'YELLOW')
    else:
        conf.end_msg("OpenMP found.", 'GREEN')

def _check_mac_osx_version(floor_version):
    '''
    Purpose:
        Converts the semantic version of the operating system to a 24-bit integer, then
        compares the result against the user-defined `floor_version`. Returns true if
        `osx_version` is greater than or equal to `floor_version`.

    Summary:
        Checks whether the operation system meets the minimum version requirements.

    Example:
        # Assume `sw_vers -ProductVersion` returns 10.11.0

        # Is 10.11.0 (0x0A0B00) greater than 10.5.0 (0x0A0500)?
        >>> _determine_mac_osx_floor(0x0A0500)
        True

        # Is 10.11.0 (0x0A0B00) equal to 10.11.0 (0x0A0B00)?
        >>> _determine_mac_osx_floor(0x0A0B00)
        True

        # Is 10.11.0 (0x0A0B00) greater than 10.12.1 (0x0A0C01)?
        >>> _determine_mac_osx_floor(0x0A0C01)
        False

    Encoding:
        OS Version      Encoded Version
        -----------     ---------------
        10.5.0      ==  0x0A0500
         ^ ^ ^             ^ ^ ^
         | | |             | | |
         major             major
           | |               | |
           minor             minor
             |                 |
             patch             patch

    '''

    assert isinstance(floor_version, int)
    s = platform.popen("/usr/bin/sw_vers -productVersion").read()

    # Extract the integer values between the '.'s
    osx_version_major, osx_version_minor, osx_version_patch = tuple(int(x) for x in s.strip().split('.'))

    # Convert major/minor/patch values into a single 24-bit integer
    osx_version = (osx_version_major & 0xff) << 16 | (osx_version_minor & 0xff) << 8 | (osx_version_patch & 0xff )

    # If the operating system version does meet or exceed the minimum
    if osx_version < floor_version:
        return False

    return True

def _determine_mac_osx_fortran_flags(conf):
    # On Mac OS-X, we need to know the specific version in order to
    # send some compile flags to the Fortran compiler.
    import platform
    conf.env.MAC_OS_NAME = None
    if platform.system() == 'Darwin' :
        conf.start_msg('Checking Mac OS-X version')

        if _check_mac_osx_version(0x0A0500):
            conf.end_msg('done', 'GREEN')
        else:
            conf.end_msg(
                "Unsupported OS X version detected (<10.5.0)",
                'RED')
            exit(1)

        conf.env.append_value('FCFLAGS', '-m64')

def _determine_sizeof_int(conf):
    conf.check(
        fragment='#include <stdio.h>\nint main() { printf("%d", sizeof(int)); return 0; }\n',
        define_name="SIZEOF_INT",
        define_ret=True,
        quote=False,
        execute=True,
        msg='Checking for sizeof(int)')

def configure(conf):
    # NOTE: All of the variables in conf.env are defined for use by
    # wscript files in subdirectories.

    # Read in options from a file.  The file is just a set of
    # commandline arguments in the same syntax.  May be spread across
    # multiple lines.
    if os.path.exists('build.cfg'):
        fd = open('build.cfg', 'r')
        for line in fd.readlines():
            tokens = line.split()
            options, args = option_parser.parse_args(tokens)
            for key, val in options.__dict__.items():
                if Options.options.__dict__.get(key) is None:
                    Options.options.__dict__[key] = val
        fd.close()

    # Load C compiler support
    conf.load('compiler_c')

    # Check for the existence of a Fortran compiler
    conf.load('compiler_fc')
    conf.check_fortran()

    # Set the location of the hstcal include directory
    conf.env.INCLUDES = [os.path.abspath('include')] # the hstcal include directory

    # A list of the local (hstcal) libraries that are typically linked
    # with the executables
    conf.env.LOCAL_LIBS = [
        'applib', 'xtables', 'hstio', 'cvos', 'CFITSIO']

    # A list of external libraries that are typically linked with the
    # executables
    conf.env.EXTERNAL_LIBS = ['m']
    if sys.platform.startswith('sunos'):
        conf.env.EXTERNAL_LIBS += ['socket', 'nsl']

    # A list of paths in which to search for external libraries
    conf.env.LIBPATH = []

    _determine_mac_osx_fortran_flags(conf)

    _setup_openmp(conf)

    _determine_sizeof_int(conf)

    # check whether the compiler supports -02 and add it to CFLAGS if it does
    if conf.options.debug:
        if conf.check_cc(cflags='-g'):
            conf.env.append_value('CFLAGS', '-g')
        if conf.check_cc(cflags='-O0'):
            conf.env.append_value('CFLAGS', '-O0')
        if conf.check_cc(cflags='-Wall'):
            conf.env.append_value('CFLAGS','-Wall')
    else:
        if conf.check_cc(cflags='-O2'):
            conf.env.append_value('CFLAGS','-O2')
        if conf.check_cc(cflags='-Wall'):
            conf.env.append_value('CFLAGS','-Wall')
        if conf.check_cc(cflags='-fstack-protector-all'):
            conf.env.append_value('CFLAGS','-fstack-protector-all')

    conf.start_msg('C compiler flags (CFLAGS)')
    conf.end_msg(' '.join(conf.env['CFLAGS']) or None)

    conf.start_msg('Fortran compiler flags (FCFLAGS)')
    conf.end_msg(' '.join(conf.env['FCFLAGS']) or None)

    conf.start_msg('Linker flags (LDFLAGS)')
    conf.end_msg(' '.join(conf.env['LDFLAGS']) or None)

    # The configuration related to cfitsio is stored in
    # cfitsio/wscript
    conf.recurse('cfitsio')

def build(bld):
    bld(name='lib', always=True)
    bld(name='test', always=True)

    targets = [x.strip() for x in bld.targets.split(',')]

    if not len(targets):
        targets = ['lib', 'test']

    if 'lib' in targets:
        bld.env.INSTALL_LIB = True
        targets.remove('lib')
    else:
        bld.env.INSTALL_LIB = None

    if 'test' in targets:
        bld.env.INSTALL_TEST = True
        targets.remove('test')
    else:
        bld.env.INSTALL_TEST = None

    bld.targets = ','.join(targets)

    # Recurse into all of the libraries
    for library in SUBDIRS:
        bld.recurse(library)

    if bld.cmd == 'clean':
        return clean(bld)

    # Add a post-build callback function
    bld.add_post_fun(post_build)

def post_build(bld):
    # WAF has its own way of dealing with build products.  We want to
    # emulate the old stsdas way of creating a flat directory full of
    # .a and .e files.  This simply runs through the build tree and
    # copies such files to the bin.* directory.
    src_root = os.path.join(
        bld.srcnode.abspath(),
        'build.' + platform.platform())
    dst_root = os.path.join(
        bld.srcnode.abspath(),
        'bin.' + platform.platform())

    if not os.path.exists(dst_root):
        os.mkdir(dst_root)

    for root, dirs, files in os.walk(src_root):
        for file in files:
            base, ext = os.path.splitext(file)
            if ext in ('.e', '.a'):
                src_path = os.path.join(root, file)
                dst_path = os.path.join(dst_root, file)
                shutil.copy(src_path, dst_path)

def clean(ctx):
    # Clean the bin.* directory
    bin_root = 'bin.' + platform.platform()
    if os.path.exists(bin_root):
        shutil.rmtree(bin_root)

def test(ctx):
    # Recurse into all of the libraries
    # Just to check that nose is installed
    try:
        import nose
    except ImportError:
        raise ImportError("nose must be installed to run hstcal tests.")

    for library in SUBDIRS:
        if library.endswith('test'):
            if os.system('nosetests %s' % library):
                raise Exception("Tests failed")

# This is a little recipe from the waf docs to produce abstract
# targets that define what to build and install.
from waflib.TaskGen import feature, before_method
@feature("*")
@before_method('process_rule')
def post_the_other(self):
    deps = getattr(self, 'depends_on', [])
    for name in self.to_list(deps):
        other = self.bld.get_tgen_by_name(name)
        other.post()
