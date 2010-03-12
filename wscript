import os, sys
import Task
from TaskGen import extension

APPNAME = 'hstcal'
VERSION = '0.1'

top = '.'
out = 'bin.' + sys.platform

SUBDIRS = [
    'applib',
    'cvos',
    'hstio',
    'pkg',
    'tables',
    ]

# Add support for simple Fortran files.  This isn't complete, but it meets the
# simple .f -> .o mapping we use here.
Task.simple_task_type(
    'fortran',
    'f77 -c ${SRC} -o ${TGT}',
    color='GREEN',
    ext_out='.o',
    ext_in='.f')

@extension('.f')
def process_fortran(self, node):
        o_node = node.change_ext('.o')
        self.create_task('fortran', [node], [o_node])
        self.add_obj_file(o_node.file())

def set_options(opt):
    opt.tool_options('compiler_cc')

    # Add a required option to specify the location of CFITSIO
    opt.add_option('--cfitsio', action='store', help="The location of CFITSIO")
    
def configure(conf):
    import Options
    
    # Check for the existence of a C compiler
    conf.check_tool('compiler_cc')

    # Store and verify the location of CFITSIO
    conf.env.CFITSIO = os.path.expanduser(Options.options.cfitsio or '')
    try:
        conf.check(
            header_name='fitsio.h',
            includes=conf.env.CFITSIO,
            compile_mode='cc',
            mandatory=True,
            msg='Checking for CFITSIO')
    except:
        conf.fatal(
            "Specify the location of CFITSIO using the --cfitsio=path commandline option")

    conf.env.INCLUDEDIR = os.path.join(
        os.path.abspath(conf.srcdir), 'include') # the hstcal include directory

        
def build(bld):
    # This is to tell waf where we keep the header files
    bld.path.find_dir('include')
    
    # Recurse into all of the libraries
    for library in SUBDIRS:
        bld.recurse(library)
    
