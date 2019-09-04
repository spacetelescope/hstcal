#! /usr/bin/env python
# encoding: utf-8
#
# Provenance:  This file was copied from waf-1.6.4/waflib/Tools/gfortran.py
# and modified to run g77.  They are similar enough that only minor changes
# were needed.

import re
from waflib import Utils
from waflib.Tools import fc, fc_config, fc_scan
from waflib.Configure import conf

@conf
def find_g77(conf):
    """Find the g77 program (will look in the environment variable 'FC')"""
    fc = conf.find_program('g77', var='FC')
    fc = conf.cmd_to_list(fc)
    conf.get_g77_version(fc)
    conf.env.FC_NAME = 'g77'

@conf
def g77_flags(conf):
    v = conf.env
    v['FCFLAGS_fcshlib']   = ['-fPIC']
    v['FCFLAGS_DEBUG'] = ['-Werror'] # why not

@conf
def g77_modifier_win32(conf):
    fc_config.fortran_modifier_win32(conf)

@conf
def g77_modifier_cygwin(conf):
    fc_config.fortran_modifier_cygwin(conf)

@conf
def g77_modifier_darwin(conf):
    fc_config.fortran_modifier_darwin(conf)

@conf
def g77_modifier_platform(conf):
    dest_os = conf.env['DEST_OS'] or Utils.unversioned_sys_platform()
    g77_modifier_func = getattr(conf, 'g77_modifier_' + dest_os, None)
    if g77_modifier_func:
        g77_modifier_func()

@conf
def get_g77_version(conf, fc):
    """Get the compiler version"""

    # ensure this is actually g77, not an imposter.
    version_re = re.compile(r"GNU\s*Fortran", re.I).search
    cmd = fc + ['--version']
    out, err = fc_config.getoutput(conf, cmd, stdin=False)
    if out: match = version_re(out)
    else: match = version_re(err)
    if not match:
        conf.fatal('Could not determine the compiler type')

    # --- now get more detailed info -- see c_config.get_cc_version
    cmd = fc  # + ['-dM', '-E', '-']
    out, err = fc_config.getoutput(conf, cmd, stdin=True)

    if out.find('__GNUC__') < 0:
        conf.fatal('Could not determine the compiler type')

    k = {}
    out = out.split('\n')
    import shlex

    for line in out:
        lst = shlex.split(line)
        if len(lst)>2:
            key = lst[1]
            val = lst[2]
            k[key] = val

    def isD(var):
        return var in k

    def isT(var):
        return var in k and k[var] != '0'

    conf.env['FC_VERSION'] = (k['__GNUC__'], k['__GNUC_MINOR__'], k['__GNUC_PATCHLEVEL__'])

def configure(conf):
    conf.find_g77()
    conf.find_ar()
    conf.fc_flags()
    conf.g77_flags()
    conf.g77_modifier_platform()
