"""Installs sample files to user specified directory"""
from os.path import expanduser, join
from shutil import copytree
from pkg_resources import Requirement, resource_filename

def install_samples():
    """Installs sample files to user specified directory"""

    print('This will install the poropyck samples.')
    print()
    print('Please enter installation directory for sample files.')
    dst = join(expanduser(input('(default - ' + expanduser('~') + '):')), 'poropyck_samples')
    src = resource_filename(Requirement.parse('poropyck'), 'poropyck/samples')
    copytree(src, dst)
