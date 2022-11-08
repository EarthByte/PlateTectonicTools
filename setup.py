## To install locally: python -m pip install .
## (If there are problems with installation of the documentation, it may be that
##  the egg file is out of sync and will need to be manually deleted - see error message
##  for details of the corrupted zip file. )
##
## To push a version through to pip.
##  - Make sure it installs correctly locally as above
##  - Update the version information in this file
##  - python setup.py sdist
##  - python -m twine upload -r testpypi dist/*  # for the test version
##  - python -m twine upload -r pypi     dist/*  # for the real version
##
## (see http://peterdowns.com/posts/first-time-with-pypi.html)


from setuptools import setup
from os import path
import io

## in development set version to none and ...
PYPI_VERSION = "0.4"

# Return the git revision as a string (from numpy)
def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout = subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', '--short', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION


if PYPI_VERSION is None:
    PYPI_VERSION = git_version()

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md')) as f:
    long_description = f.read()


if __name__ == "__main__":
    setup(name              = 'PlateTectonicTools',
          long_description  = long_description,
          long_description_content_type = 'text/markdown',
          author            = "John Cannon",
          author_email      = "john.cannon@sydney.edu.au",
          url               = "https://github.com/EarthByte/PlateTectonicTools",
          version           = PYPI_VERSION,
          description       = "Python tools for plate tectonic research",
          install_requires  = ['numpy'],
          python_requires   = '>=2.7',
          packages          = ['ptt', 'ptt.utils'],
          package_data      = {'ptt': ['Examples/notebooks/*.ipynb',
                                       'Examples/data/*.gpmlz',
                                       'Examples/data/*.rot',
                                       'utils/*.py']},
          classifiers       = ['Programming Language :: Python :: 2',
                               'Programming Language :: Python :: 2.7',
                               'Programming Language :: Python :: 3']
          )
