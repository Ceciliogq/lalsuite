# LALSuite

This is the main LALSuite development repository.
If you would like to just use a released version,
see [here](https://computing.docs.ligo.org/conda/)
for IGWN conda environments including LALSuite,
or the project pages
on [conda-forge](https://anaconda.org/conda-forge/lalsuite)
and on [PyPI](https://pypi.org/project/lalsuite/).

## Acknowledgment

We request that any academic report, publication, or other academic
disclosure of results derived from the use of this software acknowledge
the use of the software by an appropriate acknowledgment or citation.

The whole LALSuite software suite can be cited with the DOI
[10.7935/GT1W-FZ16][doi]. An example BibTeX entry could look like this:

     @misc{lalsuite,
           author         = "{LIGO Scientific Collaboration}",
           title          = "{LIGO} {A}lgorithm {L}ibrary - {LALS}uite",
           howpublished   = "free software (GPL)",
           doi            = "10.7935/GT1W-FZ16",
           year           = "2018"
     }

In addition, some codes contained in this package may be directly based
on one or several scientific papers, which should be cited when using
those specific codes; some of these can be discovered through the
documentation.

## Cloning the Repository

We now utilize [Git LFS][gitlfs] for the managament of large files and
as such `git-lfs` needs to be installed and configured to correctly
clone this repository. After installing `git-lfs` it can be configured
using:

     $ git lfs install

This only needs to be done once for each machine you access the
repository. It can then be cloned using:

     $ git clone git@git.ligo.org:lscsoft/lalsuite.git

## Building from Source

The recommended way to build LALSuite from source is in a `conda` environment.
[A recipe file](conda/environment.yml) is available with all main dependencies.
This can serve as the base for custom recipes,
or be used directly via:

     $ conda env create -f conda/environment.yml

Pulling in dependencies may take a while depending on your internet connection.
After the environment setup succeeded, you can activate it with:

     $ conda activate lalsuite-dev

You can then build the suite by executing, in order:
1. `./00boot` (once at first time)
2. `./configure` with appropriate options (see `./configure --help`)
3. `make`

After pulling updates or making your own changes,
you will usually only need to call `make` again,
as reconfiguration and re-running `00boot` should be handled automatically if needed.

If you prefer managing dependencies yourself without conda,
see [here](https://wiki.ligo.org/Computing/LALSuite#Dependencies).

## Contributing to LALSuite

The [guide to contributing to LALSuite][contributing] explains how to
report issues and contribute fixes or new features using the fork and
pull workflow. Please read and follow these directions.

## Nightly Documentation

The [LALSuite Doxygen documentation][nightlydocs] is built under
GitLab-CI every night.

## Notes on Ancient History

LALSuite was transferred to `git.ligo.org` in December 2017. Older
history has been imported, though commit hashes were rewritten during
the [Git LFS][gitlfs] conversion. Please note:

1. The `Original:` commit IDs quoted in each commit message can be used
   to compare with the [archived reference repo][oldlalsuite], old issue
   discussions on the [Redmine tracker][oldredmine], review wiki pages
   etc.

1. Commits before December 2017 may also include references to issues
   (`#number`). These refer to the corresponding [Redmine
   issue][oldredmine] (LVC-authorized access only), and any clickable
   link the internal GitLab web interface produces for those old commits
   will therefore be spurious.

## For More Information

Please visit the [LALSuite project page][project].

[doi]:          https://doi.org/10.7935/GT1W-FZ16
[gitlfs]:       https://wiki.ligo.org/Computing/GitLFS#Install_the_git_LFS_client
[contributing]: https://git.ligo.org/lscsoft/lalsuite/blob/master/CONTRIBUTING.md
[nightlydocs]:  https://lscsoft.docs.ligo.org/lalsuite
[oldlalsuite]:  https://git.ligo.org/lscsoft/lalsuite-archive
[oldredmine]:   https://bugs.ligo.org/redmine/projects/lalsuite
[project]:      https://wiki.ligo.org/Computing/LALSuite
