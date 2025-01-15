# Release Process

## Github
- create a release branch `release-1.4.0`
- Update the version here:  `common/src/version.h.release` e.g. `1.4.0`
- commit the change
- merge to master
- tag that commit `git tag v1.4.0` then `git push origin tag v1.4.0`
- https://github.com/s4hts/HTStream/releases/new create a new release here

## Bioconda
- edit `bioconda-recipes/edit/master/recipes/htstream/meta.yaml` in a fork of the bioconda repo.
  - set version `1.4.0`
  - set sha256 by downloading the file and run `sha256sum v1.4.0.tar.gz`
