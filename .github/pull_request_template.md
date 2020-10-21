# luslab-nf-modules

Thank you for contributing and submitting a pull request!

Please check you have done the following:

## Uniform coding style

- [ ] Used snake_case rather than camelCase
- [ ] Followed the output channel naming convention (e.g. `bam` for BAM files, `bed` for BED files)

## Continuous integration

- [ ] Added any new modules to `.github/workflows/modules-testing.yml`
- [ ] Added any new containers to `.github/workflows/docker-build-push.yml` and `.github/workflows/docker-linting.yml`
- [ ] Added any new workflows to `.github/workflows/workflow-testing.yml`
- [ ] Added any new test data to `test_data`
- [ ] When changing an existing module, have you bumped the version? Major.Minor.Patch - Major for breaking changes, Minor for feature updates and most things, Patch for bug fixes and small things like comment and variable name updates
