# Contributing

## Git Workflow

This project uses a version of [GitLab Flow](https://about.gitlab.com/topics/version-control/what-is-gitlab-flow/) for its workflow. This uses two primary long-lived branches, `main` and `release`, plus short-lived feature and release candidate branches:

- `main` - This is the default and main working branch of the project. Day-to-day changes start from this branch and are merged back into it when ready. `main` is protected so that it can only be merged to via Pull Request.
- `release` - This is the branch that is used for deployments of the project. It is periodically updated by a release candidate branch being created from `main`, tested, and merged into this when approved. `release` is protected so that it can only be merged to via Pull Request.
- Feature branches - These are used for day-to-day work. Feature branches are branched off of `main` and then worked on (usually by a single developer). They may be discarded if it's decided they are unneeded, or else if they're approved, they're merged back into `main` via a Pull Request. To avoid issues which can arise from feature branches being open too long and diverging from `main`, it is a good idea to aim for each to only represent a small, discrete change. A good target is approximately one day's worth of work, and only exceptionally should a feature branch be open for more than a week. The recommended naming convention for feature branches is `<initials>-<description-of-change>`, e.g. `brg-new-feature`.
- Release candidate branches - These branches are used to bridge `main` and `release`. When it's time to prepare for a release, a release candidate branch is branched off of `main` to isolate it from any further changes. It then undergoes full testing, including manual testing of the web app. If any issues are identified which require fixing, the branch can be updated until all tests pass. Once this is the case, the branch is merged into `release` as well as back into `main` (so that any fixes made in it are also reflected there). The recommended naming convention for release candidate branches is `rc-<target-version-number>`, e.g. `rc-0.1.4`.

## Editing Advice

For debugging python issues, it's recommended to install the package in editable mode via pip. This sets it up so that the python source files are used in-place rather than copied to a separate install directory, meaning that changes to them will be reflected in runs without need for a new installation. This can be done through the following command (which also installs all optional packages):

```bash
pip install --editable .[gui][test]
```
