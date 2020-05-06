# Contributing guidelines

HTStream was designed to quickly and easily add additional applications to the suite by using a centrally maintained core of fast C++ read processing code.

## Pull Request Checklist

Before sending your pull requests, make sure you followed this list.

- Read [Contributing Guidelines](CONTRIBUTING.md).
- Read [Code of Conduct](CODE_OF_CONDUCT.md).
- Check if my changes are consistent with the [guidelines](#general-guidelines-and-philosophy-for-contribution).
- Run 'make test' and ensure all the Unit tests pass.

## How to become a contributor and submit your own code

We'd love to accept your additions!

***NOTE***: Only original source code from you should be submitted for acceptance into the main repository.

### Contributing code

If you want to contribute, start working through the HTStream codebase, navigate to the [Github "issues" tab](https://github.com/ibest/HTStream/issues) and start looking through interesting issues.
[issues with the "enhancement" label](https://github.com/ibest/HTStream/labels/enhancement) are well suited for outside contributions, and we haven't had time yet to get to them.
**If you have an idea for a new feature or application**, please first create a new issue (add the 'enhancement' or 'feature' label) and get feedback from other developers and users. We would love to give guidance, save duplicate work if the feature is already under development, or save you some time if the functionality lives elsewhere.

If you decide to start on an issue, leave a comment so that other people know that you're working on it. If you want to help out, but not alone, use the issue comment thread to coordinate.

When you have additions, or improvements to HTStream, send us your pull requests! For those just getting started, Github has a
[how to](https://help.github.com/articles/using-pull-requests/).

HTStream team members will be assigned to review your pull requests.
Once the pull requests are approved and pass continuous integration checks, a HTStream team member will apply `ready to pull` label to your change.
This means we are working on getting your pull request submitted to our internal repository.
After the change has been submitted internally, your pull request will be merged automatically on GitHub.

### Contribution guidelines and standards

Before sending your pull request for
[review](https://github.com/ibest/HTStream/pulls),
make sure your changes are consistent with the guidelines.

#### General guidelines and philosophy for contribution

*   Include unit tests when you contribute new features, as they help to a) prove that your code works correctly, and b) guard against future breaking changes to lower the maintenance cost.
*   Bug fixes also generally require unit tests, because the presence of bugs usually indicates insufficient test coverage.
*   New applications are designed following the
["philosophy of Program Design in the UNIX Environment"](https://onlinelibrary.wiley.com/doi/abs/10.1002/j.1538-7305.1984.tb00055.x).
*   When you contribute a new feature/application to HTStream, the maintenance burden is (by default) transferred to the HTStream team. This means that the benefit of the contribution must be compared against the cost of maintaining the feature.
