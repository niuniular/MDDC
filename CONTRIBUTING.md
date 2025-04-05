# Contributing to MDDC

This CONTRIBUTING.md is adapted from [Peter Desmet's guide](https://gist.github.com/peterdesmet/e90a1b0dc17af6c12daf6e8b2f044e7c).

First of all, thanks for considering contributing to MDDC!

- [Repository](https://github.com/niuniular/MDDC)
- [Issues](https://github.com/niuniular/MDDC/issues)
- [New Issue](https://github.com/niuniular/MDDC/issues/new)
- [Website](https://niuniular.github.io/MDDC/)
- [Code of Conduct](https://github.com/niuniular/MDDC/blob/main/CODE_OF_CONDUCT.md)
- [Bug Report](https://github.com/niuniular/MDDC/issues/new?assignees=&labels=Bug%2CNeeds+Triage&projects=&template=bug_report.yml)
- [Documentation Improvement](https://github.com/niuniular/MDDC/issues/new?assignees=&labels=Documentation%2CNeeds+Triage&projects=&template=documentation_improvement.yml)
- [Email](mailto:raktimmu@buffalo.edu)

## Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/niuniular/MDDC/blob/main/CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

## How You Can Contribute

There are several ways you can contribute to this project.

### Share the Love ❤️

Think MDDC is useful? Let others discover it, by telling them in person, via Twitter, ResearchGate, or a blog post.

Using MDDC for a paper you are writing? Consider [citing it](https://niuniular.github.io/MDDC/authors.html#citation).

### Ask a Question ⁉️

Using MDDC and got stuck? Browse the [documentation](https://niuniular.github.io/MDDC/) to see if you can find a solution. Still stuck? Post your question as an [issue on GitHub](https://github.com/niuniular/MDDC/issues/new). While we cannot offer user support, we'll try to do our best to address it, as questions often lead to better documentation or the discovery of bugs.

Want to ask a question in private? Contact the package maintainer by [mail](mailto:raktimmu@buffalo.edu).

### Propose an Idea 💡

Have an idea for a new MDDC feature? Take a look at the [documentation](https://niuniular.github.io/MDDC/) and [issue list](https://github.com/niuniular/MDDC/issues) to see if it isn't included or suggested yet. If not, suggest your idea as an [issue on GitHub](https://github.com/niuniular/MDDC/issues/new). While we can't promise to implement your idea, it helps to:

- Explain in detail how it would work.
- Keep the scope as narrow as possible.

See below if you want to contribute code for your idea as well.

### Report a Bug 🐛

Using MDDC and discovered a bug? That's annoying! Don't let others have the same experience and report it as an [issue on GitHub](https://github.com/niuniular/MDDC/issues/new) so we can fix it. A good bug report makes it easier for us to do so, please try to give as much detail as possible at [Bug Report](https://github.com/niuniular/MDDC/issues/new?assignees=&labels=Bug%2CNeeds+Triage&projects=&template=bug_report.yml).

### Improve the Documentation 📖

Noticed a typo on the website? Think a function could use a better example? Good documentation makes all the difference, so your help to improve it is very welcome! Submit an issue here [Documentation Improvement](https://github.com/niuniular/MDDC/issues/new).

## API Documentation

The API documentation is built automatically from the docstrings of classes, functions, etc. in the source files. The docs are built with `devtools::build_manual()`. 

- Go to the `R/` or `src/` directory in the [code repository](https://github.com/niuniular/MDDC).
- Look for the file with the name of the function.
- [Propose a file change](https://help.github.com/articles/editing-files-in-another-user-s-repository/) to update the function documentation in the Roxygen comments (starting with `#'`).

## Contribute Code 📝

Care to fix bugs or implement new functionality for MDDC? Awesome! 👏 Have a look at the [issue list](https://github.com/niuniular/MDDC/issues) and leave a comment on the things you want to work on. See also the development guidelines below.

## Development Guidelines

We try to follow the [GitHub flow](https://guides.github.com/introduction/flow/) for development.

1. Fork [this repo](https://github.com/niuniular/MDDC) and clone it to your computer. To learn more about this process, see [this guide](https://guides.github.com/activities/forking/).
2. If you have forked and cloned the project before and it has been a while since you worked on it, [pull changes from the original repo](https://help.github.com/articles/merging-an-upstream-repository-into-your-fork/) to your clone by using `git pull upstream master`.
3. Open the folder on your local machine using any code editor.
4. Make your changes:
   - Write your code.
   - Test your code (bonus points for adding unit tests).
   - Document your code (see function documentation above).
   - Check your code with `pytest`.
5. Commit and push your changes.
6. Submit a [pull request](https://guides.github.com/activities/forking/#making-a-pull-request).

## Future Developments

Regular updates and bug fixes are planned to continually enhance the package's functionality and user experience. One of our primary goals is to make MDDC increasingly user-friendly, with improvements to the user experience and the layout of the outputs. User feedback is highly valued and will be a key driver of future development. This Life Cycle Statement is subject to periodic review and will be updated to reflect the evolving nature of MDDC.