# Contributing to DenMark

Thanks for helping improve DenMark!

- **Bugs and questions:** open an issue with one of the templates. A small example on the toy data makes a bug fast to fix. Open questions and ideas go to Discussions.
- **Changes:** open a pull request against `main` that says what changes and why, with a test.

## Development setup

- R: `devtools::install_dev_deps()`, then `devtools::test()` and `devtools::check()`. Model fits need CmdStan (`cmdstanr::install_cmdstan()`).
- Python: `pip install -e "python[test,plot,boundary]"` from the repository root, then `pytest python/tests`.

## Ground rules

- The R and Python packages build the same Stan data from the same input. A change to one needs the same change to the other, and the parity tests (`python/tests/test_parity.py`, with fixtures written by `data-raw/make_toy_data.R`) must pass.
- A change that alters a numerical result changes what the manuscript reports. Show the results on the toy data before and after the change in the pull request.
- The code blocks of the README and of the Get started pages run in the tests, so keep them runnable as written.
