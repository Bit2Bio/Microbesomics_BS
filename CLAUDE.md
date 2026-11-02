# Microbesomics_BS — Claude instructions

**GitHub repo**: https://github.com/Bit2Bio/Microbesomics_BS

Read `README.md` at the start of every session for full pipeline documentation:
scripts, run commands, parameters, input/output paths, and analysis results.

## Rules

- **No pipeline changes without approval.** Never modify the analysis pipeline (scripts, Snakefile, config, statistical models, filters) without explicitly proposing the change and waiting for confirmation.
- **Never appear as co-author.** Do not add `Co-Authored-By: Claude` or any Claude attribution to git commit messages.
- **Comment the WHY.** When writing or editing scripts, add comments that explain non-obvious choices: parameter values, data quirks, model design decisions, edge cases. Do not comment the obvious.
- **Absolute paths only.** All scripts use paths anchored to `/home/lorenzo/Microbesomics/`. Never use relative paths.
- **Workflow structure mirrors Atheromics_TVAS.** Scripts live in `workflow/scripts/`, pipeline is defined in `workflow/Snakefile`, parameters in `config/config.yaml`.
