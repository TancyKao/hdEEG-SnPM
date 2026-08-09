# Archive — project history (not user documentation)

These files are a **frozen record** of how the toolbox was built: planning notes, changelogs,
and the prompts handed to design assistants. They are kept for provenance and to explain *why*
decisions were made, but they are **not maintained** and may have drifted from the current
source. Do not follow them as instructions. Releases are tagged in the repository; see the
[release history](https://github.com/TancyKao/hdEEG-SnPM/releases).

If you want to *use* the tool, start at the [docs index](../README.md):

- **Learn it** → [tutorials](../tutorials/README.md)
- **Do a task** → [how-to guides](../how-to/README.md)
- **Look something up** → [reference](../reference/ANALYSIS_CATALOG.md)
- **Understand why** → [explanation](../explanation/README.md)

## What's here

| File | What it is | Superseded by |
|---|---|---|
| `PLAN.md` | The statistical-expansion plan + status (GLM presets, LMM path, spectral bridge), with the research-question-per-analysis and method-grounding tables. | Reference: [ANALYSIS_CATALOG](../reference/ANALYSIS_CATALOG.md). Explanation: [choosing-an-analysis](../explanation/choosing-an-analysis.md). |
| `MODIFICATIONS_SUMMARY.md` | v2.0 changelog (Pearson correlation + covariate control). **Cited line numbers have drifted from the current source** — treat as history. | — |
| `design-prompts/` | Prompts written to hand to a design assistant while evolving the report and GUI. Design history, not user docs. | The features they describe are now live in the GUI / report engines. |

## Note on line numbers

`MODIFICATIONS_SUMMARY.md` and the design prompts cite specific source line numbers. These
were accurate when written and are almost certainly stale now. When something here matters,
verify against the current `.m` file rather than trusting the cited line.

## Note on files these documents name

These documents were written while the repository also contained an internal verification
harness and a set of synthetic-data generators that wrote a `test_data/` tree. **None of those
are part of the public repository.** Where an archived
document told the author of the day to run one of them, the instruction has been reworded so it
does not send a reader after a file that is not here; the historical claim it was making is
unchanged. Ready-made example data now ships instead, at
[`example_data/`](https://github.com/TancyKao/hdEEG-SnPM/blob/main/example_data/README.md), and
that is what the paths in these documents have been repointed to.
