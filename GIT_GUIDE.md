GATE currently depends on a locally prepared runtime and preserved dependency
archives. Cloning this repository alone does not recreate that runtime. See
GATE_PILOT.md and the methodology source for its installation record.

## Inspect and record a change

Run these commands in PowerShell:

```powershell
Set-Location D:\survival_gwas_simulation
git status
git diff
git add scripts/report_gate.R scripts/gate_methodology.R
git diff --cached
git commit -m "Expand GATE methodology documentation"
git log --oneline -5
```

`status` shows changed files. `diff` shows changes not yet staged. `add` selects
the changes for the next snapshot. `diff --cached` reviews that selection.
`commit` records it with an explanation. `log` lists previous snapshots.
For the first commit, stage the intended source files and use a message such as
"Add simulated survival GWAS pipeline and GATE analysis".

## Connect to GitHub and push

Create an empty GitHub repository, without adding a README, licence or gitignore
there when connecting this existing project. Choose public for a visible
portfolio or private for limited access. Substitute its actual URL below:

```powershell
git remote add origin https://github.com/YOUR_USERNAME/YOUR_REPOSITORY.git
git remote -v
git push -u origin main
```

`origin` names the remote repository. `push` uploads commits. `-u` connects your
local main branch to the remote main branch, so later updates can use `git push`.
Git may open its normal authentication flow. Do not put passwords or tokens in
scripts or repository files. If the repository already has commits, inspect its
history before connecting; do not force push to overwrite it.

## Try a change on a branch

```powershell
git switch -c learning/update-readme
# Edit README.md, then inspect and commit the change.
git add README.md
git commit -m "Clarify project setup"
git switch main
git merge learning/update-readme
git push
```

A branch lets you work on a change separately. A merge brings its commits into
main. Before returning to main, commit or otherwise preserve your current edits.
On later collaborative work, fetch remote changes with `git fetch origin` and
inspect them before merging. Prefer explicit file names when staging while you
are learning, so you can see exactly what each commit contains.
