# Push to GitHub - Manual Instructions

Follow these steps to push this repository to your GitHub account.

## Option 1: Quick Setup (Automated)

Run the setup script:

```bash
cd /Users/jjanssens/Documents/Postdoc/Projects/evolution/crested_mod
./setup_github.sh
```

Then follow the instructions printed by the script.

---

## Option 2: Manual Setup (Step by Step)

### Step 1: Update README with your GitHub username

Replace `YOUR_USERNAME` with your actual GitHub username in:
- `GITHUB_README.md`
- `README.md` (if using full version docs)

You can do this with sed:
```bash
sed -i '' 's/YOUR_USERNAME/jjanssens/g' GITHUB_README.md
sed -i '' 's/YOUR_USERNAME/jjanssens/g' README.md
```

Or just open the files and search/replace manually.

### Step 2: Initialize Git

```bash
cd /Users/jjanssens/Documents/Postdoc/Projects/evolution/crested_mod

# Initialize git repository
git init

# Add all files
git add .

# Create initial commit
git commit -m "Initial commit: CREsted utilities for cross-species analysis"
```

### Step 3: Create GitHub Repository

1. Go to GitHub: https://github.com/new
2. Repository name: `crested-utils` (or your preferred name)
3. Description: `Utilities for CREsted package - cross-species predictions and analysis`
4. Choose Public or Private
5. **DO NOT** check "Initialize this repository with a README"
6. Click "Create repository"

### Step 4: Connect and Push

GitHub will show you commands. Use these:

```bash
# Add your GitHub repository as remote
git remote add origin https://github.com/YOUR_USERNAME/crested-utils.git

# Rename branch to main (if needed)
git branch -M main

# Push to GitHub
git push -u origin main
```

Replace `YOUR_USERNAME` with your actual GitHub username (e.g., `jjanssens`).

### Step 5: (Optional) Use GITHUB_README as main README

If you want to use the GitHub-specific README:

```bash
mv GITHUB_README.md README.md
git add README.md
git commit -m "Use GitHub README as main README"
git push
```

---

## Option 3: Fastest Method (Copy-Paste)

```bash
# Navigate to directory
cd /Users/jjanssens/Documents/Postdoc/Projects/evolution/crested_mod

# Update username (replace jjanssens with yours if different)
sed -i '' 's/YOUR_USERNAME/jjanssens/g' GITHUB_README.md

# Initialize and commit
git init
git add .
git commit -m "Initial commit: CREsted utilities"

# You'll need to create the repo on GitHub first, then:
git remote add origin https://github.com/jjanssens/crested-utils.git
git branch -M main
git push -u origin main
```

---

## What Will Be Pushed

Your repository will include:

### Core Modules (3)
- `crested_utils.py` - Full-featured utilities
- `minimal_predict.py` - Minimal version for zero overlap
- `insilico_mutagenesis_vect.py` - ISM utilities

### Documentation (11)
- `GITHUB_README.md` / `README.md` - Main README
- `GETTING_STARTED.md` - Quick start
- `QUICK_REFERENCE.md` - Cheat sheet
- `MINIMAL_README.md` - Minimal version docs
- `WHICH_VERSION.md` - Version selection guide
- `VISUAL_GUIDE.md` - Visual workflows
- `MIGRATION_GUIDE.md` - Migration guide
- `SUMMARY.md` - Overview
- `INDEX.md` - Navigation
- `FILE_STRUCTURE.md` - Repository structure

### Examples (3)
- `example_usage.py` - Full version examples
- `demo_minimal.py` - Minimal version demo
- `insilico_mutagenesis_vect_example.py` - ISM example

### Tests (1)
- `test_crested_utils.py` - Automated tests

### Legacy (4) - For reference
- `predict_chunked.py`
- `species_pred.py`
- `rowwise_corr.py`
- `resize_region.py`

### Configuration (3)
- `.gitignore`
- `LICENSE`
- `setup_github.sh`

**Total: ~25 files, all organized and ready!**

---

## After Pushing

1. **Add repository description** on GitHub
2. **Add topics/tags**: `crested`, `genomics`, `bioinformatics`, `python`, `single-cell`, `regulatory-genomics`
3. **Enable issues** if you want people to report bugs
4. **Add repository to your profile** (pin it if important)
5. **Share** the link with collaborators

---

## Updating the Repository

After the initial push, to add new changes:

```bash
# Make your changes, then:
git add .
git commit -m "Description of changes"
git push
```

---

## Troubleshooting

### "remote origin already exists"
```bash
git remote remove origin
git remote add origin https://github.com/YOUR_USERNAME/crested-utils.git
```

### Authentication issues
If prompted for password, you may need to use a Personal Access Token:
1. Go to GitHub → Settings → Developer settings → Personal access tokens
2. Generate new token with `repo` scope
3. Use token as password

Or set up SSH:
```bash
git remote set-url origin git@github.com:YOUR_USERNAME/crested-utils.git
```

### Large files warning
The repository should be small. If you have large data files:
1. Add them to `.gitignore`
2. Remove from git: `git rm --cached large_file.h5ad`
3. Commit: `git commit -m "Remove large files"`

---

## Recommended Repository Settings

After creating on GitHub:

1. **Description**: "Utilities for CREsted: cross-species regulatory predictions, cell type alignment, and in-silico mutagenesis"

2. **Topics**: `crested`, `bioinformatics`, `genomics`, `python`, `regulatory-genomics`, `single-cell`, `chromatin-accessibility`

3. **Website**: Link to CREsted docs: `https://crested.readthedocs.io/`

4. **README**: Will show `README.md` automatically

5. **License**: MIT (already included)

---

## Questions?

If you encounter issues:
1. Check git status: `git status`
2. Check remotes: `git remote -v`
3. Check branches: `git branch -a`

Or see GitHub's guide: https://docs.github.com/en/get-started/importing-your-projects-to-github

---

**Ready to push?** Follow Option 1 (automated) or Option 2 (manual) above!
