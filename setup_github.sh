#!/bin/bash

# Setup script for preparing GitHub repository
# Run this before pushing to GitHub

set -e  # Exit on error

echo "=========================================="
echo "Setting up crested-utils for GitHub"
echo "=========================================="

# Get GitHub username
read -p "Enter your GitHub username: " GITHUB_USER

if [ -z "$GITHUB_USER" ]; then
    echo "Error: GitHub username cannot be empty"
    exit 1
fi

echo ""
echo "GitHub username: $GITHUB_USER"
echo ""

# Replace placeholder in README
echo "Updating README with your GitHub username..."
if [ -f "GITHUB_README.md" ]; then
    sed "s/YOUR_USERNAME/$GITHUB_USER/g" GITHUB_README.md > README_temp.md
    mv README_temp.md GITHUB_README.md
    echo "✓ Updated GITHUB_README.md"
fi

# Also update the main README if it exists
if [ -f "README.md" ]; then
    sed -i.bak "s/YOUR_USERNAME/$GITHUB_USER/g" README.md
    rm README.md.bak
    echo "✓ Updated README.md"
fi

echo ""
echo "=========================================="
echo "Initializing Git repository"
echo "=========================================="

# Initialize git if not already initialized
if [ ! -d ".git" ]; then
    git init
    echo "✓ Initialized git repository"
else
    echo "✓ Git repository already initialized"
fi

# Add all files
echo ""
echo "Adding files to git..."
git add .

echo "✓ Files added"

# Create initial commit
echo ""
echo "Creating initial commit..."
git commit -m "Initial commit: CREsted utilities for cross-species analysis

- Full version with alignment and chunking (crested_utils.py)
- Minimal version for zero cell type overlap (minimal_predict.py)
- In-silico mutagenesis utilities (insilico_mutagenesis_vect.py)
- Comprehensive documentation
- Examples and tests"

echo "✓ Initial commit created"

echo ""
echo "=========================================="
echo "Next steps:"
echo "=========================================="
echo ""
echo "1. Create a new repository on GitHub:"
echo "   - Go to: https://github.com/new"
echo "   - Name: crested-utils (or your preferred name)"
echo "   - Description: Utilities for CREsted package - cross-species predictions and analysis"
echo "   - Keep it public or private as you prefer"
echo "   - DO NOT initialize with README (we already have one)"
echo ""
echo "2. Connect to your GitHub repository:"
echo "   git remote add origin https://github.com/$GITHUB_USER/crested-utils.git"
echo ""
echo "3. Push to GitHub:"
echo "   git branch -M main"
echo "   git push -u origin main"
echo ""
echo "4. (Optional) Rename GITHUB_README.md to be the main README:"
echo "   mv GITHUB_README.md README.md"
echo "   git add README.md"
echo "   git commit -m 'Use GitHub README as main README'"
echo "   git push"
echo ""
echo "=========================================="
echo "Setup complete!"
echo "=========================================="
