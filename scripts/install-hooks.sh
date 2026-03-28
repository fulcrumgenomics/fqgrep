#!/bin/bash
#
# Install git hooks for fqgrep development
#
# Usage: ./scripts/install-hooks.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
HOOKS_SOURCE="${SCRIPT_DIR}/hooks"
HOOKS_DEST="${REPO_ROOT}/.git/hooks"

echo "Installing git hooks for fqgrep..."

# Verify we're in the right place
if [[ ! -f "${REPO_ROOT}/Cargo.toml" ]]; then
    echo "Error: Must be run from the fqgrep repository"
    exit 1
fi

# Verify .git directory exists
if [[ ! -d "${REPO_ROOT}/.git" ]]; then
    echo "Error: .git directory not found"
    exit 1
fi

# Create hooks directory if it doesn't exist
mkdir -p "${HOOKS_DEST}"

# Install each hook
for hook in "${HOOKS_SOURCE}"/*; do
    if [[ -f "${hook}" ]]; then
        hook_name="$(basename "${hook}")"
        dest="${HOOKS_DEST}/${hook_name}"

        # Check if hook already exists
        if [[ -e "${dest}" ]]; then
            if [[ -L "${dest}" ]]; then
                echo "Updating symlink: ${hook_name}"
                rm "${dest}"
            else
                echo "Backing up existing ${hook_name} to ${hook_name}.backup"
                mv "${dest}" "${dest}.backup"
            fi
        else
            echo "Installing: ${hook_name}"
        fi

        # Create symlink (so updates to scripts/hooks are automatically picked up)
        ln -s "${hook}" "${dest}"

        # Ensure it's executable
        chmod +x "${hook}"
    fi
done

echo ""
echo "Git hooks installed successfully!"
echo ""
echo "The following hooks are now active:"
echo "  - pre-commit: Runs cargo ci-fmt and ci-lint before each commit"
echo ""
echo "To run tests in pre-commit hook, use:"
echo "  FQGREP_PRECOMMIT_TEST=1 git commit"
echo ""
echo "To bypass hooks (use sparingly):"
echo "  git commit --no-verify"
