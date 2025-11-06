# Contributing to CREsted Utilities

Thank you for your interest in contributing! This document provides guidelines for contributing to this project.

## 🎯 How to Contribute

### Reporting Bugs

If you find a bug:

1. **Check existing issues** to see if it's already reported
2. **Create a new issue** with:
   - Clear title
   - Description of the bug
   - Steps to reproduce
   - Expected vs actual behavior
   - Your environment (Python version, OS, CREsted version)
   - Minimal code example (if possible)

### Suggesting Enhancements

Have an idea for improvement?

1. **Check existing issues** for similar suggestions
2. **Create a new issue** with:
   - Clear description of the enhancement
   - Use case / motivation
   - Proposed implementation (optional)
   - Examples of how it would be used

### Pull Requests

Want to contribute code?

1. **Fork the repository**
2. **Create a feature branch**: `git checkout -b feature/your-feature-name`
3. **Make your changes**
4. **Test your changes**: Run `python test_crested_utils.py`
5. **Update documentation** if needed
6. **Commit with clear messages**: `git commit -m "Add feature: description"`
7. **Push to your fork**: `git push origin feature/your-feature-name`
8. **Open a Pull Request**

## 📝 Code Style

### Python

- Follow **PEP 8** style guide
- Use **type hints** for function parameters and returns
- Write **docstrings** for all functions (Google style)
- Keep functions **focused and modular**
- Add **comments** for complex logic

Example:
```python
def example_function(
    param1: str,
    param2: int = 10,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Brief description of what the function does.
    
    Parameters
    ----------
    param1 : str
        Description of param1
    param2 : int, default=10
        Description of param2
    verbose : bool, default=True
        Enable verbose output
        
    Returns
    -------
    DataFrame
        Description of return value
        
    Examples
    --------
    >>> result = example_function("test")
    >>> print(result)
    """
    # Implementation
    pass
```

### Documentation

- Use **Markdown** for documentation files
- Keep **line length < 100** characters when possible
- Use **code blocks** with language specification
- Include **examples** in documentation
- Update **relevant files** when adding features

## 🧪 Testing

### Before Submitting

1. **Run existing tests**:
   ```bash
   python test_crested_utils.py
   ```

2. **Test your additions**:
   - Write tests for new functions
   - Test with different inputs
   - Test edge cases

3. **Check examples still work**:
   ```bash
   python demo_minimal.py
   ```

### Writing Tests

Add tests to `test_crested_utils.py`:

```python
def test_your_new_function():
    """Test your new function."""
    print("Testing your_new_function()...")
    
    # Test basic functionality
    result = your_new_function(input_data)
    assert result is not None
    
    # Test edge cases
    result_empty = your_new_function([])
    assert len(result_empty) == 0
    
    print("  ✓ All tests passed")
```

## 📚 Documentation Updates

If your contribution adds or changes functionality:

1. **Update main README** (`GITHUB_README.md` / `README.md`)
2. **Update relevant guides** (e.g., `QUICK_REFERENCE.md`)
3. **Add examples** to `example_usage.py` if applicable
4. **Update docstrings** in the code
5. **Consider adding** visual examples to `VISUAL_GUIDE.md`

## 🔍 Code Review Process

After submitting a PR:

1. **Automated checks** will run (if configured)
2. **Maintainers review** your code
3. **Discussion** may happen in PR comments
4. **Changes requested** - address feedback
5. **Approval** - your code gets merged!

## 🌟 Good First Contributions

New to contributing? Try these:

- **Documentation improvements**: Fix typos, clarify explanations
- **Example additions**: Add new use case examples
- **Test coverage**: Add tests for existing functions
- **Performance**: Optimize existing functions
- **Bug fixes**: Fix issues from the issue tracker
- **Visualization**: Add plotting examples

## 💬 Communication

- **Be respectful** and constructive
- **Ask questions** if anything is unclear
- **Provide context** in issues and PRs
- **Be patient** - reviews take time

## 📋 Checklist for PRs

Before submitting:

- [ ] Code follows style guidelines
- [ ] All tests pass
- [ ] New code has tests
- [ ] Documentation is updated
- [ ] Examples work
- [ ] Commit messages are clear
- [ ] Branch is up to date with main

## 🙏 Thank You!

Every contribution, no matter how small, is valuable. Thank you for helping make this project better!

## 📞 Questions?

Not sure about something? Feel free to:
- Open an issue with your question
- Ask in a discussion thread
- Reach out to maintainers

---

*Happy contributing!* 🚀
