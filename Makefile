.PHONY: format format-check

format:
	clang-format --Werror -i -style=file */*.hpp */*.cpp

format-check:
	clang-format --dry-run --Werror -style=file */*.hpp */*.cpp
