all: build/V.pdf

build/V.pdf: FORCE V.tex build/plot.pdf | build
	TEXTINPUTS="$(call translate,build:)" \
	BIBINPUTS=build: \
	max_print_line=1048576 \
	latexmk	\
		--lualatex \
		--output-directory=build	\
		--interaction=nonstopmode	\
		--halt-on-error	\
	V.tex

build/plot.pdf: V.py .txt |build

	TEXTINPUTS=$$(pwd): python V.py

build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean
	