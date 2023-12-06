all: build/V53.pdf

build/V53.pdf: FORCE V53.tex V53.py | build
	TEXTINPUTS="$(call translate,build:)" \
	BIBINPUTS=build: \
	max_print_line=1048576 \
	python V53.py 
	latexmk	\
		--lualatex \
		--output-directory=build	\
		--interaction=nonstopmode	\
		--halt-on-error	\
	V53.tex
	


build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean
	