all: build/datenübertragung.pdf

build/datenübertragung.pdf: FORCE datenuebertragung.tex | build
	TEXTINPUTS="$(call translate,build:)" \
	BIBINPUTS=build: \
	max_print_line=1048576 \
	latexmk	\
		--lualatex \
		--output-directory=build	\
		--interaction=nonstopmode	\
		--halt-on-error	\
	datenuebertragung.tex


build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: FORCE all clean