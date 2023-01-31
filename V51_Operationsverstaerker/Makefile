all: build/main.pdf
 
# hier Python-Skripte: 
#int100.pdf: FORCE | build
#	  MATPLOTLIBRC="$$(dirname "$$(pwd)")/default/matplotlibrc" \
#	  TEXINPUTS="$$(dirname "$$(pwd)")/default:" \
#	python int.py
#umint.pfg: FORCE | build
#	  MATPLOTLIBRC="$$(dirname "$$(pwd)")/default/matplotlibrc" \
#	  TEXINPUTS="$$(dirname "$$(pwd)")/default:" \
#	python umint.py
#dif.pdf: FORCE | build
#	  MATPLOTLIBRC="$$(dirname "$$(pwd)")/default/matplotlibrc" \
#	  TEXINPUTS="$$(dirname "$$(pwd)")/default:" \
#	python dif.py

# hier weitere Abhängigkeiten für build/main.pdf deklarieren:
build/main.pdf: #int100.pdf umint.pdf dif.pdf

build/main.pdf: FORCE | build
	  TEXINPUTS=build: \
	  BIBINPUTS=build: \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=nonstopmode \
	  --halt-on-error \
	main.tex

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean
