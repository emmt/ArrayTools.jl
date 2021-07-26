SOURCES = \
    make.jl \
    makedocs.sh \
    src/broadcasting.md \
    src/indexing.md \
    src/index.md \
    src/install.md \
    src/reference.md \
    src/rubberindex.md \
    src/storage.md

default: build/index.html

clean:
	rm -f *~ */*~ */*/*~

build/index.html: $(SOURCES)
	./makedoc.sh

rebuild:
	./makedoc.sh

.PHOONY: default clean rebuild
