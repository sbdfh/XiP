CC = g++
CFLAGS = -fPIC $$(pkg-config --cflags --libs python3)
RFLAGS = -O3 -DNDEBUG
DFLAGS = -Wall -Wextra -pedantic -g
LDFLAGS = -shared

RTARGET = libxip.so
DTARGET = libxip_debug.so

SOURCE = xip.cpp
HEADERS = $(wildcard *.hpp)

DOXY = doxygen
DOXYCFG = doxygen.cfg

.PHONY: all
all: $(RTARGET)

$(RTARGET): $(SOURCE) $(HEADERS)
	$(CC) $(CFLAGS) $(RFLAGS) $(SOURCE) $(LDFLAGS) -o $(RTARGET)

.PHONY: debug
debug: $(DTARGET)

$(DTARGET): $(SOURCE) $(HEADERS)
	$(CC) $(CFLAGS) $(DFLAGS) $(SOURCE) $(LDFLAGS) -o $(DTARGET)

.PHONY: doxy
doxy:

	$(DOXY) $(DOXYCFG)

.PHONY: clean
clean:

	rm -f $(RTARGET) $(DTARGET)
